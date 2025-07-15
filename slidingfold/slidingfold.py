from jinja2 import Environment, FileSystemLoader
from Bio import SeqIO
from Bio.Seq import Seq
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import argparse
import RNA
import re
import os
import sys

version= '25/07/10'

# setting up jinja
env = Environment(loader=FileSystemLoader(__file__.rsplit(os.sep,1)[0]+os.sep), trim_blocks=True, lstrip_blocks=True)
template = env.get_template('slidingfold_jinja_template.html') 

# parsing command line arguments
parser = argparse.ArgumentParser(description="a program that makes heatmaps to analyze transient RdRp-trapping RNA structures")
parser.add_argument("input", help="location of your fasta file")
parser.add_argument("-fp", "--footprint", help="size of the polymerase footprint and active site location, separated by a colon, can be comma-separated list, e.g.: -fp 20:5,22:6, default 20:5", default='20:5', type=str, metavar="N")
parser.add_argument("-l", "--local", help='how many nts on each side to consider, can be comma-separated list, can specify ranges using hyphens, e.g.: -l 5,10-15, default 10', default='10', type=str, metavar='N')
parser.add_argument('-le','--lower_edge', help='set the lower edge for the heatmap scale', type=float, default=None)
parser.add_argument('-ue','--upper_edge', help='set the upper edge for the heatmap scale', type=float, default=None)
parser.add_argument('-rf', '--reading_frame', help='specify the reading frame to be used for translation', default=None, type=int, metavar='N')
parser.add_argument('-tag', '--output_tag', help='suffix to add to the generated file', type=str, default='')
parser.add_argument('-lp', '--lonely_pairs', help='allow lonely RNA pairs to form, default=False', action='store_true')
parser.add_argument('-rb', '--roadblock', help='use roadblock model, default=False', action='store_true')
parser.add_argument('-pdf', '--pdf', help='output a static pdf file instead of the interactive html, default=False', action='store_true')
parser.add_argument('-t', '--text', help='output a text file with all stems instead the interactive html, default=False', action='store_true')
parser.add_argument('-noAA', '--no_amino_acid', help='do not display the amino acid sequence in the pdf or html output, default=False', action='store_true')
args=parser.parse_args()

## setting up the arbitrary sequence and folding constraint
# this sequence will be inserted between the entering and exiting template, replacing
# the sequence in the RdRp footprint. The folding constraint prevents it from interacting
# with any other part of the sequence. This structure has its own dG which has to be
# subtracted from the overal dG given by ViennaRNA to only keep the dG of the entering
# and exiting template

stem_to_insert = 'AAGGGGGGAAAACCCCCCAA'
stem_to_insert_constraint = 'xx((((((xxxx))))))xx'
stem_to_insert_dG_correction=-13.5 # kcal/mol

## defining functions

# define default start positions and view lengths for samples (default zoom range)
# these are read from a TSV file in the same directory as this script
# the TSV file should have the following format:
# # Default values: <default_start_position> <default_view_length>
# # seq_id <start_position> <view_length>
# seq1 10 200
# seq2 15 250
# seq3 5 150
# The first line is a comment and defines the default values for start position and view length.
# The subsequent lines define the start position and view length for each sequence ID.
# Sequence IDs should correspond to the IDs in the FASTA file.
# The default values will be used if a sequence ID is not found in the TSV file.

# Load sample positions from TSV
SAMPLE_START_POSITIONS = {}
DEFAULT_START_POSITION = 1
DEFAULT_VIEW_LENGTH = 200

script_dir = os.path.dirname(os.path.abspath(__file__))
tsv_path = os.path.join(script_dir, 'sample_start_positions.tsv')
if os.path.exists(tsv_path):
    with open(tsv_path, 'r') as f:
        for line in f:
            if line.startswith('#'):  # Skip comments
                if line.startswith('# Default values'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        DEFAULT_START_POSITION = int(parts[1])
                        DEFAULT_VIEW_LENGTH = int(parts[2])
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                seq_id = parts[0]
                SAMPLE_START_POSITIONS[seq_id] = {
                    'start': int(parts[1]),
                    'length': int(parts[2])
                }

def find_default_range(record_id):
    """Return default zoom range based on sample configuration"""
    if record_id in SAMPLE_START_POSITIONS:
        start = SAMPLE_START_POSITIONS[record_id]['start']
        length = SAMPLE_START_POSITIONS[record_id]['length']
    else:
        start = DEFAULT_START_POSITION
        length = DEFAULT_VIEW_LENGTH
        
    return [start - 0.5, start + length - 0.5]

# peel apart comma-separated variables
def peel_apart_vars(varstring):
    split=varstring.split(',')
    peeled=tuple()
    for item in split:
        if '-' in item: # require further splitting
            peeled+=tuple(range(int(item.split('-')[0]), int(item.split('-')[1])+1))
        else:
            peeled+=(int(item),)
    return peeled

# pad variables until set length
def pad_vars(var, lenmax):
    if len(var)<lenmax: # means they are length 1
        var=tuple(var*lenmax)
    return var

# determine bp type for each position
def get_basepairs(seq, struc):
    gc, gt, at = [],[],[]

    # we do this iteratively by identifying inner basepairs, i.e., a set of ( and ) separated
    # by only ., we then identify the basepair formed by referring to the corresponding 
    # positions int he sequence, and remove this basepair from the analysis by replacing both
    # parentheses by . before iterating

    bubble_regex=re.compile("\(\.*\)")

    while bubble_regex.search(struc): # while there are () pairs remaining
        bubble=bubble_regex.search(struc) # find an inner pair
        (leftbracket, rightbracket) = (bubble.span()[0], bubble.span()[1]-1) # get positions
        # check bp sequence and keep track of positions
        if seq[leftbracket].upper() == "A" or seq[rightbracket].upper() == "A":
            at.extend([leftbracket,rightbracket])
        elif seq[leftbracket].upper() == "C" or seq[rightbracket].upper() == "C":
            gc.extend([leftbracket,rightbracket])
        else:
            gt.extend([leftbracket,rightbracket])
        # replace identified () pair by .. before iterating
        # we know only .s separate the parentheses, so we can replace the
        # whole stretch at once
        struc = struc[:leftbracket]+"."* (rightbracket - leftbracket+1)+struc[rightbracket+1:]

    # no () pairs left, time to return our results
    return (gc, gt, at)

## set up dictionaries containing values and colors for the sequence heatmaps
# colors are taken from the excellent AliView (https://ormbunkar.se/aliview/)
nts_values ={
            "A" : 0,
            "C" : 1,
            "U" : 2,
            "T" : 2,
            "G" : 3,
            "R" : 4,
            "Y" : 4,
            "N" : 4,
            } 
nts_colorscale=(
                (0,"green"),
                (0.25,'blue'),
                (0.50,'red'),
                (0.75,"gold"),
                (1, 'white'),
                )

AA_values = {
            "F" : 0,
            "L" : 1,
            "I" : 2,
            "M" : 3,
            "V" : 4,
            "S" : 5,
            "P" : 6,
            "T" : 7,
            "A" : 8,
            "Y" : 9,
            "H" : 10,
            "Q" : 11,
            "N" : 12,
            "K" : 13,
            "D" : 14,
            "E" : 15,
            "C" : 16,
            "W" : 17,
            "R" : 18,
            "G" : 19,
            "*" : 20,
            '-' : 20,
            'X'  :20,
            }

aa_colorscale=(
            (0.0, 'rgb(25,127,229)'),
            (0.05, 'rgb(119,165,214)'),
            (0.10, 'rgb(79,160,242)'),
            (0.15, 'rgb(15,84,155)'),
            (0.20, 'rgb(5,124,249)'),
            (0.25, 'rgb(2,150,2)'),
            (0.30, 'rgb(204,204,0)'),
            (0.35, 'rgb(68,201,68)'),
            (0.40, 'rgb(38,109,183)'),
            (0.45, 'rgb(20,198,198)'),
            (0.50, 'rgb(25,178,178)'),
            (0.55, 'rgb(91,237,91)'),
            (0.60, 'rgb(25,204,25)'),
            (0.65, 'rgb(229,51,25)'),
            (0.70, 'rgb(204,76,204)'),
            (0.75, 'rgb(153,63,150)'),
            (0.80, 'rgb(229,127,127)'),
            (0.85, 'rgb(2,84,168)'),
            (0.90, 'rgb(244,68,43)'),
            (0.95, 'rgb(229,153,76)'),
            (1, 'rgb(255,255,255)'),
            )

## running 

# get the folding variables from the provided arguments
footprints=tuple(int(arg.split(':')[0]) for arg in args.footprint.split(','))
activelocations=tuple(int(arg.split(':')[1]) for arg in args.footprint.split(','))
windows=peel_apart_vars(args.local)

# check if the variables have the same length
# if they don't we need to fill in the gaps to get a number of parameter sets
# that all have values for all three parameters
if not len(footprints)==len(activelocations)==len(windows):
    print('different numbers of parameters, filling in the gaps')
    fold_vars=tuple()
    # active locations and footprints are always provided simultaneously, so we can iterate through both
    for activelocation, footprint in enumerate(footprints):
        for window in windows:
            fold_vars+=((footprint, activelocations[activelocation], window),)
else:
    # zip up the vars
    fold_vars=tuple(zip(footprints, activelocations, windows))

# setting up overall variables
figs={}
orfs={}
first_id = ''

# iterate through the sequences found in the fasta file
for record in SeqIO.parse(args.input, 'fasta'):
    # setting up per-sequence variables
    dGmin=0
    results={}

    if first_id == '':
        first_id = record.id

    seqs={'forward':str(record.seq).upper(), 'reverse':str(record.seq.reverse_complement()).upper()}

    for fv_nb, fold_var in enumerate(fold_vars): # iterate through fp, al, local parameter triplets
        # setting up per-parameter variables
        results[fv_nb]={}
        footprint, activelocation, local = fold_var

        # read through sequence in forward, then in reverse, save results
        for orientation in ('forward','reverse'):
            # setting up per-orientation variables
            seq = seqs[orientation]
            results[fv_nb][orientation]= {'dG':[],'struc':[]}

            print(f'\rscanning {orientation} of {record.id} with fp={footprint}, al={activelocation}, l={local}',end='')

            # iterate through sequence
            for i in range(1,len(seq)+1):

                # we can only fold if there is enough sequence in 5' and 3' to have the full window present
                if local+activelocation <= i < len(seq)-(local+footprint-activelocation): # means we can fold

                    # we chop up the parts of the sequence relevant for folding
                    upstream = seq[i-(local+activelocation):i-activelocation]
                    downstream = seq[i+(footprint-activelocation):i+(footprint-activelocation)+local] if args.roadblock == False else ''
                    activesite = seq[i-1:i+9]

                    # insert the arbitrary sequence in between up and downstream
                    foldseq = upstream+stem_to_insert+downstream
                    # add the folding constraint, if a nucleotide is not covered, it is free to do whatever
                    constraint = '.'*local+stem_to_insert_constraint

                    # set up the ViennaRNA folding model and compound
                    md = RNA.md() # create model
                    md.noLP = 0 if args.lonely_pairs == True else 1 # (dis)allow isolated basepairs
                    fc = RNA.fold_compound(foldseq,md) # combine sequence and model
                    fc.hc_add_from_db(constraint) # add constraint to model

                    # fold using our folding compound
                    (struc, dG) = fc.mfe()

                    # reformat the structure according to nature of basepairs
                    # A-T basepairs get left as is here
                    (gc, gt, at)=get_basepairs(foldseq,struc)
                    struc_form = ''
                    for p, c in enumerate(struc):
                        if p in gc:
                            struc_form += '{' if c == '(' else '}'
                        elif p in gt:
                            struc_form += '<' if c == '(' else '>'
                        else:
                            struc_form += c
                    
                    # delete the arbitrary sequence and replace by RdRp footprint with activ site duplex
                    struc = struc_form[:local]+f'[{"-"*(activelocation-2)}{activesite}{"-"*(footprint-activelocation-10)}]'+ (struc_form[-local:] if args.roadblock == False else '')

                    # save structure and dG
                    results[fv_nb][orientation]['struc'].append(struc)
                    results[fv_nb][orientation]['dG'].append(dG-(stem_to_insert_dG_correction))
                    
                    # update minimum deltaG if necessary
                    dGmin=dG-(stem_to_insert_dG_correction) if dG-(stem_to_insert_dG_correction)<dGmin else dGmin

                else: # too close to the edges to fold, just pad the lists
                    results[fv_nb][orientation]['struc'].append('too close to edge')
                    results[fv_nb][orientation]['dG'].append('no dG calculated')
    
    print('\n', end='')

    if args.text==True:
    # if we just want a text file output, we have everything we need here
        with open(f'{args.input.rsplit(".",1)[0]}_{record.id}_stems{args.output_tag}{"_roadblock" if args.roadblock==True else ""}.tsv', 'w') as f:
            # write header
            f.write(f'# stems found using slidingfold version {version} for sequence {record.id} from {args.input}\n')
            f.write(f'# arbitrary sequence {stem_to_insert} with constraint {stem_to_insert_constraint}{" and roadblock model" if args.roadblock==True else ""}, lonely pairs were {"allowed" if args.lonely_pairs==True else "forbidden"}\n# folding parameter sets (footprint, active site location, window size) were:\n')
            for fv_nb, fold_var in enumerate(fold_vars):
                f.write(f'# {fv_nb+1}\t'+'\t'.join(map(str,fold_var))+'\n')
            f.write('\t'.join(('# folding_parameter_set','position','dG_forward','structure_forward', 'dG_reverse','structure_reverse'))+'\n')

            for fv_nb, fold_var in enumerate(fold_vars): 
                for pos in range(len(results[fv_nb][orientation]['dG'])):
                    f.write(f'{fv_nb+1}\t{pos+1:>4}\t'+
                            f'{results[fv_nb]["forward"]["dG"][pos]:<20}\t{results[fv_nb]["forward"]["struc"][pos].ljust(fold_var[0]+fold_var[2]*2)}\t'+
                            f'{results[fv_nb]["reverse"]["dG"][-(pos+1)]:<20}\t{results[fv_nb]["reverse"]["struc"][-(pos+1)].ljust(fold_var[0]+fold_var[2]*2)}\n')
                    
        print(f'wrote stems to {args.input.rsplit(".",1)[0]}_{record.id}_stems{args.output_tag}{"_roadblock" if args.roadblock==True else ""}.tsv')

    else:
    # if we want the pdf or html, we actually have to make the graphs

        if args.no_amino_acid == False:
            ## figure out the amino acid stuff
            # if the reading frame has not been specified, try to find the best one, ie the one with the longest orf
            if not args.reading_frame:
                transl=[max(list(map(len,str(Seq(seqs['forward'][i:(len(seqs['forward'][i:])//3*3)+i]).translate()).split('*')))) for i in range(3)] # translate, split, get length, keep max
                frame=transl.index(max(transl)) # get index of longest ORF
            else:
                frame=args.reading_frame-1 # converting to 0 index

            # save the reading frame for later    
            orfs[record.id]=frame

            # get the corrresponding protein sequence
            prot=str(Seq(seqs['forward'][frame:len(seqs['forward'][frame:])//3*3+frame]).translate())
            
            # pad the protein sequence so it fits with the nucleotide sequence
            # for the heatmap values we simply put the 1-letter AA code on every position of the codon
            prot_padded = '-'*frame+''.join([c+c+c for c in prot])+'-'*(len(seqs['forward'])-len(prot)*3-frame)
            # for the heatmap text we center the 1-letter code on position 2 of the codon
            prot_padded_text = ' '*frame+''.join([f' {c} ' for c in prot])+' '*(len(seqs['forward'])-len(prot)*3-frame)

        # find the range to zoom on by default
        default_range = find_default_range(record.id)

        # now make figure

        ## first figure out how many rows we have and what the row heights should be
        ## there are always 3 rows of invariant height in the middle for the sequences
        row_heights=[1.0/len(fold_vars) for i in range(len(fold_vars))]+([0.25]*3 if args.no_amino_acid==False else [0.25]*2)+[1.0/len(fold_vars) for i in range(len(fold_vars))]

        fig=make_subplots(rows=len(row_heights), cols=1, row_heights=row_heights, vertical_spacing=0, shared_xaxes=True)

        ## plotting the cRNA stems
        ## i want to go inside-out: closest to seq should be first result-dict
        for fv_nb in range(len(fold_vars)):
            footprint, activelocation, local = fold_vars[fv_nb]
            fig.add_trace(go.Heatmap(z=[results[fv_nb]['forward']['dG']], 
                                    x0=1,dx=1,
                                    colorscale='purples_r',
                                    customdata=[[(results[fv_nb]['forward']['struc'][i].replace('T','U'),
                                                f"∆G= {results[fv_nb]['forward']['dG'][i]:.2f} kcal/mol" if isinstance(results[fv_nb]['forward']['dG'][i],float) else "no ∆G calculated"
                                                ,i+1,
                                                seqs['forward'][i+1-(local+activelocation):i+1+(footprint-activelocation)+local].replace('T','U'),
                                                footprint, activelocation, local
                                                ) for i in range(len(results[fv_nb]['forward']['dG']))]],
                                    hovertemplate='footprint: %{customdata[4]} nts, active site in position %{customdata[5]}<br>%{customdata[6]} nts on each side for folding<br>position %{customdata[2]}\t%{customdata[1]}<br>%{customdata[3]}<br>%{customdata[0]}<extra></extra>',
                                    showscale=False,
                                    zmax= args.upper_edge if args.upper_edge else 0,
                                    zmin= args.lower_edge if args.lower_edge else dGmin,
                                    ), 
                            row=len(fold_vars)-fv_nb, col=1,
                            )
        ## plotting the cRNA sequence
        fig.add_trace(go.Heatmap(z=[[nts_values[c] for c in seqs['forward']]],
                                x0=1,dx=1,
                                colorscale=nts_colorscale,
                                text=[[c for c in seqs['forward'].replace('T','U')]],
                                texttemplate='%{text}',
                                customdata=[[[v+1, len(seqs['forward'])-v] for v in range(len(seq))]],
                                hovertemplate='cRNA position: %{customdata[0]}<br>vRNA position: %{customdata[1]}<extra></extra>',
                                showscale=False,
                                zmax=4,
                                zmin=0
                                ),
                    row=len(fold_vars)+1, col=1)
        
        if args.no_amino_acid == False:
            # plotting the AA sequence
            fig.add_trace(go.Heatmap(z=[[AA_values[c] for c in prot_padded]],
                                    x0=1,dx=1,
                                    colorscale=aa_colorscale,
                                    text=[[c for c in prot_padded_text]],
                                    texttemplate='%{text}',
                                    hoverinfo='skip',
                                    showscale=False,
                                    zmax=20,
                                    zmin=0
                                    ),
                        row=len(fold_vars)+2, col=1)

        ## plotting the vRNA sequence
        fig.add_trace(go.Heatmap(z=[[nts_values[c] for c in seqs['reverse'][::-1]]],
                                x0=1,dx=1,
                                colorscale=nts_colorscale,
                                text=[[c for c in seqs['reverse'][::-1].replace('T','U')]],
                                texttemplate='%{text}',
                                customdata=[[[v+1, len(seqs['forward'])-v] for v in range(len(seq))]],
                                hovertemplate='cRNA position: %{customdata[0]}<br>vRNA position: %{customdata[1]}<extra></extra>',
                                showscale=False,
                                zmax=4,
                                zmin=0
                                ),
                    row=len(fold_vars)+3 if args.no_amino_acid == False else len(fold_vars)+2, col=1)

        ## plotting the vRNA stems
        ## i want to go inside-out: closest to seq should be first result-dict
        for fv_nb in range(len(fold_vars)):
            footprint, activelocation, local = fold_vars[fv_nb]
            fig.add_trace(go.Heatmap(z=[results[fv_nb]['reverse']['dG'][::-1]], 
                                    x0=1,dx=1,
                                    colorscale='purples_r',
                                    customdata=[[(results[fv_nb]['reverse']['struc'][i].replace('T','U'), 
                                                f"∆G= {results[fv_nb]['reverse']['dG'][i]:.2f} kcal/mol" if isinstance(results[fv_nb]['reverse']['dG'][i],float) else "no ∆G calculated",
                                                i+1,
                                                seqs['reverse'][i+1-(local+activelocation):i+1+(footprint-activelocation)+local].replace('T','U'),
                                                footprint, activelocation, local
                                                ) for i in range(len(results[fv_nb]['reverse']['dG'])-1,-1,-1)]],
                                    hovertemplate='footprint: %{customdata[4]} nts, active site in position %{customdata[5]}<br>%{customdata[6]} nts on each side for folding<br>position %{customdata[2]}\t%{customdata[1]}<br>%{customdata[3]}<br>%{customdata[0]}<extra></extra>',
                                    showscale=False,
                                    zmax= args.upper_edge if args.upper_edge else 0,
                                    zmin= args.lower_edge if args.lower_edge else dGmin,
                                    ), 
                            row=len(fold_vars)+4+fv_nb if args.no_amino_acid == False else len(fold_vars)+3+fv_nb, col=1, # 4 is 3 sequence tracks +1 for 1-indexing
                            )

        # customize hover text format and figure background
        fig.update_layout(
                        hoverlabel=dict(
                            bgcolor="white",
                            font_family="Courier"),
                        template='simple_white',
                        xaxis_showticklabels=True,
                        title=f"<b>{record.id}</b>",
                        margin=dict(t=60, b=10)
                        )
        
        # move x-axis for first row to the top side of the graph
        fig.update_xaxes(side='top',
                        row=1,col=1)

        # hide all y axes:
        for row in range(len(row_heights)):
            fig.update_yaxes(visible=False,
                            row=row+1,col=1)

        # hide all x axes except first and last
        for row in range(len(row_heights)-2):
            fig.update_xaxes(visible=False,
                            row=row+2,col=1)
    
        if args.pdf==False:
            # if we have a particular range to zoom to, do so
            if default_range:
                fig.update_xaxes(range=default_range)

            # just save the finished figure for now
            figs[record.id]=fig.to_html(full_html=False, include_plotlyjs='cdn' if record.id == first_id else False)
        else:
            # pdfs just get written to disk
            fig.update_layout(width=7.5*len(seq),
                              margin=dict(t=10,r=10,b=10,l=10))
            fig.write_image(f"{args.input.rsplit('.',1)[0]}-{record.id}.pdf")
            print(f'wrote pdf file to {args.input.rsplit(".",1)[0]}-{record.id}.pdf')

if args.pdf==False and args.text==False:
    ## time to actually write the interactive html file

    ## first gather all the info we need

    # html table of folding variables
    fold_var_tablerows=[f'<tr><td>{fv_nb+1}</td><td>{fold_var[0]}</td><td>{fold_var[1]}</td><td>{fold_var[2]}</td></tr>' for fv_nb, fold_var in enumerate(fold_vars)]

    # grabbing variables that need to be writted on home page
    params=[]
    params.append(f'<b>slidingfold output for {args.input}</b>')
    params.append(f'<sub>made using version {version}</sub>')
    params.append(f'<table><tr><th>parameter set</th><th>footprint</th><th>active site location</th><th>local size</th></tr>{"".join(fold_var_tablerows)}</table>')
    params.append(f'lonely pairs were {"allowed" if args.lonely_pairs == True else "forbidden"}, inserted stem and constraint:<br><br><div class="sequence">&emsp;{stem_to_insert}<br>&emsp;{stem_to_insert_constraint}</div>')
    params.append(f'upper heatmap limit: {args.upper_edge if args.upper_edge else "auto"}, lower heatmap limit: {args.lower_edge if args.lower_edge else "auto"}')
    if args.reading_frame:
        params.append(f'reading frame was set to frame {args.readingframe}')
    else:
        params.append(f'automatically detected reading frames were:<br>&emsp;'+'<br>&emsp;'.join([f'{k}: frame {v+1}' for k, v in orfs.items()]))

    ## make the html via jinja    
    with open(f'{args.input.rsplit(".",1)[0]}_slidingmaps{args.output_tag}{"_roadblock" if args.roadblock==True else ""}.html', 'w') as fj:
        #calling jinja and passing along all our variables
        fj.write(template.render(figs=figs, params=params))
        print(f'rendered html file at {args.input.rsplit(".",1)[0]}_slidingmaps{args.output_tag}{"_roadblock" if args.roadblock==True else ""}.html')