# A script to parse SAM files containing alignments of consensuses to the reference,
# identify the corresponding UMIs from the raw reads to filter PCR duplicates,
# and realign indels to the most stable predicted RdRp-trapping RNA structure

import gzip
from collections import defaultdict, Counter
import re
import tempfile
import shutil
import os

# requires Biopython, ViennaRNA, edlib
from Bio import SeqIO
import RNA
import edlib

### Set-up
# toggle to identify UMIs and filter out PCR duplicates
filter_umis=True

# PATHS PROVIDED DOWN BELOW ARE MEANT TO BE USED WITH THE DEMO DATA
# YOU WILL NEED TO ADAPT THEM TO YOUR OWN FILE SYSTEM

## outfile path
# the outfile path is a prefix/suffix tuple that will be joined with either
# 'insertions' or 'deletions' to create 2 separate files, one per indel type
# e.g.: 'path/to/folder/' and '_realigned.tsv' would create the output files
# path/to/folder/insertions_realigned.tsv and 
# path/to/folder/deletions_realigend.tsv
out_prefix=('demo_output'+os.sep)
out_suffix=('_demo.tsv')

## defining samples
# samples are read in from a file with one tab-separated line per sample, containing the corresponding 
# NGS runs and barcodes (separated with :), e.g.: RETR	133:31	144:59
# file can include comment lines starting with #, these will be skipped
sample_file_path=os.path.join('txt_data', 'runs.tsv')

## defining references
# paths to the reference for each sample are read in from a file with one tab-separated line per sample
# the path to the fasta file has to be provided, it will then be read in.
# these names must match up with those in the samples file
# here we just keep the reference as string
reference_file=os.path.join('txt_data', 'refs.tsv')

# primers for each sample are read in from a file with one tab-separated line per sample
# cotaining the sample name followed by a tab and then the primer sequence
# these names must match up with those in the samples file
primer_file=os.path.join('txt_data','primers.tsv')

# a quick function that returns the sample-specific paths to different files
# everyone has their own way of organizing their files, but the script needs
# a way to get all file paths just form the run and barcode of each sample
# you will have to tweak this to make it fit your own datastructure
def get_paths(run, bc):
    paths=dict()

    # input files
    # first file we need is the gzipped sam file containing alignments of consensuses
    paths['sam']=os.path.join('ngs_data', run, bc, f'{run}-{bc}.sam.gz')
    # second we need the file containing the gzipped raw reads to extract UMIs from
    paths['reads']=os.path.join('ngs_data', run, bc, f'{run}-{bc}.fastq.gz')
    # the file that contains eventual blacklisted read IDs
    paths['blacklist']=os.path.join('ngs_data', run, bc, f'{run}-{bc}-blacklist.txt')

    # ouput files
    # finally an output file to write the coverage to
    paths['coverage']=os.path.join('demo_output', f'{run}-{bc}-coverage.tsv')

    return(paths)

# for fully disrupted sequences, we want to align indels the same way as the
# correponding wild type (makes it easier to compare them). This dict should
# contains disrupted:wt pairs with the sample names read it from the sample file
destruct={ 
        # for example:
        # 'RKKRdestab':'RKKR',
        '':'',
        }

outfile=(out_prefix, out_suffix)

## reading in the common input files

samples={}
with open(sample_file_path,'r') as f:
    for line in f:
        l=line.strip().split('\t')
        if not l[0].startswith('#'): 
            samples[l[0]]=tuple((tuple(run.split(':')) for run in l[1:]))

with open(reference_file, 'r') as f:
    refs={line.split('\t')[0]:{'fwd':str(SeqIO.read(line.strip().split('\t')[1], 'fasta').seq).upper(),
                               'rev':str(SeqIO.read(line.strip().split('\t')[1], 'fasta').seq.reverse_complement()).upper(),}
                               for line in f}
    
with open(primer_file, 'r') as f:
    primers={line.split('\t')[0]:line.strip().split('\t')[1] for line in f}

### FUNCTIONS

# roll an indel sequence alongside the reference to find all ways of aligning this indel. Rolling means taking the 5' nts and putting it in 3' and vice-versa
# i.e. aaa----gaaaaaaaa -> aaagaaaaaaaa -> aaagaaaaaaa -> aaagaaaaaaaa -> aaagaaaaaaaa
#      aaagaaagaaaaaaaa       gaaa           agaa          aaga           aaag
def roll_indel(indel, ref, startpos):
    indel, ref = indel.upper(), ref.upper()
    # moving left
    i = startpos
    indel_rolling, ref_rolling = indel, ref[i:i+len(indel)] # intialize at start position
    if indel_rolling != ref_rolling: # indel is not a clean duplication
        # i.e. aaaga---aaaaa -> aaagaaaaaa
        #      aaagagaaaaaaa         gaa
        return None
    while indel_rolling == ref_rolling:
        i+=1 
        ref_rolling = ref[i:i+len(indel)] # shifting by 1
        indel_rolling = indel_rolling [1:] + indel_rolling[:1] # moving 5' nts to 3'
    edge_right = i-1 # once we exit the loop we went one position too far to the right, the indel cannot align to the ref anymore
    # moving right
    i = startpos
    indel_rolling, ref_rolling = indel, ref[i:i+len(indel)]# intialize at start position
    while indel_rolling == ref_rolling:
        i += -1 
        ref_rolling = ref[i:i+len(indel)] # shifting by 1
        indel_rolling = indel_rolling [-1:] + indel_rolling[:-1] # moving  3' nts to 5'
    edge_left = i+1 # once we exit the loop we went one position too far to the right, the indel cannot align to the ref anymore
    if edge_left -2 == edge_right: # can't roll
        return (startpos, startpos)
    else:
        return (edge_left, edge_right)

# RNA folding function, based on slidingfold
def get_fold(seq, i):
    # parameters for the folding
    stem_to_insert = 'AAGGGGGGAAAACCCCCCAA'
    stem_to_insert_constraint = 'xx((((((xxxx))))))xx'
    local = 10 # length of nts to take on each side
    activelocation = 5 # location of activesite
    footprint = 20 # pol footprint size
    forbid_lonely_pairs = 1 # allow lonely pairs with 0, forbid with 1
    dG_correction = 13.50 # correction factor to add to dG (stability of the inserted stem)

    # if position is too close to the edges, we can't fold and just return an absurdly high dG
    if not local+activelocation <= i < len(seq)-(local+footprint-activelocation):
        return 9999
    
    upstream = seq[i-(local+activelocation):i-activelocation]
    downstream = seq[i+(footprint-activelocation):i+(footprint-activelocation)+local]
    foldseq = upstream+stem_to_insert+downstream # add the stem sequence into the polymerase footprint
    constraint = '.'*local+stem_to_insert_constraint # build a constraint sequence to force the inserted stem to fold, only 5' is required to be explicitly not folded
    md = RNA.md() # create folding model
    md.noLP = forbid_lonely_pairs # forbid lonely pairs
    fc = RNA.fold_compound(foldseq,md) # create folding compound
    fc.hc_add_from_db(constraint) # add our constraint
    (struc, dG) = fc.mfe() # do the folding
    return dG+dG_correction

# uses alignment via edlib to try and find a umi
def find_umi(read, primer, max_mm=0.2, max_lendiff=1):
    primer=primer.upper()
    umilen=primer.count('N')
    # split up part of the primer
    primer5, *_, primer3 = primer.split('N')

    # due to a quirk in the CirSeq alignment process we can't take
    # the orientation in the final sam file, we have to test both
    # we only need to keep the 5' most part of the sequence, since 
    # that's where the UMI should be.
    seqs=(str(read)[:len(primer)+10],str(read.reverse_complement())[:len(primer)+10])
    
    for seq in seqs:
        # we start by simply looking for perfect matches and checking they are at the right distance
        if umilen-max_lendiff <= seq.find(primer3)-(seq.find(primer5)+len(primer5)) <= umilen+max_lendiff:
            umiseq=seq[seq.find(primer5)+len(primer5):seq.find(primer3)]
            return umiseq
    
    # if that simple approach did not work, we try to align

    # add additional equalities for common 2-nts degenerate bases
    additional_equalities=[('R','A'), ('R','G'), ('Y','C'), ('Y','T')]
    
    for seq in seqs:
        # align 5' part using SHW mode, ie 3' gaps are not penalized
        aln5=edlib.align(primer5, seq, mode='SHW', additionalEqualities=additional_equalities, task='locations')
        # align 3' part using HW mode, ie 5' and 3' gaps are not penalized
        aln3=edlib.align(primer3, seq, mode='HW', additionalEqualities=additional_equalities, task='locations')

        # look if the 3' part matches well
        # i don't look at the 5' part since that region often has dels
        if aln3['editDistance']/len(primer3)<=max_mm:
            # get the end of the primer5 alignment and start of primer3 alignment
            umirange=(aln5['locations'][0][1]+1, aln3['locations'][0][0])

            # check that teh distance between both parts is OK
            if umilen-max_lendiff <= umirange[1]-umirange[0] <= umilen+max_lendiff:
                umiseq=seq[umirange[0]:umirange[1]]
                return umiseq
            
            # if its not ok, something went wrong we do one last-ditch attempt
            # we take the 10 nts in 5' of the 3' primer part, make sure there's
            # at least 10 nts in this region and check if that sequence is preceded
            # with the 3 last nts of the 5' primer part
            # this is somewhat stringent, but tolerates 5' truncations
            umiseq=seq[umirange[1]-10:umirange[1]]
            if len(umiseq)==10 and seq[umirange[1]-13:umirange[1]-10]==primer5[-3:]:
                return umiseq

    # if we made it here, we didn't find a UMI :(
    return 'XXXXXXXXXX'

# turns cigar strings into long-form cigar strings
# ie 3M1D3I would turn into MMMDIII
def unpack_cigar(cigar_str):
    if cigar_str== '*':
        return ''
    cigar_long = ''
    for item in re.findall('\d*\D',cigar_str):
        cigar_long += int(item[:-1])*item[-1:]
    return cigar_long

# make an index from a gzipped file via temp unzipping
# Biopython only supports index on unzipped files
def index_gzip(file):
    # open temp file
    with tempfile.NamedTemporaryFile(mode='w+t') as temp_zip:
        # copy the unzipped fastq file to the temp file
        shutil.copyfileobj(gzip.open(file, 'rt'), temp_zip)
        # return to top of temp file
        temp_zip.seek(0)
        # index unzipped reads in temp file
        reads=SeqIO.index(temp_zip.name, 'fastq')
    return reads

# open outfile and iterate through all our samples
with open('insertions'.join(outfile), 'w') as fins, open('deletions'.join(outfile), 'w') as fdel:
    for sample, replicates in samples.items():

        reflen=len(refs[sample.strip('_np')]['fwd'])

        # get stems for this sample
        all_dGs_f=tuple(get_fold(refs[sample.strip('_np')]['fwd'], i+1) for i in range(reflen))
        # for reverse stems, we iterate backwards to get them in forward order
        all_dGs_r=tuple(get_fold(refs[sample.strip('_np')]['rev'], reflen-1-i+1) for i in range(reflen))

        if sample.strip('_np') in destruct.keys():
            # in this case we also need stems of the structured reference
            all_dGs_f_struc=tuple(get_fold(refs[destruct[sample.strip('_np')]]['fwd'], i+1) for i in range(reflen))
            # for reverse stems, we iterate backwards to get them in forward order
            all_dGs_r_struc=tuple(get_fold(refs[destruct[sample.strip('_np')]]['rev'], reflen-1-i+1) for i in range(reflen))

        for rep_num, replicate in enumerate(replicates):
            print(f'\rscanning {sample} replicate {rep_num+1:<50}', end='')

            paths=get_paths(replicate[0],replicate[1])

            # read in the blacklist
            with open(paths['blacklist'], 'r') as f:
                blacklist=set(line[1:].strip() for line in f)

            # this set will hold previously encountered combinations of umi+start of aln+cigar
            # if a read has an already obseved combination, it will be skipped
            umi_start_cigars=set()

            # this will be used to calculate the position-by-position coverage
            coverage=[0]*reflen

            # this will hold the found indels
            # defaultdict will be by length of indel, counter will have pos+seq
            insertions=defaultdict(Counter)
            deletions=defaultdict(Counter)

            # patterns used to find indels via regexing the CIGAR string
            ins_pattern=re.compile('II*')
            del_pattern=re.compile('DD*')
            # insertions closely followed closely by deletions/the opposite
            # tend to reflect bowtie trying to max score on mismatched
            # tempaltes/references
            conta_pattern=re.compile('I.{0,5}D|D.{0,5}I')

            # if we filter the umis, we will need access to the reads,
            # so we index them here
            if filter_umis==True:
                reads=index_gzip(paths['reads'])

            # going through the sam file a first time to get a count of reads
            with gzip.open(paths['sam'], 'rt') as f:
                nreads=sum(1 for line in f)

            # going through it again, this time we look at data
            with gzip.open(paths['sam'], 'rt') as f:

                nread=0
                for line in f:
                    nread+=1
                    if nread % 10000 == 0:
                        print(f'\rscanning {sample} replicate {rep_num+1} read {nread:>6}/{nreads} ({nread/nreads*100:.1f}%)', end='')
                    # split lines up into vars
                    rid, flag, refid, start, mapq, cigar, a, b, c, alnseq, *_ = line.split('\t')

                    # skip blacklisted reads
                    if rid in blacklist:
                        continue

                    longcigar=unpack_cigar(cigar)

                    # when filtering UMIs, we consider 2 things:
                    # 1. the umi sequence
                    # 2. the start, length, and indels of the alignment
                    # since we have random fragmentation and then consensuses, this gives us 2 layers
                    # of deduplication. Here we simply check if this combination of the 2 factors has
                    # been seen before and if so we skip the read. This seems arbitrary, but since we 
                    # are not interested in substitutions, it does not matter which representative we
                    # keep, they are all by definition identical for our purposes (same start+indels)
                    if filter_umis==True:
                        umi=find_umi(reads[rid].seq, primers[sample])

                        # if we don't find a UMI, we just skip this one
                        if umi=='XXXXXXXXXX':
                            continue

                        # else we check if this combo was seen previously
                        if umi+start+cigar in umi_start_cigars:
                            # this combo was seen before, skip the line fully
                            continue
                        
                        else:
                            # never seen before, add to list and keep going
                            umi_start_cigars.add(umi+start+cigar)
                    

                    if not 'I' in longcigar and not 'D' in longcigar:
                        # for each non soft-clipped position, increase coverage by 1, adjust to 0-indexing
                        for pos in range(len(longcigar.strip('S'))):
                            coverage[int(start)-1+pos]+=1

                    else:

                        # filter out contaminants
                        # reads usually have 1 or 2 indel operations, having many more is usually sign of contamination
                        # by a closely related template. To maximize score in low similarity regions, Bowtie tends to 
                        # introduce additional indels, often a 1IxM1D which is a bit of a giveaway. 
                        # To be sure, we skip all reads that have more than 4 M operations, 4 M is 3 indels e.g.: xMxDxMxDxMxIxM
                        # splitting on M will give us a list of len n+1 where n is the number of Ms and we skip reads with 
                        # and insertion and a deletion separated by less than 4 nucleotides
                        if len(cigar.split('M'))>5:
                            continue
                        
                        if 'D' in longcigar and 'I' in longcigar:
                            if conta_pattern.search(longcigar):
                                continue

                        # we have indels in the read, we need to 
                        # 1) update the coverage
                        # 2) extract the indel information

                        # for coverage we need to correct for I characters in teh cigar, these do not consume reference
                        # a D character is considered covered here, while the nts itself is not present, the consensus does
                        # cover this position and thus should be used to calculate /1000 consensus read frequencies imo
                        for pos in range(len(longcigar.replace('I','').strip('S'))):
                            coverage[int(start)-1+pos]+=1
                    
                        for insertion in ins_pattern.finditer(longcigar):
                            # first get the 0-indexed position on reference and sequence so we can roll the edges
                            ins_len=insertion.end()-insertion.start()
                            # when calculating the start position on reference, strip any prior insertions from CIGAR
                            # since these are not present in reference
                            # -1 to 0-index start postion, -1 since the len includes start position
                            ins_pos=int(start)-1+len(longcigar[:insertion.start()].replace('I','').strip('S'))

                            # when getting the sequence of the insertion, we have to use the aligned read
                            # since by definition insertions are not in the reference
                            # however, deletions are not in the aligned read, so we have to adjust for that
                            del_adjust=longcigar[:insertion.start()].count('D')
                            ins_seq=alnseq[insertion.start()-del_adjust:insertion.end()-del_adjust]

                            insertions[ins_len].update(('-'.join((str(ins_pos), ins_seq)),))

                        for deletion in del_pattern.finditer(longcigar):
                            # first get the 0-indexed position on reference and sequence so we can roll the edges
                            del_len=deletion.end()-deletion.start()
                            # when calculating the start position on reference, strip any prior insertions from CIGAR
                            # since these are not present in reference
                            del_pos=int(start)-1+len(longcigar[:deletion.start()].replace('I','').strip('S'))
                            del_seq=refs[sample.strip('_np')]['fwd'][del_pos:del_pos+del_len]

                            deletions[del_len].update(('-'.join((str(del_pos), del_seq)),))

            
            print(f'\rrealigning {sample} replicate {rep_num+1:<50}')

            # now that we have a list of observed indels and their count
            # we can proceed to realign the indels. this way we realign
            # each unique indel only once, no matter the count
            for ins_len in sorted(insertions.keys(), reverse=True):
                for insertion, count in insertions[ins_len].most_common():

                    # separate insertion position and sequence again
                    ins_pos, ins_seq=insertion.split('-')
                    ins_pos=int(ins_pos)

                    # cool now let's roll
                    edges=roll_indel(ins_seq,refs[sample.strip('_np')]['fwd'],ins_pos)

                    if edges != None:
                        # this is a duplication, figure out best position
                        dGs_f = all_dGs_f[edges[0]:edges[1]+1]
                        # gotta correct with +ins_len-1 for reverse to account for RdRp moving other way
                        dGs_r = all_dGs_r[edges[0]+(ins_len-1):edges[1]+1+(ins_len-1)]

                        if sum(dGs_f)==0 and sum(dGs_r)==0: # no stems available for either orientation, if this is a structureless mutant, align like structured mutant
                            if sample.strip('_np') in destruct.keys(): 
                                dGs_f = all_dGs_f_struc[edges[0]:edges[1]+1]
                                dGs_r = all_dGs_r_struc[edges[0]+(ins_len-1):edges[1]+1+(ins_len-1)]

                        if min(dGs_f) < min(dGs_r): # take forward position, 1-index
                            ins_pos = edges[0]+dGs_f.index(min(dGs_f))+1
                            ins_seq= refs[sample.strip('_np')]['fwd'][ins_pos-1:ins_pos-1+ins_len]
                            ins_ori='fwd'
                        else: # take reverse position, 1-index
                            ins_pos = edges[0]+(ins_len-1)+dGs_r.index(min(dGs_r))+1
                            ins_seq= refs[sample.strip('_np')]['fwd'][ins_pos-1-(ins_len-1):ins_pos-1+1]
                            ins_ori='rev'
                        # determining insertion type by checking length and if all characters are the same
                        if len(ins_seq)==1:
                            ins_type='single-nucleotide'
                        elif len(set(ins_seq))==1: 
                            ins_type='homopolymer'
                        else:
                            ins_type='heteropolymer'
                    else:
                        # can't roll, this is complex, just update the type and ori and leave start and seq as is
                        ins_type='complex'
                        ins_ori='unk'

                    # save realigned insertion
                    ins_freq=count/coverage[int(ins_pos)-1]*1000
                    fins.write('\t'.join(map(str,(sample, rep_num+1, replicate[0], replicate[1], ins_len, ins_pos, ins_seq, ins_ori, ins_type, count, ins_freq)))+'\n')

            for del_len in sorted(deletions.keys(), reverse=True):
                for deletion, count in deletions[del_len].most_common():
                    del_pos, del_seq=deletion.split('-')
                    del_pos=int(del_pos)

                    # cool now let's roll
                    edges=roll_indel(del_seq,refs[sample.strip('_np')]['fwd'],del_pos)

                    # ! for a deletion, the pol has to skip ahead, this means we have to stop before replicating the deleted part
                    # roll_edges reports the 5' position of the rolled sequence. this is the 5' position of the deleted sequence in forward, 
                    # but the 3' position of the deleted sequence in reverse. we want to stop right before the 3', so:
                    # for forward that means going from edges[0]+del_len to edges[1]+del_len
                    # for reverse that means going from edges[0]-1        to edges[1]-1
                    dGs_f = all_dGs_f[edges[0]+del_len:edges[1]+del_len+1]
                    # we don't need a correction factor here, it's built into the range 
                    dGs_r = all_dGs_r[edges[0]-1:edges[1]]
                                                    
                    if sum(dGs_f)==0 and sum(dGs_r)==0: # no stems available for either orientation, if this is a structureless mutant, align like structured mutant
                        if sample in destruct.keys(): # there is a structures sequence we want to compare this to, fold that sequence and use those positions
                            dGs_f = all_dGs_f_struc[edges[0]+del_len:edges[1]+del_len+1]
                            dGs_r = all_dGs_r_struc[edges[0]-1:edges[1]]
                    
                    if min(dGs_f) < min(dGs_r): # take forward position
                        del_pos = edges[0]+del_len+dGs_f.index(min(dGs_f))
                        del_seq= refs[sample.strip('_np')]['fwd'][del_pos-del_len:del_pos]
                        del_ori='fwd'
                    else: # take reverse position
                        del_pos = edges[0]-1+dGs_r.index(min(dGs_r)) 
                        del_seq= refs[sample.strip('_np')]['fwd'][del_pos+1:del_pos+del_len+1]
                        del_ori='rev'

                    # determining deletion type
                    if edges[0]!=edges[1]: 
                        if len(del_seq)==1:
                            del_type='ambiguous snd'
                        else:
                            del_type='ambiguous'
                    else:
                        if len(del_seq)==1:
                            del_type='unambiguous snd'
                        else:
                            del_type='unambiguous'
                    
                    del_freq=count/coverage[int(del_pos)-1]*1000
                    fdel.write('\t'.join(map(str,(sample, rep_num+1, replicate[0], replicate[1], del_len, del_pos, del_seq, del_ori, del_type, count, del_freq)))+'\n')

            # write coverage to file
            with open(paths['coverage'], 'w') as f:
                for p, n in enumerate(coverage):
                    f.write(f'{p+1}\t{n}\n')

        
