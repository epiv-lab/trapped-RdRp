# Slidingfold

Slidingfold was designed to predict potential RdRp-trapping RNA structures. To try and predict transient rather than static RNA structures it uses a sliding-window approach inspired by the approach from [French and Pitré et al 2022](https://www.science.org/doi/10.1126/sciadv.abp8655).

Nucleotides within an identically-sized window up- and downstream of the RdRp footprint, corresponding to the exiting and entering template strands respectively, are used for sequence prediction. 

To mimic the close proximity of the entering 5’ and exiting 3’ template strand ends, the nucleotides in the RdRp footprint are replaced by an arbitrary sequence (```AAGGGGGGAAAACCCCCCAA```), constrained to fold as a stem-loop (constraint: ```xx((((((xxxx))))))xx```, where ```x```s denote nucleotides forbidden to pair, while open parentheses are forced to pair with the corresponding closed parentheses). 

Folding is then performed using the ```vrna_mfe``` function while disallowing lonely base pairs (by default). After folding, the minimum free energy value is corrected by the minimum free energy value obtained with the arbitrary stem alone.

Structures are then parsed to identify different base pairs and stored in dot-bracket format with ```()``` representing ```A-U```, ```{}``` representing ```G-C```, and ```<>``` representing ```G-U``` base pairs.

## Dependencies
Slidingfold was written and tested using Python 3.9 on MacOS 14.0 and uses the following software (all available via conda):
* [ViennaRNA](https://github.com/ViennaRNA) for RNA structure prediction (using version 2.4.18)
* [Biopython](https://biopython.org/) for parsing sequences from Fasta files (using version 1.79)
* [Plotly](https://github.com/plotly/plotly.py) for generating graphs (HTML or PDF) (using version 5.15.0)
* [Jinja](https://jinja.palletsprojects.com/en/3.1.x/) for generating interactive HTML outputs (using version 2.11.2)

Install time should be less than a minute if you know how to use conda.

The navigation menu of the interactive HTML output uses jQuery and requires and active internet connection for full functionality.
## Usage
```bash
python slidingfold.py <input> 
```

Two files are required:
* the input file, which can contain any number of sequences in fasta format. This will run an analysis with default parameters and output an interactive HTML file in the same folder as the input, with the name based on the input file.  
* the file "sample_start_positions.tsv", which indicates the parameters used for the default zoom, i.e the start and lenght of the window. This file should be in the same directory as the script and should have the following format:
# Default values: <default_start_position> <default_view_length>
# seq_id <start_position> <view_length>
seq1 10 200
seq2 15 250
seq3 5 150
The first line is a comment and defines the default values for start position and view length. The subsequent lines define the start position and view length for each sequence ID. Sequence IDs should correspond to the IDs in the FASTA file. The default values will be used if a sequence ID is not found in the TSV file.

The jinja template provided here (or your own template) must be in the same folder as the python file.

### Optional arguments:
#### Model parameters:
* ```-fp <int>:<int>``` or ```--footprint <int>:<int>``` size of the polymerase footprint and active site location, separated by a colon, can be comma-separated list, e.g.: ```-fp 20:5,22:6``` will perform analysis for a footprint of 20 with active site in position 5 as well as a footprint of 22 with active site in position 6 (default 20:5)
*  ```-l <int>``` or ```--local <int>``` size of the window to use for folding, can be comma-separated list, can specify ranges using hyphens, e.g.: ```-l 5,10-15``` will perform the analysis for window sizes 5, 10, 11, 12, 13, 14, and 15 (default 10)
* ```-lp``` or ```--lonely_pairs``` allow lonely RNA base pairs to form
* ```-rb``` or ```--roadblock``` use the roadblock model instead of the trapped-RdRp model, in effect this sets the size of the exiting template window to 0, predicting only RNA structures forming upstream of the RdRp
#### Output options
* ```-t``` or ```--text``` output a text file with all stems instead the interactive html
* ```-pdf``` or ```--pdf``` output a static pdf file instead of the interactive html
* ```-noAA``` or ```--no_amino_acid``` do not display the amino acid sequence in the pdf or html output
* ```-rf <1,2,3>``` or ```--reading_frame <1,2,3>``` specify the reading frame to be used for translation, if no reading frame is provided, the longest possible ORF for each sequence will be determined and used instead
* ```-tag <str>``` or ```--output_tag <str>``` suffix to add to the generated file
* ```-le <float>``` or ```-lower_edge <float>``` manually sets the lower edge of the heatmap scale
* ```-ue <float>``` or ```-upper_edge <float>``` manually sets the upper edge of the heatmap scale

## Demo
An example output file for the A/mallard/Netherlands/3/1999 HA is provided in the ```demo``` folder. To replicate this file yourself, run slidingfold with default parameters on the fasta input file provided in the same folder:

```python slidingfold.py <path to nl-rkkr.fasta>```

If all goes well, the resulting nl-rkkr_slidingmaps.html in the same folder as the input file should look identical to the one in the ```demo``` folder. Runtime should be a few seconds.

## Reproducing Supplementary Fig. 1
To reproduce the predictions show in Supplementary Fig.1, run the script on the fasta file Supplementary Data 1 for window sizes of 10 and 15.

```python slidingfold.py <path to supplementary_file-1.fasta> -l 10,15```