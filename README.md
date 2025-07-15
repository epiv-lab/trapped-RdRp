# The trapped-RdRp model

Trapping of the influenza virus RNA-dependent RNA polymerase (RdRp) as the molecular mechanism underlying recurrent insertions in influenza virus genomes was proposed in [Gultyaev et al 2019](https://academic.oup.com/ve/article/doi/10.1093/ve/vez034/5552811)[^1]. The model states that insertion hotspots can occur when transient RNA structures form around the influenza virus RdRp during replication of insertion-prone sequences. The physical constraint on the RdRp leads to an increase in error rates. Insertion hotspots therefore result from a synergistic combination of error-prone RNA sequences and error-enhancing RNA structures.

This repository is intended to regroup scripts written to support study of the formation of RdRp-trapping RNA structures and their impact on RdRp errors during influenza genome replication in the workgroup of [Mathilde Richard](https://www.erasmusmc.nl/en/research/researchers/richard-mathilde).

A brief overview of the scripts and their purpose:

## Slidingfold
The slidingfold script is used to predict the formation of transient RdRp-trapping RNA structures, based on the sliding window approach proposed in [French and Pitré et al 2022](https://www.science.org/doi/10.1126/sciadv.abp8655)[^2]. We slide the RdRp nucleotide-by-nucleotide along the template and use the [ViennaRNA](https://github.com/ViennaRNA) package to predict folding of small upstream and downstream sequence parts. 

Results can be output as simple text files, pdf images, or interactive graphs thanks to [Plotly](https://github.com/plotly/plotly.py) and [Jinja](https://jinja.palletsprojects.com/en/3.1.x/).

## Cirseq-analysis
In Funk et al 2024[^3], we use a CirSeq-based approach to investigate the impact of RNA sequence and structure on insertion frequencies at the HA cleavage site. The NGS data is first processed using a slightly modified version of the original [CirSeq](https://github.com/ashleyacevedo/CirSeq)[^4] script and is then further analyzed using a custom script to
* Identify UMIs and remove PCR duplicates (thanks to [Edlib](https://github.com/Martinsos/edlib))
* Identify all possible alignments of insertions/deletions
* Realign all insertion/deletions to the most stable RdRp-trapping structure possible (thanks to [ViennaRNA](https://github.com/ViennaRNA))
* Calculate position-by-position coverage of the reference

## References and notes
[^1]: A. P. Gultyaev, M. Richard, M. I. Spronken, R. C. L. Olsthoorn, R. A. M. Fouchier, Conserved structural RNA domains in regions coding for cleavage site motifs in hemagglutinin genes of influenza viruses. Virus Evol. 5, 1–11 (2019).
[^2]: H. French, E. Pitré, M. S. Oade, E. Elshina, K. Bisht, A. King, D. L. V. Bauer, A. J. W. te Velthuis, Transient RNA structures cause aberrant influenza virus replication and innate immune activation. Sci. Adv. 8, 1–11 (2022).
[^3]: M. Funk, M. I. Spronken, T. M. Bestebroer, A. C. M. de Bruin, A. P. Gultyaev, R. A. M. Fouchier, A. J. W. te Velthuis, M. Richard, Transient RNA structures underlie highly pathogenic avian influenza virus genesis.
[^4]: A. Acevedo, R. Andino, Library preparation for highly accurate population sequencing of RNA viruses. Nat. Protoc. 9, 1760–1769 (2014).