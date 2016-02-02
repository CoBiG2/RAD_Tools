# RAD Tools

This repository contains several tools built by CoBiG² members to help us deal with RAD and GBS data.
The existing tools are mentioned below.


## VCF2phy

This is a converter from VCF to phylip format. ~~It currently does not create the phylip header, you have to do it yourself.~~ Headers are created fine too.


## RAD parser

This is a demultiplexer for RAD and GBS data. It is still incomplete, because it is still lacking a way to remove the barcodes and the enzyme cut sites. Should be done pretty soon.


## FDR.R

This is a script made in R used to calculate the *Q*-values from the *p*-values csv table that is outputed from [lositan](http://popgen.net/soft/lositan/).
The input is the mentioned csv, and a new csv is outputed, with an extra *q*-values column.
That simple.


## struct_to_distruct.py

This is an *ugly* script that will transform an output from [STRUCTURE](http://pritchardlab.stanford.edu/structure.html) into an input for [DISTRUCT](http://web.stanford.edu/group/rosenberglab/distruct.html).
You **need** to have a trailing "/" at the ond of the output location, ou you'll get stragely named files in the parent directory of the one you are pointing at...


## Loci_counter.py

A simple script to add a number to each loci in *.loci* files. Outputs to *stdout*, so "pipe away"!


## structure_filter.py

Another simple script that will filter the **columns** of a .structure file according to a "subset" file that contains nothing but the wanted marker names, one per line.

## bayenv2_results_miner.py

A simple script to look mine the data files produced by Bayenv2. It will look for any SNPs which have a BF>=10 AND a spearman p-value < 0.05 (these values can be changed in the source).

## filter_loci_by_vcf.py

This script will discard any information in a .loci file that is not referenced in the provided vcf file.

## License

Everything is under the GPLv3.
