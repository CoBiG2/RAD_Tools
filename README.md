# RAD Tools

This repository contains several tools built by CoBiG² members to help us deal with RAD and GBS data.
The existing tools are mentioned below.

## Baypass_workflow.R

This script has moved [to a new repository](https://github.com/StuntsPT/pyRona)!

## geste2baypass.py

This script will convert a file from the GESTE format to BayPass input.

## VCF2phy

This is a converter from VCF to phylip format. ~~It currently does not create the phylip header, you have to do it yourself.~~ Headers are created fine too.

## FDR.R

This is a script made in R used to calculate the *Q*-values from the *p*-values csv table that is outputed from [lositan](http://popgen.net/soft/lositan/).
The input is the mentioned csv, and a new csv is outputed, with an extra *q*-values column.
That simple.


## struct_to_distruct.py

This is an *ugly* script that will transform an output from [STRUCTURE](http://web.stanford.edu/group/pritchardlab/structure.html) into an input for [DISTRUCT](http://web.stanford.edu/group/rosenberglab/distruct.html).
You **need** to have a trailing "/" at the end of the output location, or you'll get strangely named files in the parent directory of the one you are pointing at...


## Loci_counter.py

A simple script to add a number to each loci in *.loci* files. Outputs to *stdout*, so "pipe away"!


## structure_filter.py

Another simple script that will filter the **columns** of a .structure file according to a "subset" file that contains nothing but the wanted marker names, one per line.

## bayenv2_results_miner.py

A simple script to look mine the data files produced by Bayenv2. It will look for any SNPs which have a BF>=10 AND a spearman p-value < 0.05 (these values can be changed in the source).

## filter_loci_by_vcf.py

This script will discard any information in a .loci file that is not referenced in the provided vcf file.

## singleton_site_remover.py

This script will remove any singleton sites from a "phylip" formatted file.
Takes one argument - the phylip file. The output will be written to a new file
which has the same name as the original, but with the prefix
"\_no_singletons.phy".


## correct_geste.py

This script will correct some inconsistencies in GESTE files created by PGDSpider.
Sometimes, when all alleles are of the same type, the get presented in the geste files as:

```
## X    2   X 0
```

Other time they show up instead like this:

```
## X    1   X
```

This script modifies the file so that they are always presented in the first form.
It reads a GESTE file and spits the output to STDOUT, so you might want to use a shell redirect.

## vcf2dadi.py

This script will convert any VCF file into [dadi](https://bitbucket.org/gutenkunstlab/dadi) input.
It take 3 arguments (input, output and population_information).
The "population_information" argument should be a text file with one line per population.
Each line should start with the population name, followed by ":" followed by the names of the samples that belong to the referred population, separated by whitespace, like this:

    Pop1:Sample1 Sample2 Sample3
    Pop2:Sample4 Sample5 Sample6 Sample7

This script requires the [PyVCF](https://github.com/jamescasbon/PyVCF) pacakge

## License

Everything is under the GPLv3.
