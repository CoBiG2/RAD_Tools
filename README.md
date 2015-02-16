#RAD Tools


This repository contains several tools built by CoBiG² members to help us dead with RAD and GBS data.
The existing tools are mentioned below.

##VCF2phy

This is a converter from VCF to phylip format. ~~It currently does not create the phylip header, you have to do it yourself.~~ Headers are created fine too.

##RAD parser

This is a demultiplexer for RAD and GBS data. It is still incomplete, because it is still lacking a way to remove the barcodes and the enzyme cut sites. Should be done pretty soon.


##FDR.R

This is a script made in R used to calculate the *Q*-values from the *p*-values csv table that is outputed from [lositan](http://popgen.net/soft/lositan/).
The input is the mentioned csv, and a new csv is outputed, with an extra *q*-values column. 
That simple.

##License

Everything is under the GPLv3.
