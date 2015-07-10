#!/usr/bin/python3

# Copyright 2015 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of BioGenepop.
# BioGenepop is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# BioGenepop is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with BioGenepop. If not, see <http://www.gnu.org/licenses/>.


from Bio.PopGen.GenePop.EasyController import EasyController


def get_exp_obs_het(filehandle, outfile_name):
    """Doc here"""
    pop_names, loci_names = filehandle.get_basic_info()

    outfile=open(outfile_name,'w')
    double_loci = [loci_names[i//3] for i in range(len(loci_names)*3)]
    outfile.write("\t" + "\t".join(double_loci) + "\n")
    for pop in range(len(pop_names)):
        hetros = pop_names[pop] + "\t"
        loci_map = filehandle.test_hw_pop(pop, "excess")
        for locus in loci_names:
            exp_homo, obs_homo, exp_hetero, obs_hetero = filehandle.get_heterozygosity_info(pop,locus)
            if loci_map[locus] is not None:
                hetros += str(loci_map[locus][0]) + "\t" + str(exp_hetero) + "\t" + str(obs_hetero) + "\t"
            else:
                hetros += "-\t" + str(exp_hetero) + "\t" + str(obs_hetero) + "\t"
        hetros = hetros.rstrip() + "\n"
        outfile.write(hetros)
    outfile.close()

def filehandler(infile):
    """Doc here"""
    filehandle = EasyController(infile)

    return filehandle

if __name__ == "__main__":
    # Usage: python3 BioGenepop.py infile outfile
    from sys import argv
    filehandle = filehandler(argv[1])
    get_exp_obs_het(filehandle, argv[2])
