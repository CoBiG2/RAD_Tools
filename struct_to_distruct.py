#!/usr/bin/python3
# Copyright 2015 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of struct_to_distruct.
# struct_to_distruct is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# struct_to_distruct is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with struct_to_distruct. If not, see <http://www.gnu.org/licenses/>.

#Creates indfile and popfile from structure output.
#Useage: python struct_to_distruct.py "structure_output_file" "directory_to_create_popfile_and_indfile"

from sys import argv
import re

infile = open(argv[1],'r')
popfile = open(argv[2] + 'popfile','w')
indfile = open(argv[2] + 'indfile','w')

indzone = 0

for lines in infile:
    if re.match("\d+:", lines.strip()):
        popfile.write(lines)
    elif lines.strip().endswith("Inferred clusters"):
        indzone = 1
        continue
    if indzone == 1:
        indfile.write(lines)
    if indzone == 1 and lines.startswith("\n"):
        indzone = 0


infile.close()
popfile.close()
indfile.close()
