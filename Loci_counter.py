#!/usr/bin/python3

# Copyright 2015 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of Loci_counter.
# Loci_counter is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Loci_counter is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Loci_counter. If not, see <http://www.gnu.org/licenses/>.

#Usage:python3 Loci_counter.py file.loci file.counted.loci


def Loci_counter(loci_filename):
    """Counts the number of '//', adds that number to the string and prints it
    to stdout."""
    counter = 1
    infile = open(loci_filename, 'r')
    for lines in infile:
        if lines.startswith("//"):
            lines = lines.replace(" " * len(str(counter)), str(counter), 1)
            counter += 1
        print(lines, end="")


if __name__ == "__main__":
    from sys import argv
    Loci_counter(argv[1])
