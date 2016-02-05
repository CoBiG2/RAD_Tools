#!/usr/bin/python2
# Copyright 2016 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of singleton_site_remover.
# singleton_site_remover is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# singleton_site_remover is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with singleton_site_remover.  If not, see <http://www.gnu.org/licenses/>.

def read_phylip(infile_name):
    """
    Reads a phylip file and retunrs the data in columns.
    """
    row_data = []
    infile = open(infile_name, 'r')
    infile.readline()  # Skip header
    for lines in infile:
        row_data.append([lines.split()[0]] + list(lines.strip().split()[1]))

    # Transpose the data
    column_data = zip(*row_data)
    infile.close()

    return column_data


def singleton_remover(columns):
    """
    Removes any columns that contain singleton sites. Returns any columns that
    don't contain sigleton sites.
    """
    no_singleton_cols = []
    iter_cols = iter(columns)
    names = next(iter_cols)
    no_singleton_cols.append([name + "\t" for name in names])
    for column in iter_cols:
        bases = set(column)
        if "-" in bases:
            bases.discard("-")
        if "N" in bases:
            bases.discard("N")
        base_counts = [column.count(x) for x in bases]
        try:
            if min(base_counts) > 1:
                no_singleton_cols.append(column)
        except ValueError:
            pass

    return no_singleton_cols


def new_phylip_writer(columns, infile_name):
    """
    Transposes the columns back into a phylip file.
    """
    back_to_rows = zip(*columns)
    outfile = open(infile_name[:-4] + "_no_singletons.phy", 'w')
    outfile.write(str(len(back_to_rows)))
    outfile.write(" ")
    outfile.write(str(len(back_to_rows[1])))
    outfile.write("\n")
    for line in back_to_rows:
        #print("".join(line))
        outfile.write("".join(line) + "\n")


    outfile.close()


def main(infile_name):
    """
    Main function. Runs all the others, based on the given input file name.
    """
    all_cols = read_phylip(infile_name)
    new_cols = singleton_remover(all_cols)
    new_phylip_writer(new_cols, infile_name)


if __name__ == "__main__":
    # Usage: python2 singleton_site_remover.py infile.phy
    from sys import argv
    main(argv[1])
