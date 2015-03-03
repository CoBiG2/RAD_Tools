#!/usr/bin/python3

# Copyright 2015 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of structure_by_threads2.
# structure_by_threads2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# structure_by_threads2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with structure_by_threads2. If not, see <http://www.gnu.org/licenses/>.

# Usage: python3 structure_pipe.py K reps infile outpath num_of_threads
# Where "K" is the number of "Ks" to test, "reps" is the number of replicates
# and num_of_threads is the number of threads.

import subprocess
from multiprocessing import Pool

# Where is structure?
structure_bin = "/opt/structure/bin/structure"


def runprogram(iterations):
    """Run each structure job."""
    K, rep_num = iterations
    cli = [structure_bin, "-K", str(K), "-i", infile, "-o", outpath + "/K"
           + str(K) + "_rep" + str(rep_num)]
    print("Running: " + " ".join(cli))
    program = subprocess.Popen(cli, bufsize=64, shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)

    # This loop does nothing, but seems to be required for the program to fully
    # run all the jobs...
    for lines in program.stdout:
        pass


def structure_threader(Ks, replicates, threads):
    """Do the threading book-keeping to spawn jobs at the asked rate."""
    pool = Pool(threads)
    jobs = []
    # Ugly, ugly code. But since we are just dealing with counts here, we can
    # just leave it like that...
    for k in Ks:
        for rep in replicates:
            jobs.append((k, rep))

    pool.map(runprogram, jobs)
    pool.close()
    pool.join()
    print("All jobs finished.")


if __name__ == "__main__":
    from sys import argv
    # Number of K
    Ks = range(1, int(argv[1]) + 1)
    # Number of replicates
    replicates = range(1, int(argv[2]) + 1)
    infile = argv[3]
    outpath = argv[4]
    threads = int(argv[5])

    structure_threader(Ks, replicates, threads)
