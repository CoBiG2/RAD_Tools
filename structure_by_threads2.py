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

#Usage: python3 structure_pipe.py K reps infile outpath num_of_threads
#Where "K" is the number of "Ks" to test, "reps" is the number of replicates
#and num_of_threads is the number of threads.

import subprocess
from multiprocessing import Pool

##Declare some global variables:
#Where is structure?
structure_bin = "/opt/structure/bin/structure"


def RunProgram(iterations):
    """Runs external programs and deals with their output."""
    K, rep_num = iterations
    program_stdout = []
    cli = [structure_bin, "-K", str(K), "-i", infile, "-o", outpath +  "/K" + str(K) + "_rep" + str(rep_num)]
    try:
        program = subprocess.Popen(cli, bufsize=64,shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        
        for lines in program.stdout:
            lines = lines.decode("utf-8").strip()
            print(lines)
            program_stdout.append(lines)
    except:
        quit("\nERROR:Program not found... exiting. Check your configuration file.\n")



def structure_threader(Ks, replicates, threads):
    pool = Pool(threads)
    jobs = []
    for k in Ks:
        for rep in replicates:
            jobs.append((k, rep))
    
    pool.map(RunProgram, jobs)
    pool.close()
    pool.join()



if __name__ == "__main__":
    from sys import argv
    #Number of K
    Ks = range(1, int(argv[1]) + 1)
    #Number of replicates
    replicates = range(1, int(argv[2]) + 1)
    infile = argv[3]
    outpath = argv[4]
    threads = int(argv[5])
    
    structure_threader(Ks, replicates, threads)
