#!/usr/bin/python3

from sys import argv


def barcode_parser(barcodes_file):
    '''Reads and parses a barcode file and returns a dictionary with the data.
    '''
    infile = open(barcodes_file, 'r')
    codes = {}
    for lines in infile:
        lines = lines.split()
        codes[lines[1]] = lines[0]
    infile.close()

    return codes


def rev_comp(seq):
    '''Recieves a sequence and returns it's reversed and complemented
    counterpart.'''
    rv_seq = seq[::-1].translate(str.maketrans('ATGC', 'TACG'))
    return rv_seq


def qual_translator(qual):
    '''Converts a string of ASCII quals into a list of integer quals and returns
    it.'''
    values = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^\
_`abcdefghijklmnopqrstuvwxyz{|}~"
    num_qual = []
    for i in qual.strip():
        num_qual.append(values.index(i))

    return num_qual


def seq_writer(outfile_name, outfiles_location, name, seq, qual):
    '''Writes sequences in fastq format to the respective files according to the
    barcodes.'''
    outfile = open(outfiles_location + outfile_name, 'a')
    outfile.write(name)
    outfile.write(seq)
    outfile.write("+\n")
    outfile.write(seq)

    outfile.close()


def seq_match(seq, bcodes, max_mm):
    '''Accounts for differences between barcodes/cut sites and the sequence.'''
    for code in bcodes:
        mismatch = [i == j for i, j in zip(seq, code)]
        if mismatch.count(False) <= max_mm:
            return code
    else:
        return "nomatch"


def fastq_parser(fastq_file, codes):
    '''Parses a fastq file and demultiplexes the sequences into several files
    named according to the respective barcodes.'''
    fastq = open(fastq_file, 'r')
    line_no = 1
    for lines in fastq:
        if line_no == 1:
            name = lines
            line_no += 1
        elif line_no == 2:
            seq = lines
            line_no += 1
        elif line_no == 3:
            line_no += 1
        else:
            qual = lines
            num_qual = qual_translator(lines)
            line_no = 1
            if sum(num_qual)/len(num_qual) <= min_avg_qual:
                seq_writer("low_avg.fastq", outfiles_location, name, seq, qual)
            elif sorted(num_qual)[5] <= 10:
                seq_writer("too_many_poor.fastq", outfiles_location, name, seq,
                           qual)
            else:
                seq_writer(seq_match(seq, codes, max_mismatch) + ".fastq",
                           outfiles_location, name, seq, qual)
    fastq.close()


min_avg_qual = 20
outfiles_location = "/home/francisco/Desktop/fastq/split/"
max_mismatch = 2

if __name__ == "__main__":
    #Usage: python3 RAD_parser.py file.fastq file.barcodes
    barcodes = barcode_parser(argv[2])
    fastq_parser(argv[1], barcodes)
