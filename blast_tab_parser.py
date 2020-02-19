# Get the best Blast hits from a tab file based on their Bitscore and Identity % (in case of tie keeps all hits).
# e.g. python3 blast_tab_parser.py input.tab output.tab

import sys


def parse_file(input_file):
    inputFile = open(input_file, "r")
    loci = []
    matches = []
    bitscores = []
    save_lines = dict()

    for line in inputFile:
        line = line.split("\t")

        locus = line[0]

        if locus not in loci:
            loci.append(locus)

            if len(matches) > 0:
                matches_best_score = [matches[x] for x in range(len(bitscores))
                                      if bitscores[x] == max(bitscores)]

                save_lines[loci[-2]] = [max(matches_best_score), max(bitscores)]

            bitscores = []
            matches = []
            matches.append(line[2])
            bitscores.append(line[11])

        else:
            matches.append(line[2])
            bitscores.append(line[11])

    # last line exception
    matches_best_score = []
    for i in range(len(bitscores)):
        if bitscores[i] == max(bitscores):
            matches_best_score.append(matches[i])
    save_lines[loci[-1]] = [max(matches_best_score), max(bitscores)]

    inputFile.close()
    return save_lines


def write_output(input_file, saveLines, outputFile):
    inputFile = open(input_file, "r")
    output = open(outputFile, "w")

    for line in inputFile:
        line = line.split("\t")
        if line[2] == saveLines[line[0]][0] and line[11] == saveLines[line[0]][1]:
            output.write("\t".join(line))

    output.close()


if __name__ == "__main__":
    try:
        write_output(sys.argv[1], parse_file(sys.argv[1]), sys.argv[2])
    except IndexError:
        write_output(sys.argv[1], parse_file(sys.argv[1]), "best_hits.tab")
