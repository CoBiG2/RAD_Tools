#!/usr/bin/python3


def main(gestefile):
    """
    Reads a .geste file and changes "1 X" to "2 X 0" for consistency.
    Prints everything to stdout, so you ight want to use a shell redirect.
    """
    infile = open(gestefile, 'r')
    for lines in infile:
        slines = lines.split()
        try:
            if slines[2] == "1":
                lines = lines.rstrip() + " 0 \n"
        except IndexError:
            pass
        print(lines, end="")

if __name__ == "__main__":
    from sys import argv
    # Usage: python3 correct_geste.py file.geste
    main(argv[1])
