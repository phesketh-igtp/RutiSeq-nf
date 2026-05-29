#!/usr/bin/env python3

###############################################################################
# Source: https://gitlab.com/LPCDRP/illumina-blindspots.pub/-/blob/master/src/sequence_processing/illumina-single-isolate.mk?ref_type=heads
# Function:
#   This script calculates the average sequencing depth from an mpileup file and identifies genomic
#   positions with unusually low coverage relative to that average.
# Workflow:
#   1) Read mpileup file:
#      - Extract genomic position (coordinate) and coverage depth (column 4)
#      - Store values in a dictionary: {position → coverage}
#   2) Calculate mean coverage:
#      - Compute the average coverage across all positions
#      - Write mean coverage to output file: <sample>.average_coverage.txt
#   3) Define low-coverage threshold:
#      - Minimum coverage (d) = 5
#      - Reference average coverage (c_ref) = 50
#      - Compute relative threshold:
#            d_i = (d / c_ref) * mean_coverage
#      - This allows the threshold to scale with sample quality
#   4) Identify low-coverage positions:
#      - Positions where coverage ≤ d_i are classified as low coverage
#   5) Write low-coverage positions to file:
#      - Output file: <sample>.LC_positions.txt
#      - Format: position <tab> coverage
# Key assumptions / limitations:
#   a) Uses only column 4 (coverage depth) from the mpileup file
#   b) Ignores read bases and base quality information
#   c) Positions are stored as strings (not integers)
#   d) Coverage threshold is relative, not fixed
###############################################################################

import argparse

###############################################################################
# Parse arguments
###############################################################################

parser = argparse.ArgumentParser(description="""Finds positions in mpileup that have low coverage""")
parser.add_argument("-m", "--mpileup", required = True, help="mpileup file to be analyzed")
args = parser.parse_args()
mpileup = args.mpileup

###############################################################################
# Workflow
###############################################################################

def main():
    position_dict = create_position_dict(mpileup)
    avg_cov = mean_coverage(position_dict)
    coverage_file = write_mean_coverage(avg_cov)
    positions_file = lc_positions(position_dict,avg_cov)
    return positions_file, coverage_file

def create_position_dict(file):
    """Takes an mpileup file and parses the columns. Creates dictionary with positions as key and
    depth of coverage as value."""
    position_dict = {}
    with open(file,'r') as mpileup:
        for line in mpileup:
            column = line.split()
            name = column[0]
            coordinate = column[1]
            ref_base = column[2]
            reads = column[3]
            position_dict[coordinate] = float(reads)
    return position_dict

def mean_coverage(mpileup_dict):
    """Takes a dictionary with position as key and depth of coverage as value.
        Returns the mean coverage value of specific isolate."""
    # list of coverage values
    k_values = list(mpileup_dict.values())
    k_values = list(map(float, k_values)) ##k_values = map(float, k_values)
    # mean coverage of specific isolate (i)
    c_iso = sum(k_values) / len(k_values) ##c_iso = sum(k_values) / len(k_values)
    return c_iso

def write_mean_coverage(c_iso):
    """Takes average coverage of isolate and writes it to a new text file."""
    filename = mpileup.split('.')
    SRR = filename[0]
    newfilename = SRR + '.average_coverage.txt'
    with open(newfilename, 'w') as avg_coverage:
        avg_coverage.write(str(c_iso) + '\n')
    return avg_coverage

def lc_positions(mpileup_dict, c_iso):
    """Takes a dictionary with position as key and depth of coverage as value.
        Returns a text file with low coverage positions."""
    # minimum coverage
    d = float(5)
    # mean coverage for a typical sequencing run
    c_ref = float(50)
    # relative coverage
    D = d / c_ref
    # low coverage threshold for i
    d_i = D*c_iso
    #  of low coverage positions in isolate
    klc_dict = {}
    for position in mpileup_dict:
        cov = float(mpileup_dict[position])
        if cov <= d_i:
            klc_dict[position] = int(cov)
    # write each position and coverage to a new line in a text file
    low_cov_string = ''
    for position in klc_dict:
        low_cov_string += str(position) + '\t' + str(klc_dict[position]) + '\n'
    filename = mpileup.split('.')
    SRR = filename[0]
    newfilename = SRR + '.LC_positions.txt'
    with open(newfilename, 'w') as LC_positions:
        LC_positions.write(low_cov_string)
    return LC_positions

if __name__ == '__main__':
    main()


