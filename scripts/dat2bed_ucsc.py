#!/usr/bin/env python

from argparse import ArgumentParser
from binning import assign_bin


def parse_args():
    parser = ArgumentParser(
        description='Convert TRF .dat file to UCSC bin-compliant BED format'
    )
    parser.add_argument('--dat', required=True, help='Input .dat file')
    parser.add_argument('--bed', required=True, help='Output .bed file')
    return parser.parse_args()

def get_repeat_info(splitline):
    try:
        start = int(splitline[0]) - 1
        end = int(splitline[1])
        copies = splitline[2]
        motif = splitline[13]
        name = f"{copies}x{motif}"
        return start, end, name
    except (IndexError, ValueError):
        return None

def main():
    args = parse_args()
    chrom = None
    records = []

    with open(args.dat) as infile:
        for line in infile:
            if line.startswith("Sequence:"):
                chrom = line.strip().split()[1]
                if not chrom.startswith("chr"):
                    chrom = f"chr{chrom}"
                continue
            splitline = line.strip().split()
            if not splitline or not splitline[0].isdigit():
                continue
            info = get_repeat_info(splitline)
            if info:
                start, end, name = info
                bin_id = assign_bin(start, end)
                records.append((bin_id, chrom, start, end, name))

    with open(args.bed, 'w') as out:
        for bin_id, chrom, start, end, name in records:
            out.write(f"{bin_id}\t{chrom}\t{start}\t{end}\t{name}\n")

if __name__ == '__main__':
    main()
