"""
Script for converting TRF .dat to .bed(ucsc filtering style)
-----------------------------------------------------------------------
Takes .dat file created by TRF utility having repeat regions
and converts that into ucsc filtering style .bed.
#######UCSC STYLE BELOW################################################
field	    example	    SQL type	                info	    description
bin	        585	        smallint(5) unsigned	    range	    Indexing field to speed
                                                                chromosome range queries.
chrom	    chr1	    varchar(255)	            values	    Reference sequence chromosome or scaffold
chromStart	50481	    int(10) unsigned	        range	    Start position in chromosome
chromEnd	50513	    int(10) unsigned	        range	    End position in chromosome
name	    16xGT	    varchar(255)	            values	    Name of item
"""

#!/usr/bin/env python
import os
import sys
from argparse import ArgumentParser

from binning import assign_bin


def parse_args():
    """
    Argument parser that takes as an argument .dat path and as well as
    output .bed path.
    """
    parser = ArgumentParser(
        description="Convert TRF .dat file to UCSC bin-compliant BED format"
    )
    parser.add_argument("--dat", required=True, help="Input .dat file")
    parser.add_argument("--bed", required=True, help="Output .bed file")
    return parser.parse_args()


def check_input_file(path):
    """
    Check if input file exists and is non empty.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input file does not exist: {path}")
    if os.stat(path).st_size == 0:
        raise ValueError(f"Input file is empty: {path}")


def get_repeat_info(splitline):
    """
    Collects fields from .dat input line for ucsc style filtering.
    """
    try:
        start = int(splitline[0]) - 1
        end = int(splitline[1])
        copies = int(float(splitline[3]))
        motif = splitline[13]
        name = f"{copies}x{motif}"
        return start, end, name
    except (IndexError, ValueError):
        return None


def main():
    """
    Parses .dat and writes in ucsc filtering style to .bed.
    """
    args = parse_args()
    chrom = None
    records = []

    try:
        check_input_file(args.dat)
    except (FileNotFoundError, ValueError) as e:
        print(f"[dat2bed ERROR] {e}", file=sys.stderr)
        sys.exit(1)

    print("[dat2bed INFO] Starting Processing.")
    with open(args.dat, encoding="utf-8") as infile:
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

    if not records:
        print(
            f"[dat2bed ERROR] No valid repeat records found in {args.dat}",
            file=sys.stderr,
        )
        sys.exit(1)

    with open(args.bed, "w", encoding="utf-8") as out:
        for bin_id, chrom, start, end, name in records:
            out.write(f"{bin_id}\t{chrom}\t{start}\t{end}\t{name}\n")

    print("[dat2bed INFO] Completed Processing.")


if __name__ == "__main__":
    main()
