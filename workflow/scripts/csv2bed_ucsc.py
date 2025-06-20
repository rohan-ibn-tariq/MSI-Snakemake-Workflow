#!/usr/bin/env python
"""
Script for converting pytrf .csv to .bed (UCSC filtering style)
-----------------------------------------------------------------------
Takes .csv file created by pytrf utility having repeat regions
and converts that into UCSC filtering style .bed.
#######UCSC STYLE BELOW################################################
field	    example	    SQL type	                info	    description
bin	        585	        smallint(5) unsigned	    range	    Indexing field to speed
                                                                chromosome range queries.
chrom	    chr1	    varchar(255)	            values	    Reference sequence chromosome or scaffold
chromStart	50481	    int(10) unsigned	        range	    Start position in chromosome
chromEnd	50513	    int(10) unsigned	        range	    End position in chromosome
name	    16xGT	    varchar(255)	            values	    Name of item
"""

import os
import sys
import csv
from argparse import ArgumentParser

from binning import assign_bin


def parse_args():
    """
    Argument parser that takes as an argument .csv path and as well as
    output .bed path.
    """
    parser = ArgumentParser(
        description="Convert pytrf .csv file to UCSC bin-compliant BED format"
    )
    parser.add_argument("--csv", required=True, help="Input .csv file")
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


def get_repeat_info(row):
    """
    Collects fields from pytrf CSV row for UCSC style filtering.
    pytrf CSV format: [chrom, start, end, motif, motif_length, repeat_number, repeat_length]
    """
    try:
        chrom = str(row[0])
        start = int(row[1]) - 1 
        end = int(row[2])
        motif = str(row[3])
        copies = int(float(row[5]))
        
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        
        name = f"{copies}x{motif}"
        return chrom, start, end, name
    except (KeyError, ValueError, TypeError, IndexError):
        return None


def main():
    """
    Parses pytrf .csv and writes in UCSC filtering style to .bed.
    """
    args = parse_args()
    records = []

    try:
        check_input_file(args.csv)
    except (FileNotFoundError, ValueError) as e:
        print(f"[csv2bed ERROR] {e}", file=sys.stderr)
        sys.exit(1)

    print("[csv2bed INFO] Starting Processing.")
    
    try:
        with open(args.csv, 'r', encoding="utf-8") as csvfile:
            reader = csv.reader(csvfile)

            for row in reader:
                info = get_repeat_info(row)
                if info:
                    chrom, start, end, name = info
                    bin_id = assign_bin(start, end)
                    records.append((bin_id, chrom, start, end, name))
                    
    except Exception as e:
        print(f"[csv2bed ERROR] Failed to parse CSV: {e}", file=sys.stderr)
        sys.exit(1)

    if not records:
        print(
            f"[csv2bed ERROR] No valid repeat records found in {args.csv}",
            file=sys.stderr,
        )
        sys.exit(1)

    with open(args.bed, "w", encoding="utf-8") as out:
        for bin_id, chrom, start, end, name in records:
            out.write(f"{bin_id}\t{chrom}\t{start}\t{end}\t{name}\n")

    print(f"[csv2bed INFO] Completed Processing. Wrote {len(records)} records.")


if __name__ == "__main__":
    main()