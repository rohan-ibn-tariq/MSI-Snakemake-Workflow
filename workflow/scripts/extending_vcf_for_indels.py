"""
Extension Script for MSI INDELS on MASON VCF
------------------------------------------------------
Take mason variated input .vcf file and generates new
.vcf by injecting indels that are MSI specific by
looking in MS regions found in .bed. Ignores variating
region already variated by mason or variating same MS
region multiple times.
"""

#!/usr/bin/env python3

import argparse
import os
import subprocess
import random
import re
from collections import defaultdict

import pysam
from binning import assign_bin


def parse_args():
    """
    Argument parser that takes optional and required parameters for
    running the script.
    """
    parser = argparse.ArgumentParser(
        description="Inject indels into VCF at microsatellite regions"
    )
    parser.add_argument(
        "--bed", required=True, help="Repeat regions BED file with UCSC-style bins"
    )
    parser.add_argument("--vcf-in", required=True, help="Input VCF (from Mason)")
    parser.add_argument(
        "--vcf-out", required=True, help="Output VCF with injected indels"
    )
    parser.add_argument(
        "--ref-fasta",
        required=True,
        help="Reference genome FASTA file (Index Should be Present)",
    )
    parser.add_argument(
        "--seed", type=int, default=50, help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--ins-rate", type=float, default=0.15, help="MSI indel's insertion rate"
    )
    parser.add_argument(
        "--del-rate", type=float, default=0.15, help="MSI indels deletion rate"
    )
    parser.add_argument(
        "--boost-rate",
        type=float,
        default=0.6,
        help="Overall probability of injecting into eligible motifs",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable detailed logging for each injection",
    )
    return parser.parse_args()


def normalize_chromosome(chrom):
    """
    Normalizing the name without chr if chr exists for proper injection.
    """
    return (
        chrom.removeprefix("chr") if chrom.startswith("chr") else chrom
    )


def load_repeats(bedfile):
    """
    Loads MS repeats from the .bed file.
    """
    if not os.path.exists(bedfile):
        raise FileNotFoundError(f"[INJ-INDEL ERROR] BED file not found: {bedfile}")

    repeats = defaultdict(lambda: defaultdict(list))
    with open(bedfile, encoding="utf-8") as file:
        for line in file:
            if line.strip() == "":
                continue
            fields = line.strip().split("\t")
            if len(fields) < 5:
                continue
            try:
                bin_id = int(fields[0])
                chrom = normalize_chromosome(fields[1])
                start = int(fields[2])
                end = int(fields[3])
                name = fields[4]
            except (ValueError, IndexError) as e:
                print(f"[INJ-INDEL WARNING] Skipping malformed line: {line.strip()} ({e})")
                continue
            repeats[chrom][bin_id].append((start, end, name))
    return repeats


def get_bin(pos):
    """
    Get UCSC style bin number.
    """
    return assign_bin(pos - 1, pos)


def add_msi_info_fields(header):
    """
    Adds MSI-specific INFO fields to the VCF header.
    Modifies the header in place.
    """
    header.info.add("MSI_TYPE", 1, "String", "Injected indel type: INS or DEL")
    header.info.add("MSI_MOTIF", 1, "String", "Repeat motif used for indel injection")
    header.info.add("MSI_ORIG_REPEATS", 1, "Integer", "Original repeat count from BED")
    header.info.add("MSI_CHANGE", 1, "Integer", "Number of motif units inserted or deleted")


def log_injection_start(ins_rate, del_rate, boost_rate):
    """
    Logs the start of the injection process and the effective rates.
    """
    print("[INJ-INDEL INFO] Starting injection process...")
    print(f"[INJ-INDEL INFO] Effective injection rate: {(ins_rate + del_rate) * boost_rate:.2%}")
    print(f"[INJ-INDEL INFO] Insertion rate: {ins_rate * boost_rate:.2%}")
    print(f"[INJ-INDEL INFO] Deletion rate: {del_rate * boost_rate:.2%}\n")
    print("[INJ-INDEL INFO] Phase 1: Processing existing variants!")


def log_phase2_start(total_existing):
    """
    Simple Phase 1 to 2 transition log.
    """
    print(f"[INJ-INDEL INFO] Phase 1 complete. Total existing variants: {total_existing}\n")
    print("[INJ-INDEL INFO] Phase 2: Injecting MSI indels at microsatellite regions...")


def sort_vcf_with_bcftools(vcf_out):
    """
    Sort the output VCF using bcftools.
    """
    print("[INJ-INDEL INFO]Sorting output VCF...")
    temp_sorted = vcf_out + ".sorted.vcf"
    try:
        subprocess.run(["bcftools", "sort", "-o", temp_sorted, vcf_out], check=True)
        os.rename(temp_sorted, vcf_out)
        print("[INJ-INDEL INFO] VCF sorted successfully")
    except subprocess.CalledProcessError:
        print(
            f"[INJ-INDEL WARNING] Could not sort VCF with bcftools. "
            f"Please sort manually:\n  bcftools sort -o sorted.vcf {vcf_out}"
        )
    except FileNotFoundError:
        print(
            f"[INJ-INDEL WARNING] WARNING: bcftools not found. "
            f"Please sort the VCF manually:\n  bcftools sort -o sorted.vcf {vcf_out}"
        )


def log_msi_summary(
    msi_candidates,
    msi_skipped_existing,
    msi_skipped_validation,
    msi_skipped_coord_mismatch,
    injected,
    total
):
    """
    Indel Injection Summary Report
    """
    print("\n[INJ-INDEL INFO]Finished MSI injection:")
    print(f"  Total microsatellite regions examined: {msi_candidates}")
    print(f"  Skipped due to existing variants: {msi_skipped_existing}")
    print(f"  Skipped due to validation failure: {msi_skipped_validation}")
    print(f"  Skipped due to coordinate mismatch: {msi_skipped_coord_mismatch}")
    print(f"  MSI indels injected: {injected}")
    print(f"  Total variants in output: {total + injected}")


def warn_no_indels_injected():
    """
    Warning if no indels injected.
    """
    print("\n[INJ-INDEL WARNING] No indels injected. Possible reasons:")
    print("  - Rate parameters too low (try increasing --ins-rate, --del-rate, or --boost-rate)")
    print("  - All microsatellite positions already have variants")
    print("  - Microsatellite regions don't match reference sequence")


def inject_indels(
    vcf_in,
    vcf_out,
    repeats,
    config,
):
    """
    Logic to inject MS indels on varying rates and write a new .vcf output file.
    It ignores varying regions already varied by the input .vcf.
    """
    ins_rate = config["ins_rate"]
    del_rate = config["del_rate"]
    seed = config["seed"]
    boost_rate = config["boost_rate"]
    ref_fasta_path = config["ref_fasta_path"]
    verbose = config.get("verbose", False)

    if not os.path.exists(vcf_in):
        raise FileNotFoundError(f"[INJ-INDEL ERROR] VCF input not found: {vcf_in}")

    fai_path = ref_fasta_path + ".fai"
    if not os.path.exists(fai_path):
        raise FileNotFoundError(f"[INJ-INDEL ERROR] Reference FASTA index not found: {fai_path}")

    random.seed(seed)

    vcf = pysam.VariantFile(vcf_in)
    add_msi_info_fields(vcf.header)

    ref_fasta = pysam.FastaFile(ref_fasta_path)

    outdir = os.path.dirname(vcf_out)
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    out = pysam.VariantFile(vcf_out, "w", header=vcf.header)

    injected = 0
    total = 0
    existing_variant_positions = set()
    existing_variant_ranges = defaultdict(list)

    log_injection_start(ins_rate, del_rate, boost_rate)

    # Phase 1: Write all existing variants and track their positions/ranges
    for rec in vcf:
        total += 1

        ######TODO: REMOVE ############################################
        # DEBUG: Print sample info for first record
        if total == 1:
            print(f"[DEBUG] First record format: {list(rec.format)}")
            print(f"[DEBUG] Sample names: {list(rec.samples.keys())}")
            for sample_name in rec.samples:
                sample_obj = rec.samples[sample_name]
                print(f"[DEBUG] Sample '{sample_name}' object: {sample_obj}")
        
            print(f"[DEBUG] Full record string: {str(rec)}")
        
            try:
                print(f"[DEBUG] Sample values: {sample_obj.values()}")
            except:
                pass
            
            try:
                print(f"[DEBUG] Sample items: {list(sample_obj.items())}")
            except:
                pass
            
            try:
                print(f"[DEBUG] Sample keys: {list(sample_obj.keys())}")
            except:
                pass
            
            try:
                print(f"[DEBUG] Sample phased: {sample_obj.phased}")
            except:
                pass
            
            try:
                print(f"[DEBUG] Sample alleles: {sample_obj.alleles}")
            except:
                pass
        ########################TODO:DEBUG REMOVE ABOVE##################################

        out.write(rec)

        # Track positions of existing variants (using normalized chromosome names)
        norm_chrom = normalize_chromosome(rec.chrom)
        existing_variant_positions.add((norm_chrom, rec.pos))

        var_start = rec.pos
        var_end = rec.stop
        existing_variant_ranges[norm_chrom].append((var_start, var_end))

        if total % 1000 == 0:
            print(f"[INJ-INDEL INFO] Processed {total} existing variants...")

    log_phase2_start(total)

    # Phase 2: Add new MSI indels at microsatellite regions without existing variants
    msi_candidates = 0
    msi_skipped_existing = 0
    msi_skipped_validation = 0
    msi_skipped_coord_mismatch = 0

    for chrom in repeats:
        for bin_id in repeats[chrom]:
            for start, end, name in repeats[chrom][bin_id]:
                msi_candidates += 1

                match = re.match(r"^(\d+)x([ACGT]+)$", name)
                if not match:
                    print(f"WARNING: Skipping unparseable motif '{name}'")
                    continue

                total_repeats_str, motif = match.groups()
                total_repeats = int(total_repeats_str)
                motif_len = len(motif)

                # Skip if motif is too long or only one repeat
                if motif_len > 10 or total_repeats <= 1:
                    continue

                expected_length = motif_len * total_repeats
                actual_length = end - start
                if actual_length != expected_length:
                    msi_skipped_coord_mismatch += 1
                    if verbose and msi_skipped_coord_mismatch:
                        print(f"\n{'='*80}")
                        print(
                            f"  COORDINATE MISMATCH at {chrom}:{start}-{end} ({name})"
                        )
                        print(
                            f"    Expected length: {expected_length}bp"
                            f" (from {total_repeats}x{motif})"
                        )
                        print(f"    BED length: {actual_length}bp")
                        print(f"\n{'='*80}")
                    end = start + expected_length

                msi_pos = start + 1
                msi_end = end

                overlap_found = False
                for var_start, var_end in existing_variant_ranges[chrom]:
                    if not (msi_end < var_start or msi_pos > var_end):
                        overlap_found = True
                        break

                if overlap_found:
                    msi_skipped_existing += 1
                    if verbose and msi_skipped_existing:
                        print(f"\n{'='*80}")
                        print(
                            f"  Skipping MS at {chrom}:{msi_pos}-{msi_end}"
                            f" ({name}) - overlaps existing variant"
                        )
                        print(f"\n{'='*80}")
                    continue

                # Validate that the reference sequence matches expected microsatellite
                repeat_span_start = start
                repeat_span_end = start + expected_length

                try:
                    ref_seq = ref_fasta.fetch(
                        chrom, repeat_span_start, repeat_span_end
                    ).upper()
                except Exception as e:
                    print(
                        f"[INJ-INDEL WARNING]: Failed to fetch reference at {chrom}:{repeat_span_start}-{repeat_span_end}: {e}"
                    )
                    continue

                expected_seq = motif * total_repeats
                if ref_seq != expected_seq:
                    msi_skipped_validation += 1
                    if (
                        verbose and msi_skipped_validation
                    ):
                        print(f"\n{'='*80}")
                        print(
                            f"VALIDATION FAIL at {chrom}:{repeat_span_start}-{repeat_span_end}"
                        )
                        print(
                            f"    Expected: {expected_seq[:50]}{'...' if len(expected_seq) > 50 else ''}"
                        )
                        print(
                            f"    Found:    {ref_seq[:50]}{'...' if len(ref_seq) > 50 else ''}"
                        )
                        print(f"    Motif: {motif} x {total_repeats}")
                        print(f"\n{'='*80}")
                    continue

                # Decide whether to inject an indel based on rates
                r = random.random()
                if r < (ins_rate + del_rate) * boost_rate:
                    indel_type = "INS" if r < ins_rate * boost_rate else "DEL"

                    # Determine the size of the change
                    if indel_type == "INS":
                        change = random.randint(1, 10)
                    else:
                        max_del = min(total_repeats - 1, 10)
                        change = random.randint(1, max_del) if max_del >= 1 else 1

                    # Create the indel sequence
                    motif_block = motif * change

                    if verbose:
                        log_lines = []
                        log_lines.append(f"\n{'='*80}")
                        log_lines.append(f"MSI INJECTION #{injected + 1}")
                        log_lines.append(f"{'='*80}")
                        log_lines.append("Microsatellite Region (from BED):")
                        log_lines.append(
                            f"  BED coords: {chrom}:{start}-{end} (0-based, adjusted)"
                        )
                        log_lines.append(f"  VCF position: {chrom}:{msi_pos} (1-based)")
                        log_lines.append(
                            f"  Motif info: {name} = {motif} repeated {total_repeats} times"
                        )
                        log_lines.append(f"  Expected length: {expected_length}bp")

                        # Show the original microsatellite sequence
                        ms_seq = ref_fasta.fetch(
                            chrom, start, start + expected_length
                        ).upper()
                        log_lines.append("\nReference Sequence:")
                        log_lines.append(
                            f"  Full MS: {ms_seq[:60]}{'...' if len(ms_seq) > 60 else ''}"
                        )

                        # Show with context (10bp on each side)
                        context_start = max(0, start - 10)
                        context_end = start + expected_length + 10
                        upstream = ref_fasta.fetch(chrom, context_start, start).upper()
                        downstream = ref_fasta.fetch(
                            chrom, start + expected_length, context_end
                        ).upper()
                        log_lines.append(
                            f"  With context: {upstream}[{ms_seq[:30]}{'...' if len(ms_seq) > 30 else ''}]{downstream}"
                        )

                        log_lines.append("\nInjected Variant:")
                        log_lines.append(f"  Type: {indel_type}")
                        log_lines.append(
                            f"  Change: {'+' if indel_type == 'INS' else '-'}{change} {motif} units"
                        )
                        log_lines.append(
                            f"  Sequence {'added' if indel_type == 'INS' else 'removed'}: {motif_block}"
                        )

                    new_rec = vcf.header.new_record()

                    if vcf.header.contigs:
                        first_contig = list(vcf.header.contigs)[0]
                        if first_contig.startswith("chr") and not chrom.startswith(
                            "chr"
                        ):
                            new_rec.chrom = "chr" + chrom
                        else:
                            new_rec.chrom = chrom
                    else:
                        new_rec.chrom = chrom

                    new_rec.id = f"MSI_{indel_type}_{chrom}_{msi_pos}"
                    new_rec.qual = None
                    new_rec.filter.add("PASS")

                    new_rec.info["MSI_TYPE"] = indel_type
                    new_rec.info["MSI_MOTIF"] = motif
                    new_rec.info["MSI_ORIG_REPEATS"] = total_repeats
                    new_rec.info["MSI_CHANGE"] = change

                    if indel_type == "INS":
                        ref_base = ref_fasta.fetch(chrom, end - 1, end).upper()
                        new_rec.pos = end
                        new_rec.stop = new_rec.pos
                        new_rec.ref = ref_base
                        new_rec.alts = (ref_base + motif_block,)

                        if verbose:
                            log_lines.append(f"\nVCF Record:")
                            log_lines.append(f"  CHROM: {new_rec.chrom}")
                            log_lines.append(f"  POS: {new_rec.pos}")
                            log_lines.append(f"  REF: {new_rec.ref}")
                            log_lines.append(
                                f"  ALT: {new_rec.alts[0][:20]}{'...' if len(new_rec.alts[0]) > 20 else ''}"
                            )
                            log_lines.append(
                                f"  Interpretation: Insert '{motif_block}' after microsatellite at position {end}"
                            )
                            log_lines.append(
                                f"  New repeat count: {total_repeats + change}"
                            )
                    else:
                        anchor_pos = (
                            end - change * motif_len
                        )
                        ref_start = anchor_pos - 1
                        ref_end = end

                        if ref_start < 0:
                            continue

                        ref_seq = ref_fasta.fetch(chrom, ref_start, ref_end).upper()
                        expected_len = 1 + (change * motif_len)
                        if len(ref_seq) != expected_len:
                            print(
                                f"WARNING: Reference length mismatch for deletion at {chrom}:{ref_start}-{ref_end}"
                            )
                            continue

                        new_rec.pos = ref_start + 1
                        new_rec.stop = new_rec.pos + len(ref_seq) - 1
                        new_rec.ref = ref_seq
                        new_rec.alts = (ref_seq[0],)

                        if verbose:
                            log_lines.append("\nVCF Record:")
                            log_lines.append(f"  CHROM: {new_rec.chrom}")
                            log_lines.append(f"  POS: {new_rec.pos}")
                            log_lines.append(
                                f"  REF: {new_rec.ref[:20]}{'...' if len(new_rec.ref) > 20 else ''}"
                            )
                            log_lines.append(f"  ALT: {new_rec.alts[0]}")
                            log_lines.append(
                                f"  Interpretation: Delete '{motif_block}' from end of microsatellite region (starting at {anchor_pos + 1})"
                            )
                            log_lines.append(
                                f"  New repeat count: {total_repeats - change}"
                            )

                    if verbose:
                        for line in log_lines:
                            print(line)

                    out.write(new_rec)
                    injected += 1

                    if not verbose and injected % 100 == 0:
                        print(f"[INJ-INDEL INFO]: Injected {injected} MSI indels...")

                    existing_variant_positions.add((chrom, new_rec.pos))
                    existing_variant_ranges[chrom].append((new_rec.pos, new_rec.stop))

    out.close()
    ref_fasta.close()
    vcf.close()
    sort_vcf_with_bcftools(vcf_out)
    log_msi_summary(
        msi_candidates,
        msi_skipped_existing,
        msi_skipped_validation,
        msi_skipped_coord_mismatch,
        injected,
        total,
    )

    if injected == 0:
        warn_no_indels_injected()


def main():
    """
    Main execution function with three calls technically doing three different jobs.
    1. Argument Parsing.
    2. Loading .bed data to Data structure
    3. Injecting MSI Indels.
    """
    args = parse_args()
    repeats = load_repeats(args.bed)
    config = {
        "ins_rate": args.ins_rate,
        "del_rate": args.del_rate,
        "seed": args.seed,
        "boost_rate": args.boost_rate,
        "ref_fasta_path": args.ref_fasta,
        "verbose": args.verbose,
    }
    inject_indels(
        args.vcf_in,
        args.vcf_out,
        repeats,
        config,
    )


if __name__ == "__main__":
    main()
