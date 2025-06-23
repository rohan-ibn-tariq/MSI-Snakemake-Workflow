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
import random
import re
import subprocess
import sys
from collections import defaultdict

import pysam
from binning import assign_bin


def validate_probability(value, param_name, allow_zero=True):
    """Validate probability is in valid range [0, 1] or (0, 1]."""
    try:
        fvalue = float(value)

        if allow_zero:
            if fvalue < 0 or fvalue > 1:
                raise argparse.ArgumentTypeError(
                    f"{param_name} must be in range [0, 1]. Got: {fvalue}"
                )
        else:
            if fvalue <= 0 or fvalue > 1:
                raise argparse.ArgumentTypeError(
                    f"{param_name} must be in range (0, 1]. Got: {fvalue}"
                )

        return fvalue
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"{param_name} must be a number. Got: {value}"
        ) from exc


def validate_mono_penalty(value):
    """Validate mono_penalty is in valid range (0, 1]."""
    return validate_probability(value, "mono_penalty", allow_zero=False)


def validate_rate(value):
    """Validate rate is in valid range [0, 1]."""
    return validate_probability(value, "rate", allow_zero=True)


def validate_combined_rates(ins_rate, del_rate):
    """Validate that combined insertion and deletion rates don't exceed 1."""
    total_rate = ins_rate + del_rate
    if total_rate > 1:
        raise ValueError(
            f"Combined ins_rate + del_rate ({total_rate:.3f}) cannot exceed 1.0. "
            f"Got ins_rate={ins_rate}, del_rate={del_rate}"
        )


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
        "--ins-rate",
        type=validate_rate,
        default=0.15,
        help="MSI indel's insertion rate",
    )
    parser.add_argument(
        "--del-rate", type=validate_rate, default=0.15, help="MSI indels deletion rate"
    )
    parser.add_argument(
        "--boost-rate",
        type=validate_rate,
        default=0.4,
        help="Overall probability of injecting into eligible motifs",
    )
    parser.add_argument(
        "--mono-penalty",
        type=validate_mono_penalty,
        default=0.05,
        help="Probability multiplier for mononucleotide repeats "
        "(0 < value ≤ 1, default: 0.05 = 20x less likely)",
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
    return chrom.removeprefix("chr") if chrom.startswith("chr") else chrom


def get_motif_bias_multiplier(motif_length, mono_penalty):
    """
    Return bias multiplier based on motif length.
    Heavily penalizes mononucleotide repeats, keeps others normal.
    """
    if motif_length == 1:
        return mono_penalty

    return 1.0


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
                print(
                    f"[INJ-INDEL WARNING] Skipping malformed line: {line.strip()} ({e})"
                )
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
    header.info.add(
        "MSI_CHANGE", 1, "Integer", "Number of motif units inserted or deleted"
    )


def log_injection_start(ins_rate, del_rate, boost_rate, mono_penalty):
    """
    Logs the start of the injection process and the effective rates.
    """
    print("[INJ-INDEL INFO] Starting injection process...")
    print(
        f"[INJ-INDEL INFO] Effective injection rate: {(ins_rate + del_rate) * boost_rate:.2%}"
    )
    print(f"[INJ-INDEL INFO] Insertion rate: {ins_rate * boost_rate:.2%}")
    print(f"[INJ-INDEL INFO] Deletion rate: {del_rate * boost_rate:.2%}")
    print(
        f"[INJ-INDEL INFO] Mononucleotide penalty: "
        f"{mono_penalty:.3f} ({1/mono_penalty:.0f}x less likely)"
    )
    print(
        "[INJ-INDEL INFO] Using MSI tumor genotype distribution for human diploid cells:"
    )
    print("  Heterozygous (0/1): 45% - One allele affected by MSI")
    print("  Homozygous variant (1/1): 45% - Both alleles affected")
    print("  Reference (0/0): 10% - Some cells unaffected\n")
    print("[INJ-INDEL INFO] Phase 1: Processing existing variants!")


def log_phase2_start(total_existing):
    """
    Simple Phase 1 to 2 transition log.
    """
    print(
        f"[INJ-INDEL INFO] Phase 1 complete. Total existing variants: {total_existing}\n"
    )
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
    total,
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
    print(
        "  - Rate parameters too low (try increasing --ins-rate, --del-rate, or --boost-rate)"
    )
    print("  - All microsatellite positions already have variants")
    print("  - Microsatellite regions don't match reference sequence")


def assign_msi_genotype():
    """
    Assign realistic genotype for MSI variants in human tumor samples.
    """
    r = random.random()

    if r < 0.45:
        return (0, 1)
    if r < 0.90:
        return (1, 1)
    return (0, 0)


def ensure_fasta_index(ref_fasta_path):
    """
    Create FASTA index if it doesn't exist using pysam.
    """
    fai_path = ref_fasta_path + ".fai"

    if not os.path.exists(fai_path):
        print(f"[INJ-INDEL INFO] Creating FASTA index: {fai_path}")
        pysam.faidx(ref_fasta_path)
        print("[INJ-INDEL INFO] FASTA index created successfully")
    else:
        print(f"[INJ-INDEL INFO] Using existing FASTA index: {fai_path}")


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
    mono_penalty = config["mono_penalty"]

    if not os.path.exists(vcf_in):
        raise FileNotFoundError(f"[INJ-INDEL ERROR] VCF input not found: {vcf_in}")

    ensure_fasta_index(ref_fasta_path)

    random.seed(seed)

    vcf = pysam.VariantFile(vcf_in)  # pylint: disable=no-member
    add_msi_info_fields(vcf.header)

    ref_fasta = pysam.FastaFile(ref_fasta_path)  # pylint: disable=no-member

    outdir = os.path.dirname(vcf_out)
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    # pylint: disable=no-member
    out = pysam.VariantFile(vcf_out, "w", header=vcf.header)

    injected = 0
    total = 0
    existing_variant_positions = set()
    existing_variant_ranges = defaultdict(list)

    log_injection_start(ins_rate, del_rate, boost_rate, mono_penalty)

    # Phase 1: Write all existing variants and track their positions/ranges
    for rec in vcf:
        total += 1

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

                # Skip if motif is too long or less than 5 repeats or more than 50 repeats
                if motif_len > 6 or total_repeats < 5 or total_repeats > 50:
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
                except (OSError, ValueError, KeyError) as e:
                    print(
                        f"[INJ-INDEL WARNING]: Failed to fetch reference at "
                        f"{chrom}:{repeat_span_start}-{repeat_span_end}: {e}"
                    )
                    continue

                expected_seq = motif * total_repeats
                if ref_seq != expected_seq:
                    msi_skipped_validation += 1
                    if verbose and msi_skipped_validation:
                        print(f"\n{'='*80}")
                        print(
                            f"VALIDATION FAIL at {chrom}:{repeat_span_start}-{repeat_span_end}"
                        )
                        print(
                            f"    Expected: "
                            f"{expected_seq[:50]}{'...' if len(expected_seq) > 50 else ''}"
                        )
                        print(
                            f"    Found:    {ref_seq[:50]}{'...' if len(ref_seq) > 50 else ''}"
                        )
                        print(f"    Motif: {motif} x {total_repeats}")
                        print(f"\n{'='*80}")
                    continue

                # Decide whether to inject an indel based on rates
                r = random.random()
                motif_bias = get_motif_bias_multiplier(motif_len, mono_penalty)
                effective_rate = (ins_rate + del_rate) * boost_rate * motif_bias

                if r < effective_rate:
                    # if r < (ins_rate + del_rate) * boost_rate:
                    indel_type = "INS" if r < ins_rate * boost_rate else "DEL"

                    # Determine the size of the change
                    if indel_type == "INS":
                        max_ins = min(6, 50 - total_repeats)
                        if max_ins < 1:
                            continue
                        change = random.randint(1, max_ins)
                    else:
                        max_del = min(6, total_repeats - 5)
                        if max_del < 1:
                            continue
                        change = random.randint(1, max_del)

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
                            f"  With context: {upstream}"
                            f"[{ms_seq[:30]}{'...' if len(ms_seq) > 30 else ''}]{downstream}"
                        )

                        log_lines.append("\nInjected Variant:")
                        log_lines.append(f"  Type: {indel_type}")
                        log_lines.append(
                            f"  Change: {'+' if indel_type == 'INS' else '-'}{change} {motif} units"
                        )
                        log_lines.append(
                            f"  Sequence {'added' if indel_type == 'INS' else 'removed'}:"
                            f" {motif_block}"
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

                        for sample in vcf.header.samples:
                            new_rec.samples[sample]["GT"] = assign_msi_genotype()

                        if verbose:
                            log_lines.append("\nVCF Record:")
                            log_lines.append(f"  CHROM: {new_rec.chrom}")
                            log_lines.append(f"  POS: {new_rec.pos}")
                            log_lines.append(f"  REF: {new_rec.ref}")
                            log_lines.append(
                                f"  ALT: {new_rec.alts[0][:20]}"
                                f"{'...' if len(new_rec.alts[0]) > 20 else ''}"
                            )
                            log_lines.append(
                                f"  GT: {new_rec.samples[list(vcf.header.samples)[0]]['GT']}"
                            )
                            log_lines.append(
                                f"  Interpretation: Insert "
                                f"'{motif_block}' after microsatellite at position {end}"
                            )
                            log_lines.append(
                                f"  New repeat count: {total_repeats + change}"
                            )
                    else:
                        anchor_pos = end - change * motif_len
                        ref_start = anchor_pos - 1
                        ref_end = end

                        if ref_start < 0:
                            continue

                        ref_seq = ref_fasta.fetch(chrom, ref_start, ref_end).upper()
                        expected_len = 1 + (change * motif_len)
                        if len(ref_seq) != expected_len:
                            print(
                                f"WARNING: Reference length mismatch for "
                                f"deletion at {chrom}:{ref_start}-{ref_end}"
                            )
                            continue

                        new_rec.pos = ref_start + 1
                        new_rec.stop = new_rec.pos + len(ref_seq) - 1
                        new_rec.ref = ref_seq
                        new_rec.alts = (ref_seq[0],)

                        for sample in vcf.header.samples:
                            new_rec.samples[sample]["GT"] = assign_msi_genotype()

                        if verbose:
                            log_lines.append("\nVCF Record:")
                            log_lines.append(f"  CHROM: {new_rec.chrom}")
                            log_lines.append(f"  POS: {new_rec.pos}")
                            log_lines.append(
                                f"  REF: {new_rec.ref[:20]}{'...' if len(new_rec.ref) > 20 else ''}"
                            )
                            log_lines.append(f"  ALT: {new_rec.alts[0]}")
                            log_lines.append(
                                f"  GT: {new_rec.samples[list(vcf.header.samples)[0]]['GT']}"
                            )
                            log_lines.append(
                                f"  Interpretation: Delete '{motif_block}' "
                                f"from end of microsatellite region (starting at {anchor_pos + 1})"
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
    # sort_vcf_with_bcftools(vcf_out)
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

    try:
        validate_combined_rates(args.ins_rate, args.del_rate)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    repeats = load_repeats(args.bed)
    config = {
        "ins_rate": args.ins_rate,
        "del_rate": args.del_rate,
        "seed": args.seed,
        "boost_rate": args.boost_rate,
        "mono_penalty": args.mono_penalty,
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
