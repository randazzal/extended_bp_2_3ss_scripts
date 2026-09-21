#!/usr/bin/env python3

import gzip
import argparse
import os


# ============================================================
# BARCODE DEFINITIONS
# Replace these with the barcodes from this experiment.
#
# Format:
#     "barcode_sequence": "sample_name"
# ============================================================


BARCODE_MAP = {
    "TACACCGAA": "AGAP1_bpmut50",
    "ACCTGAACC": "AGAP1_bpmut51",
    "GTTAGGTTG": "AGAP1_mini52",
    "AACTCCAAC": "AGAP1_mini53",
    "GGAAGAGCA": "AGAP1_Rup41",
    "TGCTGACAT": "AGAP1_Rup42",
    "CGTCGTAAG": "AGAP1_Rmut43",
    "TGCTAGCTC": "AGAP1_Rmut44",
    "CCTGAAGAT": "AGAP1_struct54",
    "CGATATGCC": "AGAP1_struct55",
    "AATCATGCC": "AGAP1_Vup46",
    "TCGTGACGA": "AGAP1_Vup47",
    "GTGCAACGT": "AGAP1_Vmut48",
    "GACTCACAT": "AGAP1_Vmut49",
    "TTGCCAGTC": "AGAP1_wt37",
    "ATCGCAATC": "AGAP1_wt38",

    "AGGAAGAGT": "ASAP2_bpmut64",
    "CGATGGATA": "ASAP2_bpmut65",
    "GAGCCTGTA": "ASAP2_mini66",
    "TCAGGATAC": "ASAP2_mini67",
    "GCCAAGTCA": "ASAP2_Rup56",
    "GAAGAGCGT": "ASAP2_Rup57",
    "ACCGCACAT": "ASAP2_Rmut58",
    "CACCGCTTA": "ASAP2_Rmut59",
    "ACCATGTGG": "ASAP2_struct68",
    "TACTGAAGG": "ASAP2_struct69",
    "CTCAAGACC": "ASAP2_Vup60",
    "CACATAACC": "ASAP2_Vup61",
    "ATCTACCGT": "ASAP2_Vmut62",
    "GGACCTATG": "ASAP2_Vmut63",
    "TCTACGTCC": "ASAP2_wt39",
    "CAGCAGACA": "ASAP2_wt40",

    "TGTTCCATC": "RALGAPA2_bpmut72",
    "GACGGATAA": "RALGAPA2_bpmut73",
    "TACCTGTTG": "RALGAPA2_wt35",
    "AAGCGCATG": "RALGAPA2_wt36",
    "ATGGTCCAT": "VPS13B_bpmut70",
    "GTGTCAGTT": "VPS13B_bpmut71",
    "TCAACCGCT": "VPS13B_wt33",
    "GTCGATAGA": "VPS13B_wt34"
}


BARCODE_LEN = 9

# Python slice:
#     sequence[-126:-117]
#
BARCODE_START_FROM_END = 126
BARCODE_END_FROM_END = 117

MAX_MISMATCHES = 2


def open_fastq(filename, mode="rt"):
    """Open regular or gzipped FASTQ file."""
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    return open(filename, mode)


def read_fastq(handle):
    """Yield FASTQ records as four-line lists."""
    while True:
        header = handle.readline()

        if not header:
            break

        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()

        if not quality:
            raise RuntimeError("Incomplete FASTQ record encountered.")

        yield [header, sequence, plus, quality]


def hamming_distance(seq1, seq2):
    """Return number of mismatched positions."""
    return sum(a != b for a, b in zip(seq1, seq2))


def extract_barcode(sequence):
    """
    Extract the 9-bp barcode from R2.

    Barcode position:
        29 bases from the 3' end through 21 bases from the 3' end.

    Equivalent Python slice:
        sequence[-29:-20]
    """

    sequence = sequence.strip()

    if len(sequence) < BARCODE_START_FROM_END:
        return None

    barcode = sequence[-BARCODE_START_FROM_END:-BARCODE_END_FROM_END]

    if len(barcode) != BARCODE_LEN:
        return None

    return barcode


def assign_barcode(sequence):
    """
    Assign a read to a barcode.

    Returns:
        sample_name, mismatch_count

    or:
        None, None

    Reads with multiple equally good barcode matches are
    considered ambiguous and are returned as unassigned.
    """

    observed_barcode = extract_barcode(sequence)

    if observed_barcode is None:
        return None, None

    matches = []

    for known_barcode, sample in BARCODE_MAP.items():

        distance = hamming_distance(
            observed_barcode,
            known_barcode
        )

        if distance <= MAX_MISMATCHES:
            matches.append(
                (distance, sample, known_barcode)
            )

    # No barcode within the allowed mismatch threshold
    if not matches:
        return None, None

    # Sort by mismatch count
    matches.sort(key=lambda x: x[0])

    best_distance = matches[0][0]

    # Keep only the best matches
    best_matches = [
        match for match in matches
        if match[0] == best_distance
    ]

    # More than one barcode has the same best distance
    if len(best_matches) > 1:
        return None, None

    _, sample, _ = best_matches[0]

    return sample, best_distance


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Split paired-end FASTQ reads according to a 9-bp "
            "barcode located 29-21 bases from the 3' end of R2."
        )
    )

    parser.add_argument(
        "-1",
        "--r1",
        required=True,
        help="Input R1 FASTQ(.gz)"
    )

    parser.add_argument(
        "-2",
        "--r2",
        required=True,
        help="Input R2 FASTQ(.gz)"
    )

    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Output directory"
    )

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # --------------------------------------------------------
    # Create output files
    # --------------------------------------------------------

    outputs = {}

    for sample in BARCODE_MAP.values():

        r1_out = gzip.open(
            os.path.join(
                args.outdir,
                f"{sample}_R1.fastq.gz"
            ),
            "wt"
        )

        r2_out = gzip.open(
            os.path.join(
                args.outdir,
                f"{sample}_R2.fastq.gz"
            ),
            "wt"
        )

        outputs[sample] = (r1_out, r2_out)

    # Unassigned reads
    outputs["unassigned"] = (
        gzip.open(
            os.path.join(
                args.outdir,
                "unassigned_R1.fastq.gz"
            ),
            "wt"
        ),
        gzip.open(
            os.path.join(
                args.outdir,
                "unassigned_R2.fastq.gz"
            ),
            "wt"
        )
    )

    # --------------------------------------------------------
    # Counters
    # --------------------------------------------------------

    counts = {
        sample: 0
        for sample in BARCODE_MAP.values()
    }

    counts["unassigned"] = 0

    mismatch_counts = {
        sample: {
            0: 0,
            1: 0,
            2: 0
        }
        for sample in BARCODE_MAP.values()
    }

    total = 0

    # --------------------------------------------------------
    # Process paired FASTQ files
    # --------------------------------------------------------

    with open_fastq(args.r1, "rt") as r1_handle, \
         open_fastq(args.r2, "rt") as r2_handle:

        r1_iter = read_fastq(r1_handle)
        r2_iter = read_fastq(r2_handle)

        while True:

            try:
                r1_record = next(r1_iter)
            except StopIteration:
                r1_record = None

            try:
                r2_record = next(r2_iter)
            except StopIteration:
                r2_record = None

            # Both files finished
            if r1_record is None and r2_record is None:
                break

            # One file finished before the other
            if r1_record is None or r2_record is None:
                raise RuntimeError(
                    "R1 and R2 contain different numbers of reads."
                )

            total += 1

            r2_sequence = r2_record[1].strip()

            sample, distance = assign_barcode(
                r2_sequence
            )

            # ------------------------------------------------
            # Unassigned read
            # ------------------------------------------------

            if sample is None:

                sample = "unassigned"

            # ------------------------------------------------
            # Write paired reads
            # ------------------------------------------------

            r1_out, r2_out = outputs[sample]

            r1_out.writelines(r1_record)
            r2_out.writelines(r2_record)

            counts[sample] += 1

            # Record mismatch distribution
            if sample != "unassigned" and distance is not None:
                mismatch_counts[sample][distance] += 1

    # --------------------------------------------------------
    # Close output files
    # --------------------------------------------------------

    for r1_out, r2_out in outputs.values():
        r1_out.close()
        r2_out.close()

    # --------------------------------------------------------
    # Print summary
    # --------------------------------------------------------

    print()
    print("=" * 60)
    print("Barcode demultiplexing complete")
    print("=" * 60)

    print(f"Total read pairs: {total:,}")
    print()

    for sample in BARCODE_MAP.values():

        n = counts[sample]

        if total > 0:
            pct = 100 * n / total
        else:
            pct = 0

        print(
            f"{sample:20s} "
            f"{n:12,d} "
            f"({pct:6.2f}%)"
        )

        if n > 0:
            print(
                f"    0 mismatches: {mismatch_counts[sample][0]:,}"
            )
            print(
                f"    1 mismatch:  {mismatch_counts[sample][1]:,}"
            )
            print(
                f"    2 mismatches: {mismatch_counts[sample][2]:,}"
            )

    n = counts["unassigned"]

    if total > 0:
        pct = 100 * n / total
    else:
        pct = 0

    print(
        f"{'unassigned':20s} "
        f"{n:12,d} "
        f"({pct:6.2f}%)"
    )

    print("=" * 60)


if __name__ == "__main__":
    main()
