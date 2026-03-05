from collections import defaultdict

import argparse

# usage python3 find_UGC.py --fasta input_sequences.fasta > output.txt
parser = argparse.ArgumentParser(description="Process microexon sequences.")
parser.add_argument("--fasta", required=True, help="Input multi-FASTA file")
args = parser.parse_args()

fasta_file = args.fasta

counts = {}

with open(fasta_file) as f:
    seq_id = None
    seq = []

    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if seq_id is not None:
                sequence = "".join(seq).upper()
                counts[seq_id] = sequence.count("TGC")
            seq_id = line[1:]
            seq = []
        else:
            seq.append(line)

    # last sequence
    if seq_id is not None:
        sequence = "".join(seq).upper()
        counts[seq_id] = sequence.count("TGC")

# print results
print("Sequence_ID\tTGC_count")
for k, v in counts.items():
    print(f"{k}\t{v}")
