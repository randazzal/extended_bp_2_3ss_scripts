import sys
import glob
from collections import defaultdict
import os
import math

def parse_dp(path):
    """Parse an RNAstructure .dp file with -log10(prob) and return summed pairing probabilities per position."""
    position_probs = defaultdict(float)

    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            try:
                pos1 = int(parts[0])
                pos2 = int(parts[1])
                neglog10 = float(parts[2])
                prob = (10 ** (-neglog10)) * 100  # convert back to probability (%)
            except ValueError:
                continue

            # Skip self-pairs if present (pos1 == pos2)
            if pos1 == pos2:
                continue

            # Add probability for both positions
            position_probs[pos1] += prob
            position_probs[pos2] += prob

    return position_probs

def main():
    dp_files = glob.glob("*.dp")
    if not dp_files:
        print("No .dp files found in current directory.")
        sys.exit(1)

    for dp in dp_files:
        probs = parse_dp(dp)

        # Get filename (with extension). Use os.path.splitext(dp)[0] if you want without .dp
        filename = os.path.splitext(os.path.basename(dp))[0]

        outname = os.path.splitext(dp)[0] + "_sumProbs.tsv"
        with open(outname, "w") as out:
            out.write("Position\tSumPairingProbability(%)\tDP_File\n")
            for pos in sorted(probs.keys()):
                out.write(f"{pos}\t{probs[pos]:.6f}\t{filename}\n")

        print(f"Processed {dp} -> {outname}")

if __name__ == "__main__":
    main()
