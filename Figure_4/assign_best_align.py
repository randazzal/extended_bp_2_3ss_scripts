#!/usr/bin/env python3

import argparse
import gzip
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def norm_id(name):
    name = name.split()[0]
    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]
    return name


def read_fastq_gz(path):
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline()
            p = fh.readline()
            q = fh.readline()
            if not q:
                raise ValueError("Truncated FASTQ: %s" % path)
            h = h.rstrip("\n")
            s = s.rstrip("\n")
            p = p.rstrip("\n")
            q = q.rstrip("\n")
            if not h.startswith("@") or not p.startswith("+"):
                raise ValueError("Invalid FASTQ format in %s" % path)
            yield norm_id(h[1:]), h[1:], s, q


def strip_fasta_suffix(path):
    name = path.name
    if name.endswith(".gz"):
        name = name[:-3]
    for suf in [".fasta", ".fa", ".fna", ".fas"]:
        if name.lower().endswith(suf):
            name = name[:-len(suf)]
            break
    return name


def open_text_fasta(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "r", errors="replace")


def fasta_iter(path):
    header = None
    seq_chunks = []
    with open_text_fasta(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)


def build_combined_reference(ref_dir, combined_fasta):
    ref_dir = Path(ref_dir)
    combined_fasta = Path(combined_fasta)

    ref_files = sorted(
        [
            p for p in ref_dir.iterdir()
            if p.is_file() and (
                p.name.lower().endswith(".fa") or
                p.name.lower().endswith(".fasta") or
                p.name.lower().endswith(".fna") or
                p.name.lower().endswith(".fas") or
                p.name.lower().endswith(".fa.gz") or
                p.name.lower().endswith(".fasta.gz") or
                p.name.lower().endswith(".fna.gz") or
                p.name.lower().endswith(".fas.gz")
            )
        ]
    )

    if not ref_files:
        sys.exit("No FASTA files found in %s" % ref_dir)

    with open(str(combined_fasta), "w") as out:
        for ref_path in ref_files:
            ref_name = strip_fasta_suffix(ref_path)
            for seq_id, seq in fasta_iter(str(ref_path)):
                out.write(">%s__%s\n" % (ref_name, seq_id))
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i+60] + "\n")

    return ref_files


def ensure_bwa_index(combined_fasta):
    if os.path.exists(str(combined_fasta) + ".bwt"):
        return
    cmd = ["bwa", "index", str(combined_fasta)]
    ret = subprocess.call(cmd)
    if ret != 0:
        sys.exit("bwa index failed on %s" % combined_fasta)


def get_as(fields):
    for f in fields[11:]:
        if f.startswith("AS:i:"):
            return int(f[5:])
    return 0


def run_bwa_mem_and_collect_best(ref_fasta, r1, r2, threads):
    cmd = [
        "bwa", "mem",
        "-t", str(threads),
        str(ref_fasta),
        str(r1),
        str(r2),
    ]

    p = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True
    )

    # best_scores[read_id][mate][ref_name] = best AS observed
    best_scores = defaultdict(lambda: {1: defaultdict(int), 2: defaultdict(int)})

    for line in p.stdout:
        if not line or line[0] == "@":
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            continue

        qname = norm_id(fields[0])
        flag = int(fields[1])
        rname = fields[2]

        if flag & 4:
            continue
        if flag & 2048:
            continue

        if flag & 64:
            mate = 1
        elif flag & 128:
            mate = 2
        else:
            mate = 1

        if rname == "*":
            continue

        ref_group = rname.split("__", 1)[0]
        ascore = get_as(fields)

        if ascore > best_scores[qname][mate][ref_group]:
            best_scores[qname][mate][ref_group] = ascore

    stderr = p.stderr.read()
    ret = p.wait()
    if ret != 0:
        sys.stderr.write(stderr)
        sys.exit("bwa mem failed")

    return best_scores


def choose_best_ref(mate1_scores, mate2_scores, min_pair_as):
    refs = set(mate1_scores.keys()) | set(mate2_scores.keys())
    if not refs:
        return "unassigned", 0

    best_ref = "unassigned"
    best_score = -1

    for ref in refs:
        s1 = mate1_scores.get(ref, 0)
        s2 = mate2_scores.get(ref, 0)
        pair_as = s1 + s2
        if pair_as > best_score or (pair_as == best_score and ref < best_ref):
            best_score = pair_as
            best_ref = ref

    if best_score < min_pair_as:
        return "unassigned", best_score

    return best_ref, best_score


def main():
    ap = argparse.ArgumentParser(
        description="Assign paired-end reads to the best-matching reference FASTA using bwa mem."
    )
    ap.add_argument("r1", help="R1 FASTQ.gz")
    ap.add_argument("r2", help="R2 FASTQ.gz")
    ap.add_argument("ref_dir", help="Directory containing reference FASTA files")
    ap.add_argument("out_dir", help="Output directory")
    ap.add_argument("--threads", type=int, default=4, help="Threads for bwa mem")
    ap.add_argument("--min-pair-as", type=int, default=0, help="Minimum combined AS to assign")
    args = ap.parse_args()

    r1 = Path(args.r1)
    r2 = Path(args.r2)
    ref_dir = Path(args.ref_dir)
    out_dir = Path(args.out_dir)

    if not r1.exists():
        sys.exit("Missing R1: %s" % r1)
    if not r2.exists():
        sys.exit("Missing R2: %s" % r2)
    if not ref_dir.is_dir():
        sys.exit("Not a directory: %s" % ref_dir)

    out_dir.mkdir(parents=True, exist_ok=True)

    combined_fasta = out_dir / "combined_references.fasta"
    build_combined_reference(ref_dir, combined_fasta)
    ensure_bwa_index(combined_fasta)

    sys.stderr.write("Running bwa mem...\n")
    best_scores = run_bwa_mem_and_collect_best(combined_fasta, r1, r2, args.threads)

    assignments = {}
    for read_id, mate_dict in best_scores.items():
        mate1 = mate_dict[1]
        mate2 = mate_dict[2]
        best_ref, best_score = choose_best_ref(mate1, mate2, args.min_pair_as)
        assignments[read_id] = (best_ref, best_score)

    handles = {}

    def get_handle(ref_name, mate):
        key = (ref_name, mate)
        if key not in handles:
            path = out_dir / ("%s_%s.fastq.gz" % (ref_name, mate))
            handles[key] = gzip.open(str(path), "wt")
        return handles[key]

    report = open(str(out_dir / "assignments.tsv"), "w")
    report.write("read_id\tbest_reference\tpair_as\n")

    r1_iter = read_fastq_gz(r1)
    r2_iter = read_fastq_gz(r2)

    for rec1, rec2 in zip(r1_iter, r2_iter):
        id1, header1, seq1, qual1 = rec1
        id2, header2, seq2, qual2 = rec2

        if id1 != id2:
            raise ValueError("R1/R2 names do not match: %s != %s" % (id1, id2))

        best_ref, best_score = assignments.get(id1, ("unassigned", 0))
        if best_score < args.min_pair_as:
            best_ref = "unassigned"

        w1 = get_handle(best_ref, "R1")
        w2 = get_handle(best_ref, "R2")

        w1.write("@%s\n%s\n+\n%s\n" % (header1, seq1, qual1))
        w2.write("@%s\n%s\n+\n%s\n" % (header2, seq2, qual2))

        report.write("%s\t%s\t%d\n" % (id1, best_ref, best_score))

    report.close()
    for h in handles.values():
        h.close()

    sys.stderr.write("Done.\n")


if __name__ == "__main__":
    main()
