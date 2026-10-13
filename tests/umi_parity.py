#!/usr/bin/env python3
"""Compare Rust and Python exact-UMI counts on seeded synthetic alignments."""

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SUMMARY_FIELDS = (
    "input_reads", "output_molecules", "removed_reads", "exact_acgt_families",
    "multi_transcript_families", "mixed_insert_sequence_families",
    "umi_corrected_families", "acgt_fragment_contexts", "ambiguous_umi_reads_kept",
    "exact_molecules_including_ambiguous", "additional_molecules_collapsed_by_correction",
    "max_exact_umis_in_one_context", "family_size_histogram",
    "observed_umis_per_family_histogram", "umi_base_counts_by_position",
    "top_umis", "observed_umi_strings",
)


def write_inputs(folder, reads, seed):
    rng = random.Random(seed)
    models = {
        "AT1G00010.1": ("chr1", "+", [(101, 110), (201, 220)]),
        "AT1G00010.2": ("chr1", "+", [(51, 60), (101, 110), (201, 220)]),
        "AT1G00010.3": ("chr1", "+", [(101, 112), (203, 220)]),
        "AT1G00020.1": ("chr1", "-", [(101, 110), (201, 220)]),
        "AT1G00020.2": ("chr1", "-", [(51, 60), (101, 110), (201, 220)]),
        "AT2G00010.1": ("chr2", "+", [(101, 110), (201, 220)]),
    }
    exons = folder / "exon.info"
    exons.write_text("".join(
        f"{chrom}\t{a}\t{b}\t{strand}\tParent=transcript:{tx};\n"
        for tx, (chrom, strand, parts) in models.items() for a, b in parts
    ))
    umis = ["AAAAAAAA", "AAAAAAAT", "CCCCCCCC", "GGGGGGGG", "ACGTACGT", "AAAANAAA", "NNNNNNNN"]
    patterns = [("8=", 8), ("3=1X4=", 8), ("3=1I5=", 9), ("3=1D4=", 7), ("2=2I2=2D2=", 8)]
    rows = []
    for number in range(reads):
        tx = rng.choice(list(models))
        length = sum(b - a + 1 for a, b in models[tx][2])
        start = rng.randrange(1, length - 8 + 2)
        cigar, query = rng.choice(patterns)
        sequence = "".join(rng.choices("ACGTN", k=query))
        umi = rng.choice(umis)
        rows.append(f"read_{number:08d}_{umi}\t{tx}\t{start}\t{cigar}\t{sequence}\n")
    rng.shuffle(rows)
    alignments = folder / "reads.tsv"
    alignments.write_text("".join(rows))
    return exons, alignments


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path,
                        default=Path(os.environ.get("NJU_NATIVE", ROOT / "target/release/nju-seq")))
    parser.add_argument("--out", type=Path, help="New result directory; must not already exist")
    parser.add_argument("--reads", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260926)
    args = parser.parse_args()
    if args.reads < 0:
        parser.error("--reads must be nonnegative")
    if args.out is None:
        temporary = tempfile.TemporaryDirectory(prefix="nju-umi-parity-")
        folder = Path(temporary.name)
    else:
        folder = args.out.resolve()
        folder.mkdir(parents=True, exist_ok=False)
    exons, alignments = write_inputs(folder, args.reads, args.seed)
    commands = {
        "python": [
            sys.executable, str(ROOT / "mrna_analysis/umi_dedup.py"),
            "--input", str(alignments), "--exon-info", str(exons), "--sample", "test",
            "--method", "unique", "--output", str(folder / "python.tsv"),
            "--families", str(folder / "python.families.tsv.gz"),
            "--summary", str(folder / "python.json"),
        ],
        "rust": [
            str(args.binary.resolve()), "umi", "--input", str(alignments), "--exons", str(exons),
            "--sample", "test", "--output", str(folder / "rust.tsv"),
            "--families", str(folder / "rust.families.tsv.gz"),
            "--summary", str(folder / "rust.json"),
        ],
    }
    for name, command in commands.items():
        with (folder / f"{name}.log").open("w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    expected = json.loads((folder / "python.json").read_text())
    actual = json.loads((folder / "rust.json").read_text())
    python_rows = (folder / "python.tsv").read_bytes()
    rust_rows = (folder / "rust.tsv").read_bytes()
    python_families = gzip.decompress((folder / "python.families.tsv.gz").read_bytes())
    rust_families = gzip.decompress((folder / "rust.families.tsv.gz").read_bytes())
    assert python_rows == rust_rows, "Retained alignment TSV differs"
    assert python_families == rust_families, "Family TSV differs"
    for key in SUMMARY_FIELDS:
        assert expected[key] == actual[key], (key, expected[key], actual[key])
    family_rows = [row.split(b"\t") for row in rust_families.splitlines()[1:]]
    assert sum(int(row[7]) for row in family_rows) == args.reads
    assert len(family_rows) == actual["output_molecules"]
    report = {
        "seed": args.seed, "reads": args.reads, "output_molecules": actual["output_molecules"],
        "retained_tsv_byte_equal": True, "families_tsv_byte_equal": True,
        "all_summary_statistics_equal": True, "independent_family_count_conservation": True,
        "compared_fields": SUMMARY_FIELDS,
        "retained_tsv_sha256": hashlib.sha256(rust_rows).hexdigest(),
        "family_tsv_sha256": hashlib.sha256(rust_families).hexdigest(),
        "result_directory": str(folder), "commands": commands,
    }
    (folder / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
