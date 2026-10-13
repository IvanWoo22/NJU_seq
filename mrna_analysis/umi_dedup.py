#!/usr/bin/env python3
"""Deduplicate one NJU-seq library after almostuniquematch, before count.pl.

Input: five columns (read ID, transcript, 1-based start, =/X/I/D CIGAR, sequence).
Each input read ID must occur exactly once. A molecular context is the chromosome,
strand and complete exonic genomic path spanned by R1, including both endpoints.
This collapses equivalent transcript representations without merging stop sites.
UMIs containing N are retained individually; malformed UMIs/coordinates fail.
"""

import argparse
import bisect
from collections import Counter
from dataclasses import dataclass
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import sys
import tempfile
import time


@dataclass(frozen=True)
class Transcript:
    chrom: str
    strand: str
    exons: tuple
    offsets: tuple
    length: int

    def blocks(self, start, length):
        """Return ascending genomic, inclusive intervals of a transcript range."""
        if start < 1 or length < 1 or start + length - 1 > self.length:
            raise ValueError(f"range {start}+{length} exceeds transcript length {self.length}")
        lo, hi = start - 1, start - 1 + length
        if self.strand == "-":
            lo, hi = self.length - hi, self.length - lo
        result = []
        i = bisect.bisect_right(self.offsets, lo) - 1
        while lo < hi:
            a, b = self.exons[i]
            take = min(hi - lo, b - a + 1 - (lo - self.offsets[i]))
            g_start = a + lo - self.offsets[i]
            result.append((g_start, g_start + take - 1))
            lo += take
            i += 1
        return tuple(result)

    def position(self, position):
        return self.blocks(position, 1)[0][0]


def load_annotation(path, refstr="Parent=transcript:", transid="AT"):
    pattern = re.compile(re.escape(refstr) + "(" + re.escape(transid) + r"\w+\.[0-9]+)")
    raw = {}
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5:
                raise ValueError(f"exon.info line {number}: expected five fields")
            chrom, a, b, strand, attributes = fields
            match = pattern.search(attributes)
            if not match:
                continue
            transcript = match[1]
            a, b = int(a), int(b)
            if a < 1 or b < a or strand not in ("+", "-"):
                raise ValueError(f"exon.info line {number}: invalid exon")
            # Match the existing merge.pl chromosome normalization.
            chrom = chrom.replace("chr", "", 1)
            if transcript not in raw:
                raw[transcript] = (chrom, strand, [])
            if raw[transcript][:2] != (chrom, strand):
                raise ValueError(f"conflicting annotation for {transcript}")
            raw[transcript][2].append((a, b))
    result = {}
    for transcript, (chrom, strand, ranges) in raw.items():
        merged = []
        for a, b in sorted(ranges):
            if merged and a <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(b, merged[-1][1]))
            else:
                merged.append((a, b))
        offsets, length = [], 0
        for a, b in merged:
            offsets.append(length)
            length += b - a + 1
        result[transcript] = Transcript(chrom, strand, tuple(merged), tuple(offsets), length)
    if not result:
        raise ValueError("no transcripts matched the annotation options")
    return result


def cigar_lengths(cigar):
    ops = re.findall(r"([1-9][0-9]*)([=XID])", cigar)
    if not ops or "".join(n + op for n, op in ops) != cigar:
        raise ValueError(f"unsupported CIGAR {cigar!r}; require end-to-end --xeq (=,X,I,D)")
    reference = sum(int(n) for n, op in ops if op != "I")
    query = sum(int(n) for n, op in ops if op != "D")
    errors = sum(int(n) for n, op in ops if op != "=")
    return reference, query, errors


@dataclass
class UmiFamily:
    count: int
    rank: tuple
    row: str
    transcript: str
    sequence: str
    multiple_transcripts: bool = False
    multiple_sequences: bool = False

    def add(self, rank, row, transcript, sequence):
        self.count += 1
        self.multiple_transcripts |= transcript != self.transcript
        self.multiple_sequences |= sequence != self.sequence
        if rank < self.rank:
            self.rank, self.row = rank, row


def collect(input_path, annotation, umi_length=8):
    buckets, ambiguous = {}, []
    seen = set()
    umi_counts = Counter()
    base_counts = [Counter() for _ in range(umi_length)]
    digest = hashlib.sha256()
    input_reads = 0
    with open(input_path) as handle:
        for number, line in enumerate(handle, 1):
            digest.update(line.encode())
            try:
                row = line.rstrip("\n")
                fields = row.split("\t")
                if len(fields) != 5:
                    raise ValueError("expected five fields")
                read_id, reference, start, cigar, sequence = fields
                if read_id in seen:
                    raise ValueError(f"read ID repeated: {read_id}; input must be almostuniquematch output")
                seen.add(read_id)
                if "_" not in read_id:
                    raise ValueError("read ID has no _UMI suffix")
                umi = read_id.rsplit("_", 1)[1]
                if len(umi) != umi_length or re.fullmatch(r"[ACGTN]+", umi) is None:
                    raise ValueError(f"invalid {umi_length}-nt UMI {umi!r}")
                tx_match = re.match(r"([^|\s]+)", reference)
                tx = tx_match[1]
                if tx not in annotation:
                    raise ValueError(f"transcript missing from annotation: {tx}")
                ref_length, query_length, errors = cigar_lengths(cigar)
                if ref_length < 2 or len(sequence) < 3 or len(sequence) != query_length:
                    raise ValueError("CIGAR/sequence length mismatch or fragment too short")
                model = annotation[tx]
                start = int(start)
                blocks = model.blocks(start, ref_length)
                context = (model.chrom, model.strand, blocks)
                site = model.position(start + ref_length - 2)
                # Lexicographic ties provide a reproducible representative.
                rank = (errors, read_id, row)
                umi_counts[umi] += 1
                for i, base in enumerate(umi):
                    base_counts[i][base] += 1
                input_reads += 1
                if "N" in umi:
                    ambiguous.append((context, site, umi, row))
                else:
                    if context not in buckets:
                        buckets[context] = (site, {})
                    previous_site, families = buckets[context]
                    if previous_site != site:
                        raise ValueError("same fragment path maps to conflicting stop sites")
                    if umi in families:
                        families[umi].add(rank, row, tx, sequence)
                    else:
                        families[umi] = UmiFamily(1, rank, row, tx, sequence)
            except (ValueError, KeyError) as exc:
                raise ValueError(f"{input_path}:{number}: {exc}") from exc
    return buckets, ambiguous, {
        "input_reads": input_reads,
        "input_sha256": digest.hexdigest(),
        "umi_base_counts_by_position": [dict(sorted(c.items())) for c in base_counts],
        "top_umis": [{"umi": u, "reads": n} for u, n in umi_counts.most_common(20)],
        "observed_umi_strings": len(umi_counts),
    }


def cluster_umis(families, method, clusterer=None):
    if method == "unique":
        return [(umi,) for umi in sorted(families)]
    # Stable insertion order also fixes equal-count root traversal in UMI-tools.
    counts = {umi.encode(): families[umi].count
              for umi in sorted(families, key=lambda u: (-families[u].count, u))}
    groups = clusterer(counts, threshold=1)
    normalized = [tuple(sorted((u.decode() for u in g),
                              key=lambda u: (-families[u].count, u))) for g in groups if g]
    assert sorted(u for g in normalized for u in g) == sorted(families)
    return sorted(normalized)


def write_result(output, families_output, buckets, ambiguous, method, sample, clusterer=None):
    stats = Counter()
    size_histogram = Counter()
    distinct_histogram = Counter()
    maximum_umis = 0
    with open(output, "w") as out, gzip.open(families_output, "wt", compresslevel=3) as report:
        report.write("FamilyID\tSample\tChromosome\tStrand\tGenomicBlocks\tCandidatePosition"
                     "\tRootUMI\tReads\tObservedUMIs\tRepresentativeRead\tStatus\n")

        def emit(context, site, umi, row, count, observed, status):
            stats["output_molecules"] += 1
            size_histogram[count] += 1
            distinct_histogram[len(observed)] += 1
            chrom, strand, blocks = context
            block_text = ",".join(f"{a}-{b}" for a, b in blocks)
            read_id = row.split("\t", 1)[0]
            out.write(row + "\n")
            report.write(f"{sample}:{stats['output_molecules']}\t{sample}\t{chrom}\t{strand}"
                         f"\t{block_text}\t{site}\t{umi}\t{count}\t{','.join(observed)}"
                         f"\t{read_id}\t{status}\n")

        for context in sorted(buckets):
            site, families = buckets[context]
            maximum_umis = max(maximum_umis, len(families))
            stats["exact_acgt_families"] += len(families)
            for group in cluster_umis(families, method, clusterer):
                root = group[0]
                selected = families[root]
                count = sum(families[u].count for u in group)
                stats["multi_transcript_families"] += (
                    any(families[u].multiple_transcripts for u in group)
                    or len({families[u].transcript for u in group}) > 1)
                stats["mixed_insert_sequence_families"] += (
                    any(families[u].multiple_sequences for u in group)
                    or len({families[u].sequence for u in group}) > 1)
                stats["umi_corrected_families"] += len(group) > 1
                emit(context, site, root, selected.row, count, group, "acgt")
        for context, site, umi, row in sorted(ambiguous, key=lambda x: (x[0], x[3])):
            emit(context, site, umi, row, 1, (umi,), "ambiguous_umi_kept")
    # Explicit zero values make empty input and exact-only runs unambiguous.
    for key in ("output_molecules", "exact_acgt_families", "multi_transcript_families",
                "mixed_insert_sequence_families", "umi_corrected_families"):
        stats[key] += 0
    return {
        **dict(stats),
        "acgt_fragment_contexts": len(buckets),
        "ambiguous_umi_reads_kept": len(ambiguous),
        "exact_molecules_including_ambiguous": stats["exact_acgt_families"] + len(ambiguous),
        "additional_molecules_collapsed_by_correction":
            stats["exact_acgt_families"] + len(ambiguous) - stats["output_molecules"],
        "max_exact_umis_in_one_context": maximum_umis,
        "family_size_histogram": dict(sorted(size_histogram.items())),
        "observed_umis_per_family_histogram": dict(sorted(distinct_histogram.items())),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--exon-info", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--families", type=Path, required=True, help="gzip-compressed family TSV")
    parser.add_argument("--summary", type=Path, required=True, help="JSON summary")
    parser.add_argument("--sample", required=True, help="one library only; never pool NC/TR libraries")
    parser.add_argument("--method", choices=("unique", "directional"), default="unique")
    parser.add_argument("--umi-length", type=int, default=8)
    parser.add_argument("--refstr", default="Parent=transcript:")
    parser.add_argument("--transid", default="AT")
    args = parser.parse_args(argv)
    if args.umi_length < 1 or re.fullmatch(r"[A-Za-z0-9_.-]+", args.sample) is None:
        parser.error("positive UMI length and a simple, nonempty sample identifier required")
    targets = (args.output, args.families, args.summary)
    if len({p.resolve() for p in targets}) != 3:
        parser.error("output, families and summary must be different files")
    for path in targets:
        if path.exists():
            parser.error(f"refusing to overwrite {path}")
        if not path.parent.is_dir():
            parser.error(f"output directory does not exist: {path.parent}")
    start = time.monotonic()
    clusterer, umi_version = None, None
    if args.method == "directional":
        from umi_tools import UMIClusterer, __version__
        clusterer = UMIClusterer(cluster_method="directional")
        umi_version = __version__
    annotation = load_annotation(args.exon_info, args.refstr, args.transid)
    buckets, ambiguous, summary = collect(args.input, annotation, args.umi_length)
    with tempfile.TemporaryDirectory(prefix=".umi_dedup_", dir=args.output.parent) as temp:
        temp = Path(temp)
        stats = write_result(temp / "reads.tmp", temp / "families.tsv.gz", buckets, ambiguous,
                             args.method, args.sample, clusterer)
        summary.update(stats)
        summary.update({
            "sample": args.sample, "method": args.method, "umi_length": args.umi_length,
            "input": str(args.input.resolve()), "exon_info": str(args.exon_info.resolve()),
            "exon_info_sha256": hashlib.sha256(args.exon_info.read_bytes()).hexdigest(),
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "python_version": sys.version.split()[0], "umi_tools_version": umi_version,
            "key": "library + chromosome + strand + complete R1 genomic exon path + UMI",
            "ambiguous_umi_policy": "retain each N-containing read separately",
            "correction_distance": 1 if args.method == "directional" else 0,
            "removed_reads": summary["input_reads"] - stats["output_molecules"],
            "elapsed_seconds": round(time.monotonic() - start, 3),
            "peak_rss_kib_linux": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        })
        (temp / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
        for source, target in zip(("reads.tmp", "families.tsv.gz", "summary.json"), targets):
            # Link publishes complete files and fails if another process created the target.
            os.link(temp / source, target)
    print(json.dumps({k: summary[k] for k in ("sample", "method", "input_reads",
                                            "output_molecules", "removed_reads")}), file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, ImportError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
