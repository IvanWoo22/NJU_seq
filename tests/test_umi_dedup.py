"""Boundary checks for molecular identity and count compatibility."""

import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("nju_umi_dedup", ROOT / "mrna_analysis/umi_dedup.py")
umi = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = umi
SPEC.loader.exec_module(umi)


class DedupTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.annotation = self.directory / "exon.info"
        lines = []
        models = {
            "AT1G00010.1": ("+", [(101, 110), (201, 220)]),
            "AT1G00010.2": ("+", [(51, 60), (101, 110), (201, 220)]),
            "AT1G00010.3": ("+", [(101, 112), (203, 220)]),
            "AT1G00020.1": ("-", [(101, 110), (201, 220)]),
        }
        for transcript, (strand, exons) in models.items():
            for start, end in exons:
                lines.append(f"chr1\t{start}\t{end}\t{strand}\tParent=transcript:{transcript};")
        self.annotation.write_text("\n".join(lines) + "\n")
        self.models = umi.load_annotation(self.annotation)

    def row(self, name, tag="AAAAAAAA", transcript="AT1G00010.1", start=8, cigar="8=", sequence="ACGTACGT"):
        return f"{name}_{tag}\t{transcript}\t{start}\t{cigar}\t{sequence}\n"

    def collect(self, rows):
        path = self.directory / "reads.tmp"
        path.write_text("".join(rows))
        return umi.collect(path, self.models)

    def run_rows(self, rows, method="unique"):
        buckets, ambiguous, summary = self.collect(rows)
        clusterer = None
        if method == "directional":
            from umi_tools import UMIClusterer
            clusterer = UMIClusterer(cluster_method="directional")
        output = self.directory / "out.tmp"
        stats = umi.write_result(output, self.directory / "families.gz", buckets, ambiguous,
                                 method, "sample", clusterer)
        return output.read_text(), stats

    def test_coordinate_mapping_against_explicit_bases(self):
        for model in self.models.values():
            bases = [p for a, b in model.exons for p in range(a, b + 1)]
            if model.strand == "-":
                bases.reverse()
            for start in range(1, len(bases) + 1):
                self.assertEqual(model.position(start), bases[start - 1])
                for length in range(1, len(bases) - start + 2):
                    observed = [p for a, b in model.blocks(start, length) for p in range(a, b + 1)]
                    self.assertEqual(observed, sorted(bases[start - 1:start - 1 + length]))

    def test_exact_duplicates_across_isoforms_are_one_molecule(self):
        output, stats = self.run_rows([
            self.row("z"), self.row("a", transcript="AT1G00010.2", start=18),
            self.row("b", tag="CCCCCCCC"),
        ])
        self.assertEqual(stats["output_molecules"], 2)
        self.assertEqual(stats["multi_transcript_families"], 1)
        self.assertIn("a_AAAAAAAA", output)
        self.assertNotIn("z_AAAAAAAA", output)

    def test_different_endpoints_strands_and_splice_paths_are_retained(self):
        rows = [
            self.row("original"),
            self.row("different_start", start=7, cigar="9=", sequence="AACGTACGT"),
            self.row("different_end", cigar="9=", sequence="ACGTACGTA"),
            self.row("other_strand", transcript="AT1G00020.1", start=16),
            self.row("other_path", transcript="AT1G00010.3"),
        ]
        output, stats = self.run_rows(rows)
        self.assertEqual(stats["output_molecules"], 5)
        self.assertEqual(len(output.splitlines()), 5)

    def test_ambiguous_N_umis_are_not_collapsed(self):
        _, stats = self.run_rows([self.row("a", "AAAANAAA"), self.row("b", "AAAANAAA")])
        self.assertEqual(stats["output_molecules"], 2)
        self.assertEqual(stats["ambiguous_umi_reads_kept"], 2)

    def test_directional_correction_and_nonmerge_for_balanced_counts(self):
        rows = [self.row(f"a{i}") for i in range(4)] + [self.row("b", "AAAAAAAT")]
        _, exact = self.run_rows(rows)
        output, directional = self.run_rows(rows, "directional")
        self.assertEqual(exact["output_molecules"], 2)
        self.assertEqual(directional["output_molecules"], 1)
        self.assertEqual(directional["family_size_histogram"], {5: 1})
        self.assertIn("a0_AAAAAAAA", output)
        rows += [self.row(f"b{i}", "AAAAAAAT") for i in range(3)]
        _, balanced = self.run_rows(rows, "directional")
        self.assertEqual(balanced["output_molecules"], 2)

    def test_representative_uses_fewest_error_bases(self):
        output, _ = self.run_rows([self.row("a", cigar="7=1X"), self.row("z")])
        self.assertIn("z_AAAAAAAA", output)
        self.assertNotIn("a_AAAAAAAA", output)

    def test_input_order_does_not_change_output(self):
        rows = [self.row(f"a{i}") for i in range(4)] + [
            self.row("b", "AAAAAAAT"), self.row("c", "CCCCCCCC"), self.row("n", "NAAAAAAA")]
        for method in ("unique", "directional"):
            first, _ = self.run_rows(rows, method)
            second, _ = self.run_rows(list(reversed(rows)), method)
            self.assertEqual(first, second)

    def test_malformed_inputs_fail(self):
        cases = [
            [self.row("a", "AAA")],
            [self.row("a", "AAAAAAAZ")],
            [self.row("a"), self.row("a")],
            [self.row("a", transcript="AT9G99999.1")],
            [self.row("a", start=28)],
            [self.row("a", cigar="8M")],
            [self.row("a", cigar="7=")],
            ["not five columns\n"],
        ]
        for rows in cases:
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                self.collect(rows)

    def test_indel_reference_and_query_lengths(self):
        self.assertEqual(umi.cigar_lengths("3=2I4=1X2D"), (10, 10, 5))
        # An insertion with the same reference span remains in the same context.
        _, stats = self.run_rows([
            self.row("insert", cigar="4=1I4=", sequence="ACGTAACGT"), self.row("exact")])
        self.assertEqual(stats["output_molecules"], 1)

    def test_empty_input_has_explicit_zero_summary(self):
        output, stats = self.run_rows([])
        self.assertEqual(output, "")
        self.assertEqual(stats["output_molecules"], 0)

    def test_annotation_unions_duplicate_and_overlapping_exons(self):
        with self.annotation.open("a") as handle:
            handle.write("chr1\t105\t110\t+\tParent=transcript:AT1G00010.1;\n")
        model = umi.load_annotation(self.annotation)["AT1G00010.1"]
        self.assertEqual(model.length, 30)
        self.assertEqual(model.blocks(8, 8), ((108, 110), (201, 205)))

    def test_perl_annotation_parser_does_not_reuse_previous_transcript(self):
        lines = []
        for transcript in ("AT1G00010.1", "AT1G00010.2"):
            for a, b in ((101, 110), (201, 210)):
                lines.append(f"1\t{a}\t{b}\t-\tParent=transcript:{transcript};\n")
            if transcript.endswith(".1"):
                lines.append("1\t401\t500\t-\tParent=transcript:at1g_other;\n")
        annotation = "".join(lines)
        self.annotation.write_text(annotation)
        perl = os.environ.get("NJU_PERL", "perl")
        counts = self.directory / "count.tmp"
        counts.write_text("AT1G00010.1\t5\tA\tC\tG\t0\t2\n")
        result = subprocess.run([
            perl, str(ROOT / "mrna_analysis/merge.pl"), "--refstr", "Parent=transcript:",
            "--transid", "AT", "--geneid", "AT", "-i", str(counts), "--stdout",
        ], input=annotation, text=True, capture_output=True, check=True)
        self.assertEqual(result.stdout, "1\t206\t-\tA\tC\tG\tAT1G00010\t0\t2\n")
        self.assertEqual(result.stderr, "")
        reads = self.row("same", start=5) + self.row("same", transcript="AT1G00010.2", start=5)
        result = subprocess.run([
            perl, str(ROOT / "mrna_analysis/dedup.pl"), "--refstr", "Parent=transcript:",
            "--transid", "AT", "--info", str(self.annotation),
        ], input=reads, text=True, capture_output=True, check=True)
        self.assertEqual(len(result.stdout.splitlines()), 1)
        self.assertEqual(result.stderr, "")


if __name__ == "__main__":
    unittest.main()
