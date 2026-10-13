"""Behavioral regression tests for repaired original Perl/Shell scripts."""
import gzip
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
PERL = os.environ.get("NJU_PERL", "perl")


class ScriptRepairs(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="nju-repairs-")
        self.addCleanup(self.tmp.cleanup)
        self.d = Path(self.tmp.name)

    def put(self, name, text):
        p = self.d / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
        return p

    def run_script(self, path, args=(), data="", ok=True):
        p = subprocess.run([PERL, str(ROOT / path), *map(str, args)],
                           input=data, text=True, capture_output=True, timeout=15)
        if ok:
            self.assertEqual(p.returncode, 0, p.stderr)
        else:
            self.assertNotEqual(p.returncode, 0)
        return p

    def test_multimatch_distinct_keys_names_and_empty_input(self):
        for kind in ("mrna", "rrna"):
            script = f"{kind}_analysis/multimatch_judge.pl"
            for rows in (
                "q\tAT1G03620.1\t11\t3=\tAAA\nq\tAT1G03620.11\t1\t3=\tAAA\n",
                "q.1\tTX\t1\t3=\tAAA\nqa1\tTX\t1\t3=\tAAA\n",
            ):
                self.assertEqual(len(self.run_script(script, data=rows).stdout.splitlines()), 2)
            p = self.run_script(script)
            self.assertEqual((p.stdout, p.stderr), ("", ""))

    def test_rrna_insertions_not_classified_as_perfect(self):
        rows = "q\tTX1\t1\t3=1I3=\tAAAAAAA\nq\tTX2\t1\t6=\tAAAAAA\n"
        self.assertEqual(self.run_script("rrna_analysis/multimatch_judge.pl", data=rows).stdout,
                         rows.splitlines()[1] + "\n")

    def test_quality_rejects_unsupported_cigar(self):
        for kind in ("mrna", "rrna"):
            self.run_script(f"{kind}_analysis/matchquality_judge.pl", data="q\tTX\t1\t18M\t" + "A"*18 + "\n", ok=False)

    def test_count_indels_use_reference_positions(self):
        f = self.put("align.tsv", "q\tTX\t10\t3=1I1=\tACGTA\n")
        self.assertEqual(self.run_script("mrna_analysis/count.pl", [f]).stdout,
                         "TX\t12\tC\tG\tA\t0\t1\n")
        f.write_text("q\tTX\t10\t3=1D1=\tACGA\n")
        self.assertEqual(self.run_script("mrna_analysis/count.pl", [f]).stdout,
                         "TX\t13\tG\tN\tA\t0\t1\n")

    def test_score_header_empty_replicate_and_duplicates(self):
        nc = self.put("NC/count.tsv", "1\t100\t+\tA\tC\tG\tG1\t0\t1\n")
        tr = self.put("TR/count.tsv", "1\t100\t+\tA\tC\tG\tG1\t0\t40\n")
        # NC total must allow the unchanged fold-change criterion to pass.
        nc.write_text(nc.read_text() + "1\t101\t+\tA\tC\tG\tG1\t0\t100\n")
        tr.write_text(tr.read_text() + "1\t101\t+\tA\tC\tG\tG1\t0\t1\n")
        rows = self.run_script("mrna_analysis/score.pl", [nc, tr]).stdout.splitlines()
        self.assertEqual(len(rows), 2)
        self.assertEqual([len(r.split("\t")) for r in rows], [15, 15])
        tr2 = self.put("TR2/count.tsv", "")
        rows = self.run_script("mrna_analysis/score.pl", [nc, tr, tr2]).stdout.splitlines()
        self.assertEqual(len(rows), 1)
        self.assertEqual(len(rows[0].split("\t")), 18)
        # Joint score is an average threshold: TR2 need not pass 30 on its own.
        tr2.write_text("1\t100\t+\tA\tC\tG\tG1\t0\t25\n1\t101\t+\tA\tC\tG\tG1\t0\t1\n")
        rows = self.run_script("mrna_analysis/score.pl", [nc,tr,tr2]).stdout.splitlines()
        self.assertEqual([len(row.split("\t")) for row in rows],[18,18])
        self.assertEqual(rows[1].split("\t")[-4:],["197","98.5","63","31.5"])
        tr.write_text(tr.read_text()*2)
        self.run_script("mrna_analysis/score.pl", [nc, tr], ok=False)

    def test_coordinate_dedup_preserves_opposite_strands(self):
        annotation = self.put("exon.info", "1\t100\t102\t+\tParent=transcript:AT1G00010.1;\n1\t98\t100\t-\tParent=transcript:AT1G00020.1;\n")
        rows = "q\tAT1G00010.1\t1\t3=\tAAA\nq\tAT1G00020.1\t1\t3=\tAAA\n"
        args = ["--refstr", "Parent=transcript:", "--transid", "AT", "--info", annotation]
        self.assertEqual(self.run_script("mrna_analysis/dedup.pl", args, rows).stdout, rows)
        self.run_script("mrna_analysis/dedup.pl", args, rows.replace("\t1\t", "\t4\t"), ok=False)

    def test_set_membership_duplicates_empty_and_current_scores(self):
        row = "1\t100\t+\tA\tC\tG\tG1\t2\t30\t40\t50\t20\t20\t50\t50\n"
        a = self.put("a.tsv", row * 3)
        b = self.put("b.tsv", "")
        c = self.put("c.tsv", row.replace("100", "200"))
        self.assertEqual(self.run_script("tool/common.pl", [a,b,c]).stdout, "")
        for script in ("specific_point.pl", "specific_point_matrix.pl", "remove_common_point.pl"):
            p = self.run_script("mrna_analysis/"+script, [a,b])
            self.assertEqual(p.stdout.count("1\t100\t+"), 1)
        b.write_text(row)
        c.write_text(row)
        self.assertEqual(self.run_script("tool/common.pl", [a,b,c]).stdout,
                         "1\t100\t+\tA\tC\tG\tG1\t150\n")
        self.run_script("mrna_analysis/specific_point_withsocre.pl", [a], ok=False)

    def test_top_intersection_terminates_at_requested_number(self):
        rows = "".join(f"1\t{n}\t+\tA\tC\tG\tG1\t1\t2\t3\t4\t{score}\n" for n,score in [(100,30),(200,20),(300,10)])
        paths = [self.put(f"{i}.tsv", rows) for i in range(3)]
        p = self.run_script("mrna_analysis/extract_point_neo.pl", [*paths,1,1])
        self.assertEqual(len(p.stdout.splitlines()), 1)
        self.run_script("mrna_analysis/extract_point_neo.pl", [*paths,1,0], ok=False)

    def test_motif_bounds_chr_alias_and_last_column(self):
        fa = self.put("genome.fa", ">chr1\nAACCGGTT\n")
        sites = self.put("sites.tsv", "chr1\t1\t+\nchr1\t4\t+\n")
        p = self.run_script("mrna_analysis/motif_nm.pl", [fa,sites,2,2])
        self.assertEqual(p.stdout, "chr1\t4\t+\tACCGG\n")
        motifs = self.put("motifs.tsv", "miR\taccgg\n")
        annotated = self.put("annotated.tsv", "1\t4\t+\tA\tC\tG\tG1\t2\t30\t40\t50\t20\t20\t50\t50\tACCGG\n")
        self.assertTrue(self.run_script("mrna_analysis/motif_compare.pl", [motifs,annotated]).stdout.endswith("\t1\tmiR\n"))
        short = self.run_script("mrna_analysis/motif_mirna.pl", data=">short\nACG\n>long\nAACCGGTT\n")
        self.assertEqual(short.stdout, "long\tAACCGGUU\n")

    def test_region_coordinate_default_and_legacy_option(self):
        regions = self.put("r7.tsv", "1\tG1\t+\tTX1\t1-10\t11-30\t31-40\n")
        point = "1\t11\t+\tA\tC\tG\tG1\n"
        p = self.run_script("mrna_analysis/judge_transcriptregion_2.pl", [regions], point)
        self.assertEqual(p.stdout.split("\t")[5:9:3], ["0", "1"])
        p = self.run_script("mrna_analysis/judge_transcriptregion_2.pl", ["--legacy-shift", regions], point)
        self.assertEqual(p.stdout.split("\t")[5:9:3], ["1", "0"])
        p = self.run_script("mrna_analysis/main_transcript_4.pl", data="1\t+\t11\tcds\t0.1\t0\t20\t0\t2\t1\n")
        self.assertAlmostEqual(sum(float(r.split("\t")[2]) for r in p.stdout.splitlines()), 1)

    def test_exon_zero_distance_and_oriented_first_last(self):
        model = self.put("exons.tsv", "TX1.1\tG1\t1\t100-120,180-200\t100,120,180,200\n")
        sites = self.put("sites.tsv", "1\t100\t+\tA\tC\tG\tG1\n")
        bar, porta = self.d/"bar.tsv", self.d/"porta.tsv"
        self.run_script("mrna_analysis/exon_distance_2.pl", [model,sites,bar,porta])
        self.assertEqual(sum(float(r.split("\t")[2]) for r in bar.read_text().splitlines()), 1)
        rows = "110\tG1\tTX1\t100-120,180-200\t-\n190\tG1\tTX1\t100-120,180-200\t-\n"
        p = self.run_script("mrna_analysis/exon_distance_patch.pl", ["--which","first"], rows)
        self.assertEqual(p.stdout.split("\t")[0], "190")

    def test_tn5_terminal_tag_and_exact_boundary(self):
        fa = self.put("ref.fa", ">TX\nAAACCCCCC\n")
        sam = self.put("a.sam", "q_sub_GGGGGGGG\t0\tTX\t1\t60\t3=\t*\t0\t0\tAAA\tIII\n")
        self.assertEqual(self.run_script("mrna_analysis/tn5_check.pl", [fa,sam]).stdout,
                         "q_sub\tGGGGGGGG\tCCCCCC\t6\n")

    def test_fastq_delete_validates_and_preserves_unlisted_ids(self):
        names = self.put("names.tsv", "a\n")
        fq = self.put("reads.fq", "@a\nAC\n+\nII\n@b\nGT\n+\nII\n")
        output = self.d / "out.fq"
        self.run_script("tool/delete_fastq.pl", ["-n", names, "-i", fq, "-o", output])
        self.assertEqual(output.read_text(), "@b\nGT\n+\nII\n")
        fq.write_text("@a\nAC\n+\n")
        self.run_script("tool/delete_fastq.pl", ["-n",names,"-i",fq,"-o",output], ok=False)

    def test_reference_counts_are_literal_and_mutations_handle_uncovered(self):
        fa = self.put("ref.fa", ">TX.1\nAAAA\n")
        align = self.put("align.tsv", "q\tTXa1\t1\t3=\tAAA\n")
        p = self.run_script("rrna_analysis/readend_count.pl", [fa,align,"TX.1"])
        self.assertEqual(sum(int(r.split("\t")[3]) for r in p.stdout.splitlines()), 0)
        align.write_text("q\tTX.1\t1\t3M\tACA\n")
        p = self.run_script("rrna_analysis/mutation_count.pl", [fa,align,"TX.1"])
        rows = [r.split("\t") for r in p.stdout.splitlines()]
        self.assertEqual(rows[1][7:9], ["C","mut"])
        self.assertEqual(rows[3][7:9], ["N","uncovered"])

    def test_fetch_fasta_matches_titles_not_sequence(self):
        fa = self.put("a.fa", ">keep\nAUUU\n>other\nkeep\n")
        p = self.run_script("tool/fetch_fasta.pl", ["--file",fa,"--string","keep","--rna2dna"])
        self.assertEqual(p.stdout, ">keep\nATTT\n")

    def test_grouped_gff_flush_and_nonnadjacent_overlap(self):
        def gff(kind,a,b,attrs,strand="+"):
            return f"1\ttest\t{kind}\t{a}\t{b}\t.\t{strand}\t.\t{attrs}\n"
        gene = gff("gene",100,200,"ID=gene:G1;")
        models = "".join(gff("transcript",a,b,f"transcript_id=TX1.{i};") + gff("exon",a,b,f"transcript_id=TX1.{i};")
                         for i,a,b in [(1,100,200),(2,110,190)])
        alter, unique = self.d/"alter.tsv", self.d/"unique.tsv"
        self.run_script("mrna_analysis/judge_altersplice.pl", ["--transwording","transcript","--geneid","ID=gene:","--alter",alter,"--unique",unique], gene+models)
        self.assertTrue(alter.read_text().startswith("1\tG1\t+"))
        ranges = self.put("range.tsv", "1\tG1\t+\t100-200\t-\t-\n")
        annotation = gene+gff("mRNA",100,200,"ID=transcript:TX1.1;")+gff("CDS",120,180,"Parent=transcript:TX1.1;")
        p = self.run_script("mrna_analysis/judge_differentregion.pl", ["--transwording","mRNA","--geneid","ID=gene:","--transid","ID=transcript:","--ymlinput",ranges], annotation)
        self.assertEqual(p.stdout.split("\t")[-2], "120-180")
        overlapping = gff("gene",1,100,"ID=gene:A;")+gff("gene",5,6,"ID=gene:B;","-")+gff("gene",10,20,"ID=gene:C;")
        p = self.run_script("mrna_analysis/filter_overlapgene.pl", data=overlapping)
        self.assertEqual(p.stdout, gff("gene",5,6,"ID=gene:B;","-"))

    def test_exon_and_codon_producers_accept_current_annotation(self):
        ranks = self.put("ranks.tsv", "TX1.1\tG1\t1\n")
        codon = self.run_script("mrna_analysis/codon_distance_1.pl", [ranks], "1\t180\t182\t+\tParent=transcript:TX1.1;\n")
        cf = self.put("codon.tsv", codon.stdout)
        annotation = "1\ttest\tmRNA\t100\t200\t.\t+\t.\tID=transcript:TX1.1;\n" + "".join(f"1\ttest\texon\t{a}\t{b}\t.\t+\t.\tParent=transcript:TX1.1;\n" for a,b in [(100,120),(180,200)])
        p = self.run_script("mrna_analysis/exon_distance_1.pl", ["--transid","Parent=transcript:",ranks], annotation)
        self.assertEqual(p.stdout.split("\t")[3], "100-120,180-200")
        model = self.run_script("mrna_analysis/codon_distance_2.pl", [cf], annotation)
        mf = self.put("model.tsv", model.stdout)
        sites = self.put("sites.tsv", "1\t119\t+\tA\tC\tG\tG1\n")
        p = self.run_script("mrna_analysis/codon_distance_3.pl", [mf,sites,self.d/"bar.tsv",0,0])
        self.assertEqual(p.stdout.split("\t")[5], "-2")

    def test_unordered_exons_and_missing_region_gene(self):
        lines = "1\tt\tmRNA\t1\t10\t.\t+\t.\tID=transcript:T1.1;\n1\tt\tmRNA\t20\t30\t.\t+\t.\tID=transcript:T1.2;\n1\tt\texon\t1\t10\t.\t+\t.\tParent=transcript:T1.1;\n"
        p = self.run_script("mrna_analysis/filter_nonsenseexon.pl", ["--transwording","mRNA","--transid","ID=transcript:","--exonid","Parent=transcript:"], lines)
        self.assertIn("Parent=transcript:T1.1;", p.stdout)
        model = self.put("r8.tsv", "1\tT1.1\tG1\t1\t+\t-\t1-10\t-\n")
        sites = self.put("sites.tsv", "1\t5\t+\tA\tC\tG\tMissing\n")
        p = self.run_script("mrna_analysis/main_transcript_3.pl", [model,sites,self.d/"annotation.tsv"])
        self.assertEqual(p.stdout, "")

    def test_paired_fastq_desynchronization_fails_and_empty_is_valid(self):
        r1, r2, out = self.d/"R1.fq.gz", self.d/"R2.fq.gz", self.d/"out.fq.gz"
        for path,name,seq in ((r1,"a","AC"),(r2,"b","GT")):
            with gzip.open(path,"wt") as f: f.write(f"@{name}\n{seq}\n+\nII\n")
        self.run_script("quality_control/pe_consistency.pl", [r1,r2,out], ok=False)
        for path in (r1,r2):
            with gzip.open(path,"wt") as f: pass
        p = self.run_script("quality_control/pe_consistency.pl", [r1,r2,out])
        self.assertIn("Proportion:\t0.00%",p.stdout)

    def test_rrna_score_zero_coverage_and_coordinate_mismatch(self):
        text = "".join(f"{p}\tA\t0\t0\n" for p in range(1,61))
        nc = self.put("NC/rrna.tsv",text)
        tr = self.put("TR/rrna.tsv",text)
        p = self.run_script("rrna_analysis/score.pl", [nc,tr])
        self.assertEqual(len(p.stdout.splitlines()),60) # header and L-1 evaluable positions
        self.assertNotIn("nan",p.stdout.lower())
        tr.write_text(text.replace("1\tA", "1\tC", 1))
        self.run_script("rrna_analysis/score.pl", [nc,tr], ok=False)

    def test_splice_union_coordinates_and_unknown_chromosome(self):
        ase = self.put("ase.tsv", "1\t1\t10\tx\tG1\t+\n1\t1\t10\tx\tG1\t+\n")
        cse = self.put("cse.tsv", "1\t11\t20\tx\tG1\t+\n")
        sites = self.put("sites.tsv", "1\t5\t+\n1\t15\t+\n2\t100\t+\n")
        self.assertEqual(self.run_script("presentation/judge_splice_exon.pl",[ase,cse,sites]).stdout,
                         "1\n1\n10\n10\n")

    def test_shell_failures_and_cleanup_do_not_touch_unrelated_logs(self):
        env = dict(os.environ)
        names = self.put("names.tsv", "a\n")
        p = subprocess.run(["bash",str(ROOT/"tool/extract_fastq.sh"),str(names),str(self.d/"missing1"),str(self.d/"o1"),str(self.d/"missing2"),str(self.d/"o2")],env=env,capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0)
        stub = self.put("bin/Rscript", "#!/bin/sh\nexit 0\n")
        stub.chmod(0o755)
        env["PATH"] = str(stub.parent)+os.pathsep+env["PATH"]
        table = self.put("points.tsv", "1\t100\t+\t20\n")
        collateral = self.put("plot.unrelated.log", "keep me\n")
        args = ["A",str(table),"4","B",str(table),"4","C",str(table),"4","plot","1"]
        p = subprocess.run(["bash",str(ROOT/"presentation/point_venn.sh"),*args],cwd=self.d,env=env,capture_output=True,text=True)
        self.assertEqual(p.returncode,0,p.stderr)
        self.assertTrue(collateral.exists())

    def test_delete_fastq_wrapper_uses_selected_perl(self):
        selected = shutil.which(PERL)
        self.assertIsNotNone(selected)
        stub = self.put('bin/perl', '#!/bin/sh\nexit 79\n')
        stub.chmod(0o755)
        env = dict(os.environ, NJU_PERL=selected,
                   PATH=str(stub.parent) + os.pathsep + os.environ['PATH'])
        names = self.put('names.tsv', 'remove/1\nremove/2\n')
        args = ['bash', str(ROOT / 'tool/extract_fastq.sh'), str(names)]
        for mate in (1, 2):
            source, out = self.d / f'R{mate}.fq.gz', self.d / f'kept{mate}.fq'
            with gzip.open(source, 'wt') as handle:
                handle.write(f'@remove/{mate}\nACGT\n+\nIIII\n'
                             f'@keep/{mate}\nTGCA\n+\nIIII\n')
            args.extend([str(source), str(out)])
        result = subprocess.run(args, env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        for mate in (1, 2):
            self.assertEqual((self.d / f'kept{mate}.fq').read_text(),
                             f'@keep/{mate}\nTGCA\n+\nIIII\n')

    def test_legacy_score_join_valid_schema_and_current_gene_annotation(self):
        base = ["1","100","+","A","C","G","G1","2","3"] + ["1"]*11 + ["x"]*6
        base[19] = "5"
        a = self.put("a.tsv", "\t".join(base)+"\n")
        base[19] = "6"
        b = self.put("b.tsv", "\t".join(base)+"\n")
        p = self.run_script("mrna_analysis/specific_point_withsocre.pl",[a,b])
        row = [r for r in p.stdout.splitlines() if not r.startswith("#")][0]
        self.assertEqual(row.split("\t")[-1], "11")
        sites = self.put("sites.tsv", "Chromosome\tPosition\tDirection\tLeftBase\tBase\tRightBase\tGeneID\n1\t100\t+\tA\tC\tG\tG1\n")
        gff = "1\tt\tgene\t1\t200\t.\t+\t.\tID=gene:G1;Name=Example;\n"
        p = self.run_script("tool/add_gene_name.pl", ["--id","ID=gene:","--name","Name=","--file",sites,"--col",8], gff)
        self.assertEqual(p.stdout.splitlines()[0].split("\t")[-1], "GeneName")
        self.assertEqual(p.stdout.splitlines()[1].split("\t")[-1], "Example")

    def test_representative_two_mane_tags_exclude_untagged_competitor(self):
        lines = "".join(f"1\tt\tmRNA\t1\t{end}\t.\t+\t.\tID=transcript:{tx};Parent=gene:G1;{tag}\n"
                        for tx,end,tag in [("T1.1",20,"tag=MANE_Select;"),("T1.2",30,"tag=MANE_Select;"),("T1.3",10,"tag=appris_principal_1;")])
        p = self.run_script("mrna_analysis/representative_transcript.pl", ["--geneid","Parent=gene:","--transid","ID=transcript:"], lines)
        self.assertEqual(p.stdout, "T1.1\tG1\n")

    def test_transcript_selectors_ignore_child_features(self):
        lines = "##gff-version 3\n" + "".join(
            f"1\tt\t{feature}\t1\t{end}\t.\t+\t.\tgene_id=G1;transcript_id={tx};tag=basic;\n"
            for feature,tx,end in [("mRNA","T1",20),("transcript","T2",30),
                                   ("CDS","T2",2),("exon","T3",1)])
        args = ["--geneid","gene_id=","--transid","transcript_id="]
        for script,expected in [("main_transcript_1.pl","T2\tG1\t1\nT1\tG1\t2\n"),
                                ("representative_transcript.pl","T1\tG1\n"),
                                ("filter_transcript.pl","T1\nT2\n")]:
            result = self.run_script("mrna_analysis/"+script, args, lines)
            self.assertEqual(result.stdout, expected)
            self.assertEqual(result.stderr, "")

    def test_relative_region_positions_mirror_on_both_strands(self):
        sites = "".join(f"1\t{pos}\t{strand}\tA\tC\tG\t{gene}\n"
                        for strand,gene in [("+","G1"),("-","G2")] for pos in [1,10])
        r7 = self.put("r7.tsv", "1\tG1\t+\tT1\t-\t1-10\t-\n1\tG2\t-\tT2\t-\t1-10\t-\n")
        r9 = self.put("r9.tsv", "1\tG1\t+\t-\t-\t-\t-\t1-10\t-\n1\tG2\t-\t-\t-\t-\t-\t1-10\t-\n")
        r8 = self.put("r8.tsv", "1\tT1\tG1\t1\t+\t-\t1-10\t-\n1\tT2\tG2\t1\t-\t-\t1-10\t-\n")
        points = self.put("points.tsv",sites)
        outputs = [self.run_script("mrna_analysis/judge_transcriptregion_4.pl",[r7],sites).stdout,
                   self.run_script("mrna_analysis/stat_differentregion_3.pl",[r9],sites).stdout,
                   self.run_script("mrna_analysis/main_transcript_3.pl",[r8,points,self.d/"annotated.tsv"]).stdout]
        expected = {("+",1):0.1,("+",10):1,("-",10):0.1,("-",1):1}
        for output in outputs:
            actual = {(row[1],int(row[2])):float(row[4]) for row in map(lambda s:s.split("\t"),output.splitlines())}
            self.assertEqual(actual,expected)

    def test_rrna_legacy_ttest_requires_explicit_interface_and_matching_rows(self):
        a = self.put("a.tsv", "1\tA\n")
        b = self.put("b.tsv", "1\tA\n")
        self.run_script("rrna_analysis/Nm_ttest.pl", [a,b], ok=False)
        self.run_script("rrna_analysis/Nm_ttest.pl", ["--legacy-three-replicates",a,b], ok=False)

    def test_empty_alignment_log_is_defined(self):
        log = self.put("sample.with.dots.bowtie2.log", "0 reads; of these:\n")
        p = self.run_script("tool/stat_alignment.pl", [log])
        self.assertIn("sample.with.dots\t0\t0\tNA\tNA", p.stdout)

    def test_almostunique_wrapper_keeps_its_distinct_legacy_policy(self):
        align = self.put("align.tsv", "a\tTX1\t1\t4=\tAAAA\nb\tTX1\t1\t4=\tAAAA\nb\tTX2\t1\t4=\tCCCC\n")
        fq = self.d/"reads.fq.gz"
        with gzip.open(fq,"wt") as f:f.write("@a\nTTTT\n+\nIIII\n@b\nAAAA\n+\nIIII\n")
        output = self.d/"legacy.tsv"
        sentinel = self.put("unique.tmp", "unrelated scratch\n")
        p = subprocess.run(["bash",str(ROOT/"mrna_analysis/almostunique.sh"),str(align),str(fq),str(self.d),str(output)],capture_output=True,text=True)
        self.assertEqual(p.returncode,0,p.stderr)
        self.assertEqual(sorted(r.split("\t")[0] for r in output.read_text().splitlines()),["a","b"])
        self.assertEqual(sentinel.read_text(),"unrelated scratch\n")
        direct = self.d/"direct.tsv"
        self.run_script("mrna_analysis/almostuniquematch.pl",[fq,align,direct])
        self.assertEqual([r.split("\t")[0] for r in direct.read_text().splitlines()],["b"])

    def test_density_consumers_and_gene_labeler(self):
        row="G1\t1\t+\t1-10\t10\t2\t11-20\t10\t0\t21-30\t10\t1\tF\n"
        for script in ["judge_transcriptregion_3.pl","stat_altersplice_2.pl"]:
            self.assertEqual(self.run_script("mrna_analysis/"+script,data=row).stdout,"G1\t0.2\t0\t0.1\n")
        region=self.put("regions.tsv","1\tG1\t+\tTX1\t1-10\t11-20\t21-30\n")
        splice=self.put("splice.tsv","1\tG1\t+\t1-30\t-\t-\n")
        point="1\t11\t+\tA\tC\tG\tG1\n"
        self.assertTrue(self.run_script("mrna_analysis/judge_region.pl",[splice,region],point).stdout.endswith("\tCSE\tCDS\n"))


if __name__ == "__main__":
    unittest.main()
