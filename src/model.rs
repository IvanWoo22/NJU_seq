//! Input records and transcript-to-genome coordinate models.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Alignment {
    pub read: String,
    pub reference: String,
    pub start: u64,
    pub cigar: String,
    pub sequence: String,
}

impl Alignment {
    pub fn parse(row: &str) -> Result<Self, String> {
        let fields: Vec<_> = row.split('\t').collect();
        if fields.len() != 5 {
            return Err("expected five alignment fields".into());
        }
        if fields[0].is_empty() || fields[1].is_empty() {
            return Err("empty read ID or reference".into());
        }
        if fields[2].is_empty() || !fields[2].bytes().all(|b| b.is_ascii_digit()) {
            return Err(format!("invalid alignment start {:?}", fields[2]));
        }
        let start = fields[2].parse::<u64>().map_err(|e| e.to_string())?;
        if start == 0 {
            return Err("alignment start must be positive".into());
        }
        Ok(Self {
            read: fields[0].into(),
            reference: fields[1].into(),
            start,
            cigar: fields[3].into(),
            sequence: fields[4].into(),
        })
    }

    pub fn to_tsv(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}",
            self.read, self.reference, self.start, self.cigar, self.sequence
        )
    }
}

pub fn parse_cigar(cigar: &str) -> Result<Vec<(u64, char)>, String> {
    let bytes = cigar.as_bytes();
    let mut offset = 0;
    let mut result = Vec::new();
    while offset < bytes.len() {
        let begin = offset;
        if !matches!(bytes[offset], b'1'..=b'9') {
            return Err(format!(
                "unsupported CIGAR {cigar:?}; require positive =,X,I,D operations"
            ));
        }
        while offset < bytes.len() && bytes[offset].is_ascii_digit() {
            offset += 1;
        }
        let n = cigar[begin..offset]
            .parse::<u64>()
            .map_err(|e| e.to_string())?;
        if offset == bytes.len() || !matches!(bytes[offset], b'=' | b'X' | b'I' | b'D') {
            return Err(format!(
                "unsupported CIGAR {cigar:?}; require end-to-end --xeq (=,X,I,D)"
            ));
        }
        result.push((n, bytes[offset] as char));
        offset += 1;
    }
    if result.is_empty() {
        return Err("empty CIGAR".into());
    }
    Ok(result)
}

#[derive(Clone, Debug)]
pub struct Transcript {
    pub chrom: String,
    pub strand: char,
    pub gene: String,
    pub exons: Vec<(u64, u64)>,
    pub offsets: Vec<u64>,
    pub length: u64,
}

impl Transcript {
    pub fn new(
        chrom: String,
        strand: char,
        gene: String,
        mut exons: Vec<(u64, u64)>,
    ) -> Result<Self, String> {
        if chrom.is_empty() || !matches!(strand, '+' | '-') || exons.is_empty() {
            return Err("transcript requires chromosome, strand and at least one exon".into());
        }
        exons.sort_unstable();
        let mut merged: Vec<(u64, u64)> = Vec::new();
        for (a, b) in exons {
            if a == 0 || b < a {
                return Err(format!("invalid exon {a}-{b}"));
            }
            if let Some(last) = merged.last_mut()
                && a <= last.1.saturating_add(1)
            {
                last.1 = last.1.max(b);
                continue;
            }
            merged.push((a, b));
        }
        let mut offsets = Vec::with_capacity(merged.len());
        let mut length = 0u64;
        for &(a, b) in &merged {
            offsets.push(length);
            length = length
                .checked_add(b - a + 1)
                .ok_or("transcript length overflow")?;
        }
        Ok(Self {
            chrom,
            strand,
            gene,
            exons: merged,
            offsets,
            length,
        })
    }

    /// Inclusive blocks, sorted in ascending genomic order on either strand.
    pub fn blocks(&self, start: u64, length: u64) -> Result<Vec<(u64, u64)>, String> {
        if start == 0 || length == 0 || start > self.length || length > self.length - start + 1 {
            return Err(format!(
                "range {start}+{length} exceeds transcript length {}",
                self.length
            ));
        }
        let (mut lo, mut hi) = (start - 1, start - 1 + length);
        if self.strand == '-' {
            (lo, hi) = (self.length - hi, self.length - lo);
        }
        let mut i = self.offsets.partition_point(|&offset| offset <= lo) - 1;
        let mut result = Vec::new();
        while lo < hi {
            let (a, b) = self.exons[i];
            let within_exon = lo - self.offsets[i];
            let take = (hi - lo).min(b - a + 1 - within_exon);
            let genomic_start = a + within_exon;
            result.push((genomic_start, genomic_start + (take - 1)));
            lo += take;
            i += 1;
        }
        Ok(result)
    }

    pub fn map_position(&self, position: u64) -> Result<u64, String> {
        Ok(self.blocks(position, 1)?[0].0)
    }
}

#[derive(Clone, Debug)]
pub struct ExonIndex {
    pub transcripts: HashMap<String, Transcript>,
}

impl ExonIndex {
    pub fn load(path: &Path) -> Result<Self, String> {
        Self::load_with(path, "auto", "")
    }

    pub fn load_with(path: &Path, refstr: &str, transid: &str) -> Result<Self, String> {
        let reader =
            BufReader::new(File::open(path).map_err(|e| format!("{}: {e}", path.display()))?);
        Self::from_reader(reader, refstr, transid).map_err(|e| format!("{}: {e}", path.display()))
    }

    pub fn from_reader<R: BufRead>(
        mut reader: R,
        refstr: &str,
        transid: &str,
    ) -> Result<Self, String> {
        if refstr.is_empty() {
            return Err("empty annotation reference prefix".into());
        }
        type RawTranscript = (String, char, String, Vec<(u64, u64)>);
        let mut raw: HashMap<String, RawTranscript> = HashMap::new();
        let mut line = String::new();
        let mut number = 0;
        loop {
            line.clear();
            if reader.read_line(&mut line).map_err(|e| e.to_string())? == 0 {
                break;
            }
            number += 1;
            let fields: Vec<_> = line.trim_end_matches(['\r', '\n']).split('\t').collect();
            if fields.len() != 5 {
                return Err(format!("exon.info line {number}: expected five fields"));
            }
            let annotation = if refstr == "auto" {
                supported_annotation(fields[4])?
            } else {
                annotation_transcript(fields[4], refstr, transid).map(|tx| {
                    let gene = tx
                        .rsplit_once('.')
                        .map_or(tx.as_str(), |(g, _)| g)
                        .to_owned();
                    (tx, gene)
                })
            };
            let Some((tx, gene)) = annotation else {
                continue;
            };
            let a = fields[1]
                .parse::<u64>()
                .map_err(|e| format!("exon.info line {number}: {e}"))?;
            let b = fields[2]
                .parse::<u64>()
                .map_err(|e| format!("exon.info line {number}: {e}"))?;
            let strand = match fields[3] {
                "+" => '+',
                "-" => '-',
                _ => return Err(format!("exon.info line {number}: invalid strand")),
            };
            if a == 0 || b < a {
                return Err(format!("exon.info line {number}: invalid exon"));
            }
            // The maintained Perl mainline strips a chromosome prefix only.
            let chrom = fields[0]
                .strip_prefix("chr")
                .unwrap_or(fields[0])
                .to_string();
            let item = raw
                .entry(tx.clone())
                .or_insert_with(|| (chrom.clone(), strand, gene.clone(), Vec::new()));
            if item.0 != chrom || item.1 != strand || item.2 != gene {
                return Err(format!("conflicting annotation for {tx}"));
            }
            item.3.push((a, b));
        }
        let mut transcripts = HashMap::new();
        for (id, (chrom, strand, gene, exons)) in raw {
            transcripts.insert(id, Transcript::new(chrom, strand, gene, exons)?);
        }
        if transcripts.is_empty() {
            return Err("no transcripts matched the annotation options".into());
        }
        Ok(Self { transcripts })
    }
}

/// Supported local five-column annotations: Ensembl Arabidopsis and GENCODE
/// human. Attribute boundaries and versioned IDs are literal, not regex guesses.
fn supported_annotation(attributes: &str) -> Result<Option<(String, String)>, String> {
    let get = |key: &str| {
        attributes.split(';').find_map(|a| {
            let (k, v) = a.split_once('=')?;
            (k == key).then_some(v)
        })
    };
    let Some(parent) = get("Parent") else {
        return Ok(None);
    };
    if let Some(tx) = parent.strip_prefix("transcript:")
        && tx.starts_with("AT")
    {
        let tx = annotation_transcript(attributes, "Parent=transcript:", "AT")
            .ok_or_else(|| format!("Invalid Arabidopsis transcript: {parent}"))?;
        let gene = tx.rsplit_once('.').unwrap().0.to_owned();
        return Ok(Some((tx, gene)));
    }
    let tx = parent.strip_prefix("transcript:").unwrap_or(parent);
    if !tx.starts_with("ENST") {
        return Ok(None);
    }
    human_stable_id(tx, "ENST")?;
    let gene = get("gene_id").ok_or_else(|| format!("Missing gene_id for {tx}"))?;
    let stable_gene = human_stable_id(gene, "ENSG")?;
    // GENCODE exon gene_id may be unsuffixed even for a _PAR_Y Parent.
    // Preserve the annotation's gene identity, and the Parent's full identity.
    if gene.ends_with("_PAR_Y") && !tx.ends_with("_PAR_Y") {
        return Err(format!("Inconsistent PAR_Y transcript/gene: {tx}, {gene}"));
    }
    Ok(Some((tx.to_owned(), stable_gene)))
}

/// Remove only the numeric version, keeping the GENCODE PAR_Y identity.
pub fn human_stable_id(id: &str, prefix: &str) -> Result<String, String> {
    let (base, suffix) = id
        .strip_suffix("_PAR_Y")
        .map_or((id, ""), |b| (b, "_PAR_Y"));
    let (stem, version) = base
        .split_once('.')
        .ok_or_else(|| format!("Invalid human ID: {id}"))?;
    let number = stem.strip_prefix(prefix).unwrap_or("");
    if number.is_empty()
        || !number.bytes().all(|b| b.is_ascii_digit())
        || version.is_empty()
        || !version.bytes().all(|b| b.is_ascii_digit())
    {
        return Err(format!("Invalid human ID: {id}"));
    }
    Ok(format!("{stem}{suffix}"))
}

fn annotation_transcript(attributes: &str, refstr: &str, transid: &str) -> Option<String> {
    for (offset, _) in attributes.match_indices(refstr) {
        let tail = &attributes[offset + refstr.len()..];
        let Some(rest) = tail.strip_prefix(transid) else {
            continue;
        };
        let stem_len = rest
            .char_indices()
            .take_while(|(_, c)| c.is_alphanumeric() || *c == '_')
            .map(|(i, c)| i + c.len_utf8())
            .last()
            .unwrap_or(0);
        if stem_len == 0 {
            continue;
        }
        let Some(version) = rest[stem_len..].strip_prefix('.') else {
            continue;
        };
        let version_len = version.bytes().take_while(|b| b.is_ascii_digit()).count();
        if version_len == 0 {
            continue;
        }
        return Some(format!(
            "{transid}{}.{}",
            &rest[..stem_len],
            &version[..version_len]
        ));
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn human_gene_ids_strands_splice_boundaries_and_par_copies() {
        let input = concat!(
            "chrX\t101\t105\t+\tParent=ENST00001.2;gene_id=ENSG00009.3\n",
            "chrX\t201\t205\t+\tParent=ENST00001.2;gene_id=ENSG00009.3\n",
            "chrY\t301\t305\t+\tParent=ENST00001.2_PAR_Y;gene_id=ENSG00009.3\n",
            "chrY\t401\t405\t+\tParent=ENST00001.2_PAR_Y;gene_id=ENSG00009.3\n",
            "chr2\t501\t505\t-\tParent=ENST00002.4;gene_id=ENSG00008.1\n",
            "chr2\t601\t605\t-\tParent=ENST00002.4;gene_id=ENSG00008.1\n",
        );
        let index = ExonIndex::from_reader(input.as_bytes(), "auto", "").unwrap();
        assert_eq!(index.transcripts.len(), 3);
        let x = &index.transcripts["ENST00001.2"];
        let y = &index.transcripts["ENST00001.2_PAR_Y"];
        assert_eq!(x.gene, "ENSG00009");
        assert_eq!(y.gene, "ENSG00009");
        assert_eq!(y.chrom, "Y");
        assert_eq!(x.blocks(4, 4).unwrap(), vec![(104, 105), (201, 202)]);
        assert_eq!(y.blocks(4, 4).unwrap(), vec![(304, 305), (401, 402)]);
        let reverse = &index.transcripts["ENST00002.4"];
        assert_eq!(reverse.map_position(1).unwrap(), 605);
        assert_eq!(reverse.map_position(6).unwrap(), 505);
        assert_eq!(reverse.blocks(4, 4).unwrap(), vec![(504, 505), (601, 602)]);
        assert_eq!(
            human_stable_id("ENSG00009.3_PAR_Y", "ENSG").unwrap(),
            "ENSG00009_PAR_Y"
        );
    }

    #[test]
    fn human_invalid_or_conflicting_annotation_fails() {
        for input in [
            "chr1\t1\t3\t+\tParent=ENST00001.1\n",
            "chr1\t1\t3\t+\tParent=ENST00001.1;gene_id=ENST00001.1\n",
            "chr1\t1\t3\t+\tParent=ENST00001.1_bad;gene_id=ENSG00001.1\n",
            "chr1\t1\t3\t+\tParent=ENST00001.1;gene_id=ENSG00001.1\nchr1\t5\t7\t+\tParent=ENST00001.1;gene_id=ENSG00002.1\n",
            "chr1\t1\t3\t+\tParent=ENST00001.1;gene_id=ENSG00001.1\nchrY\t1\t3\t+\tParent=ENST00001.1;gene_id=ENSG00001.1\n",
            "chr1\t1\t3\t+\tOtherParent=ENST00001.1;gene_id=ENSG00001.1\n",
        ] {
            assert!(
                ExonIndex::from_reader(input.as_bytes(), "auto", "").is_err(),
                "{input}"
            );
        }
    }

    #[test]
    fn positions_and_ranges_match_explicit_bases() {
        for strand in ['+', '-'] {
            let tx = Transcript::new(
                "1".into(),
                strand,
                "AT1".into(),
                vec![(101, 110), (201, 220)],
            )
            .unwrap();
            let mut bases: Vec<u64> = (101..=110).chain(201..=220).collect();
            if strand == '-' {
                bases.reverse();
            }
            for start in 1..=tx.length {
                assert_eq!(tx.map_position(start).unwrap(), bases[start as usize - 1]);
                for length in 1..=tx.length - start + 1 {
                    let observed: Vec<_> = tx
                        .blocks(start, length)
                        .unwrap()
                        .iter()
                        .flat_map(|&(a, b)| a..=b)
                        .collect();
                    let mut expected =
                        bases[(start - 1) as usize..(start - 1 + length) as usize].to_vec();
                    expected.sort_unstable();
                    assert_eq!(observed, expected);
                }
            }
            assert!(tx.blocks(0, 1).is_err());
            assert!(tx.blocks(30, 2).is_err());
        }
    }

    #[test]
    fn exon_union_and_unmatched_annotation() {
        let input = "chr1\t101\t110\t+\tParent=transcript:AT1G00010.1;\nchr1\t105\t115\t+\tParent=transcript:AT1G00010.1;\nchr1\t116\t120\t+\tParent=transcript:AT1G00010.1;\nchr1\t999\t999\t+\tParent=transcript:at_other;\n";
        let index = ExonIndex::from_reader(input.as_bytes(), "Parent=transcript:", "AT").unwrap();
        let model = &index.transcripts["AT1G00010.1"];
        assert_eq!(model.chrom, "1");
        assert_eq!(model.gene, "AT1G00010");
        assert_eq!(model.exons, vec![(101, 120)]);
        assert_eq!(model.length, 20);
    }

    #[test]
    fn strict_cigar_and_alignment_parsing() {
        assert_eq!(
            parse_cigar("3=2I4=1X2D").unwrap(),
            vec![(3, '='), (2, 'I'), (4, '='), (1, 'X'), (2, 'D')]
        );
        for text in ["", "8M", "0=", "08=", "8=1", "1S8=", "8=oops"] {
            assert!(parse_cigar(text).is_err(), "{text}");
        }
        let row = "read_AAAAAAAA\tAT1.1\t8\t8=\tACGTACGT";
        assert_eq!(Alignment::parse(row).unwrap().to_tsv(), row);
        assert!(Alignment::parse("read\tAT1.1\t0\t8=\tACGTACGT").is_err());
        assert!(Alignment::parse("not five fields").is_err());
    }
}
