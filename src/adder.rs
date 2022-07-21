use crate::{
    compact_printer::CompactHomologies,
    external,
    matching::solve_matching_problem,
    score_calc::ScoringCtxt,
    structures::{AdderPayload, CrucibleCtxt},
};
use ahash::AHashMap;
use anyhow::bail;
use itertools::Itertools;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use seq_io::{fasta::OwnedRecord, BaseRecord};
use std::{
    cell::RefCell,
    fs::File,
    io::{BufReader, BufWriter},
    path::PathBuf,
    sync::Arc,
};
use thread_local::ThreadLocal;
use tracing::info;

pub struct AdderContext {
    base_dir: PathBuf,
    hmm_ctxt: CrucibleCtxt,
    queries: Vec<OwnedRecord>,
    transposed_scores: Vec<Vec<(u32, f64)>>,
}

impl AdderContext {
    pub fn from_scoring_ctxt(
        base_dir: &PathBuf,
        scorer: ScoringCtxt,
        payload: AdderPayload,
    ) -> anyhow::Result<Self> {
        let transposed = payload.transpose(&scorer.hmm_ctxt);
        let (queries, hmm_ctxt) = (scorer.queries, scorer.hmm_ctxt);
        Ok(Self {
            base_dir: base_dir.to_owned(),
            hmm_ctxt,
            queries,
            transposed_scores: transposed,
        })
    }
    pub fn manual_construction(base_dir: &PathBuf) -> anyhow::Result<Self> {
        let hmm_ctxt_path = base_dir.join("melt.json");
        let scores_path = base_dir.join("scores.json");
        let queries_path = base_dir.parent().unwrap().join("queries.fasta");
        let hmm_ctxt = serde_json::from_reader(BufReader::new(File::open(&hmm_ctxt_path)?))?;
        let transposed = AdderPayload::from_path(&scores_path)?.transpose(&hmm_ctxt);
        let queries_failiable: Result<Vec<_>, _> =
            seq_io::fasta::Reader::new(File::open(&queries_path)?)
                .records()
                .into_iter()
                .collect();
        let queries = queries_failiable?;
        Ok(Self {
            base_dir: base_dir.to_owned(),
            hmm_ctxt,
            queries,
            transposed_scores: transposed,
        })
    }

    pub fn base_alignment_path(&self) -> PathBuf {
        self.base_dir.parent().unwrap().join("backbone.aln.fasta")
    }

    pub fn default_output_path(&self) -> PathBuf {
        self.base_dir.parent().unwrap().join("merged.afa")
    }
}

#[derive(Debug)]
pub struct Subweights {
    pub weights: Vec<AHashMap<(u32, u32), f64>>,
}

impl Subweights {
    pub fn from_ctxt(ctxt: &AdderContext) -> Self {
        let n = ctxt.queries.len();
        Subweights {
            weights: vec![AHashMap::new(); n],
        }
    }

    pub fn merge_in(&mut self, rhs: Subweights) {
        for (w, r) in self.weights.iter_mut().zip(rhs.weights) {
            for (k, v) in r {
                w.entry(k).and_modify(|e| *e += v).or_insert(v);
            }
        }
    }
}

impl AdderContext {
    pub fn hmm_path(&self, hmm_id: u32) -> PathBuf {
        self.base_dir
            .join("subsets")
            .join(format!("{}.hmm", hmm_id))
    }

    pub fn process_one_hmm(
        &self,
        hmm_id: u32,
        subweights: &mut Subweights,
        // hits: &[(u32, f64)],
    ) -> anyhow::Result<()> {
        let metadata = &self.hmm_ctxt.metadata[hmm_id as usize];
        let hits = &self.transposed_scores[hmm_id as usize];
        if hits.is_empty() {
            return Ok(());
        }
        let queries_for_hmm = hits
            .iter()
            .map(|&(seq_id, _)| &self.queries[seq_id as usize]);
        let hmm_path = self.hmm_path(hmm_id);
        let raw_afa = queries_for_hmm
            .chunks(500)
            .into_iter()
            .flat_map(|c| external::hmmalign(&hmm_path, c).expect("hmmalign failed"))
            .collect_vec();
        let raw_afa_view: &[u8] = &raw_afa;
        let mut reader = seq_io::fasta::Reader::new(raw_afa_view);
        let mut record_id = 0usize;
        while let Some(unverified_record) = reader.next() {
            let record = unverified_record?;
            let seq_id = hits[record_id].0;
            let seq_weight = hits[record_id].1;
            let mut residue_ix = 0u32; // which character of the query are we at?
            let mut column_ix = 0u32; // which column of the consensus are we at?
            for &c in record.seq_lines().flatten() {
                match c {
                    b'.' => {
                        continue;
                    }
                    b'-' => {
                        column_ix += 1;
                    }
                    _ if c.is_ascii_uppercase() => {
                        let weight_delta =
                            seq_weight * metadata.chars_cnt[column_ix as usize] as f64;
                        let global_column = metadata.column_poitions[column_ix as usize];
                        subweights.weights[seq_id as usize]
                            .entry((residue_ix, global_column as u32))
                            .and_modify(|w| *w += weight_delta)
                            .or_insert(weight_delta);
                        residue_ix += 1;
                        column_ix += 1;
                    }
                    _ if c.is_ascii_lowercase() => {
                        residue_ix += 1;
                    }
                    _ => {
                        bail!("Unexpected character in alignment: {}", c as char);
                    }
                }
            }
            assert_eq!(column_ix, metadata.column_poitions.len() as u32);
            record_id += 1;
        }
        Ok(())
    }
}

pub fn unoptimized_process_transposed_payload(ctxt: &AdderContext) -> anyhow::Result<Subweights> {
    let tls = Arc::new(ThreadLocal::new());
    (0..ctxt.hmm_ctxt.num_hmms())
        .into_par_iter()
        .for_each(|hmm_id| {
            let local = tls.clone();
            let subweights = local.get_or(|| RefCell::new(Subweights::from_ctxt(ctxt)));
            let mut borrowed = subweights.borrow_mut();
            ctxt.process_one_hmm(hmm_id as u32, &mut borrowed)
                .expect("Failed to run hmmalign.");
        });
    let mut subweights = Subweights::from_ctxt(ctxt);
    Arc::try_unwrap(tls).unwrap().into_iter().for_each(|s| {
        subweights.merge_in(s.into_inner());
    });
    Ok(subweights)
}

pub fn oneshot_add_queries(basedir: &PathBuf) -> anyhow::Result<()> {
    let ctxt = AdderContext::manual_construction(basedir)?;
    let default_output_path = ctxt.default_output_path();
    let base_alignment_path = ctxt.base_alignment_path();
    add_queries(ctxt, &default_output_path, &base_alignment_path)
}

pub fn add_queries(
    ctxt: AdderContext,
    outfile: &PathBuf,
    base_alignment_path: &PathBuf,
) -> anyhow::Result<()> {
    let subweights = unoptimized_process_transposed_payload(&ctxt)?;
    let m = ctxt.hmm_ctxt.metadata[0].column_poitions.len();
    let dp_solutions: Vec<Vec<i32>> = subweights
        .weights
        .into_par_iter()
        .enumerate()
        .map(|(i, w)| {
            let n = ctxt.queries[i].seq.len();
            solve_matching_problem((n, m), w)
        })
        .collect();
    let mut c_homologies =
        CompactHomologies::new(ctxt.hmm_ctxt.num_consensus_columns(), dp_solutions);
    c_homologies.append_consensus_column_hits();
    let hydrated_homologies = c_homologies.transl();
    info!(
        "output homologies formatted, output alignment will have {} columns",
        hydrated_homologies.num_columns
    );
    let mut output_writer = BufWriter::new(File::create(outfile)?);
    hydrated_homologies.write_all_sequences(
        &ctxt.queries,
        base_alignment_path,
        &mut output_writer,
    )?;
    Ok(())
}
