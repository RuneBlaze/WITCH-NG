use std::{
    cmp::Reverse,
    fs::File,
    io::{BufReader, BufWriter},
    path::PathBuf,
};

use ahash::AHashMap;
use itertools::Itertools;
use ordered_float::NotNan;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use seq_io::fasta::OwnedRecord;
use tracing::{debug, info};

use crate::{
    external::hmmsearch,
    structures::{AdderPayload, CrucibleCtxt},
};

pub struct ScoringCtxt {
    pub base_dir: PathBuf,
    pub hmm_ctxt: CrucibleCtxt,
    pub queries: Vec<OwnedRecord>,
    pub seq_ids: AHashMap<String, u32>,
}

#[derive(Debug, Clone)]
pub struct BitscoreTracker {
    pub hmm_ids: Vec<u32>,
    pub bitscores: Vec<f64>,
}

impl Default for BitscoreTracker {
    fn default() -> Self {
        Self {
            hmm_ids: Vec::new(),
            bitscores: Vec::new(),
        }
    }
}

impl BitscoreTracker {
    pub fn calc_adjusted_scores(&self, ctxt: &ScoringCtxt) -> impl Iterator<Item = (u32, f64)> {
        let hmm_sizes = self
            .hmm_ids
            .iter()
            .map(|i| ctxt.hmm_ctxt.metadata[*i as usize].num_seqs())
            .collect_vec();
        let mut converted = self
            .hmm_ids
            .iter()
            .zip(self.bitscores.iter())
            .map(|(hmm_id, score_i)| {
                // let hmm_size
                let size_i = ctxt.hmm_ctxt.metadata[*hmm_id as usize].num_seqs();
                let exponents = self
                    .bitscores
                    .iter()
                    .zip(hmm_sizes.iter())
                    .map(|(b, s)| b - score_i + (*s as f64 / size_i as f64).log2());
                let denominator = exponents.map(|e| 2.0f64.powf(e)).sum::<f64>();
                (Reverse(NotNan::new(1.0 / denominator).unwrap()), *hmm_id)
            })
            .collect_vec();
        if converted.len() > 10 {
            converted.select_nth_unstable(9);
        }
        converted.truncate(10);
        converted.into_iter().map(|(s, c)| (c, s.0.into_inner()))
    }
}

impl ScoringCtxt {
    pub fn manual_construction(base_dir: &PathBuf) -> anyhow::Result<Self> {
        let hmm_ctxt_path = base_dir.join("melt.json");
        let queries_path = base_dir.parent().unwrap().join("queries.fasta");
        let hmm_ctxt = serde_json::from_reader(BufReader::new(File::open(&hmm_ctxt_path)?))?;
        let queries_failiable: Result<Vec<_>, _> =
            seq_io::fasta::Reader::new(File::open(&queries_path)?)
                .records()
                .into_iter()
                .collect();
        let queries = queries_failiable?;
        let mut seq_ids: AHashMap<String, u32> = AHashMap::new();
        for (i, q) in queries.iter().enumerate() {
            seq_ids.insert(String::from_utf8(q.head.clone())?, i as u32);
        }
        Ok(Self {
            base_dir: base_dir.to_owned(),
            hmm_ctxt,
            queries,
            seq_ids,
        })
    }

    pub fn from_ehmms_ctxt(
        base_dir: PathBuf,
        hmm_ctxt: CrucibleCtxt,
        queries_path: &PathBuf,
    ) -> anyhow::Result<Self> {
        let queries_failiable: Result<Vec<_>, _> =
            seq_io::fasta::Reader::new(File::open(&queries_path)?)
                .records()
                .into_iter()
                .collect();
        let queries = queries_failiable?;
        info!("read {} query sequences", queries.len());
        let mut seq_ids: AHashMap<String, u32> = AHashMap::new();
        for (i, q) in queries.iter().enumerate() {
            seq_ids.insert(String::from_utf8(q.head.clone())?, i as u32);
        }
        Ok(Self {
            base_dir,
            hmm_ctxt,
            queries,
            seq_ids,
        })
    }

    pub fn hmm_path(&self, hmm_id: u32) -> PathBuf {
        self.base_dir
            .join("subsets")
            .join(format!("{}.hmm", hmm_id))
    }
    pub fn scores_path(&self) -> PathBuf {
        self.base_dir.join("scores.json")
    }

    pub fn produce_payload(&self) -> anyhow::Result<AdderPayload> {
        let h = self.hmm_ctxt.num_hmms();
        let q = self.queries.len();
        let mut score_trackers = vec![BitscoreTracker::default(); q];
        let hmmsearch_results: Vec<(u32, u32, f64)> = (0..h)
            .into_par_iter()
            .flat_map_iter(|i| {
                debug!("scoring hmm {}", i);
                let hmm_path = self.hmm_path(i as u32);
                let search_res = hmmsearch(&hmm_path, self.queries.iter(), &self.seq_ids)
                    .expect("hmmsearch failed");
                search_res.into_iter().map(move |(b, c)| (i as u32, b, c))
            })
            .collect();
        for (hmm_id, seq_id, score) in hmmsearch_results {
            score_trackers[seq_id as usize].hmm_ids.push(hmm_id);
            score_trackers[seq_id as usize].bitscores.push(score);
        }
        let new_scores: Vec<Vec<(u32, f64)>> = score_trackers
            .par_iter()
            .map(|st| st.calc_adjusted_scores(self).collect_vec())
            .collect();
        Ok(AdderPayload {
            sequence_tophits: new_scores,
        })
    }
}

pub fn oneshot_score_queries(basedir: &PathBuf) -> anyhow::Result<()> {
    let ctxt = ScoringCtxt::manual_construction(basedir)?;
    let payload = ctxt.produce_payload()?;
    let mut w = BufWriter::new(File::create(ctxt.scores_path())?);
    serde_json::to_writer(&mut w, &payload.sequence_tophits)?;
    Ok(())
}
