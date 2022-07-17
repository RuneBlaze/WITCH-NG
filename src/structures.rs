use std::{fs::File, io::BufReader, path::Path};

use ndarray::{Array, Ix2};
use serde::{Deserialize, Serialize};

/// data structure for keeping track of singleton columns efficiently
pub struct CompactHomologies {}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct HmmMeta {
    pub sequence_range: (usize, usize),
    pub chars_cnt: Vec<u32>,
    pub column_poitions: Vec<usize>,
}

impl HmmMeta {
    pub fn new(
        sequence_range: (usize, usize),
        chars_cnt: Vec<u32>,
        column_poitions: Vec<usize>,
    ) -> Self {
        Self {
            sequence_range,
            chars_cnt,
            column_poitions,
        }
    }
}

pub struct TaxaHierarchy {
    pub reordered_taxa: Vec<usize>,
    pub decomposition_ranges: Vec<(usize, usize)>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CrucibleCtxt {
    pub version: u32,
    pub metadata: Vec<HmmMeta>,
}

impl CrucibleCtxt {
    pub fn new(metadata: Vec<HmmMeta>) -> Self {
        Self {
            version: 0,
            metadata,
        }
    }

    pub fn retrieve_nchars_noalloc(
        nchars_partial_sum: &Array<u32, Ix2>,
        sequence_range: (usize, usize),
        buf: &mut [u32],
    ) {
        let shape = nchars_partial_sum.shape();
        let (start, end) = sequence_range;
        let k = shape[1];
        for i in 0..k {
            buf[i] = nchars_partial_sum[(end, i)] - nchars_partial_sum[(start, i)];
        }
    }

    pub fn num_hmms(&self) -> usize {
        self.metadata.len()
    }
}

pub struct AdderPayload {
    /// a list of top hits tuple of HMM id and adjusted bitscore for each sequence
    pub sequence_tophits: Vec<Vec<(u32, f64)>>,
}

impl AdderPayload {
    pub fn from_path<P>(path: P) -> anyhow::Result<Self>
    where
        P: AsRef<Path>,
    {
        let tophits = serde_json::from_reader(BufReader::new(File::open(path)?))?;
        Ok(Self {
            sequence_tophits: tophits,
        })
    }

    /// consumes self and returns a mapping from HMM id to sequence and adjusted bitscores for hmmalign
    pub fn transpose(self, ctxt: &CrucibleCtxt) -> Vec<Vec<(u32, f64)>> {
        let n = ctxt.num_hmms();
        let mut res: Vec<Vec<(u32, f64)>> = vec![vec![]; n];
        for (seq_id, score_pairs) in self.sequence_tophits.into_iter().enumerate() {
            for (hmm_id, score) in score_pairs {
                res[hmm_id as usize].push((seq_id as u32, score));
            }
        }
        res
    }
}
