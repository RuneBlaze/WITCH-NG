use ahash::{AHashMap, RandomState};
use itertools::Itertools;
use ndarray::{Array, ShapeBuilder};
use std::collections::{BTreeSet};
use tracing::debug;

struct CompressedSimilarities {
    pub shape: (usize, usize),
    pub weights: AHashMap<(u32, u32), f64>,
    pub transl_row1: Vec<u32>,
    pub transl_row2: Vec<u32>,
}

fn compress_coordinates(weights: AHashMap<(u32, u32), f64>) -> CompressedSimilarities {
    let mut row1: BTreeSet<u32> = BTreeSet::new();
    let mut row2: BTreeSet<u32> = BTreeSet::new();
    for w in weights.keys() {
        row1.insert(w.0);
        row2.insert(w.1);
    }
    // compressed coordinate to global coordinate
    let compressed1 = row1.into_iter().collect_vec();
    let compressed2 = row2.into_iter().collect_vec();
    // global coordinate to compressed coordinate
    let mut to_compressed1: AHashMap<u32, u32> = AHashMap::new();
    let mut to_compressed2: AHashMap<u32, u32> = AHashMap::new();
    compressed1.iter().enumerate().for_each(|(i, &x)| {
        to_compressed1.insert(x, i as u32);
    });
    compressed2.iter().enumerate().for_each(|(i, &x)| {
        to_compressed2.insert(x, i as u32);
    });
    let new_weights: AHashMap<(u32, u32), f64, RandomState> =
        AHashMap::from_iter(weights.iter().map(|(&(x, y), &w)| {
            let (x, y) = (to_compressed1[&x], to_compressed2[&y]);
            ((x, y), w)
        }));
    let shape = (compressed1.len(), compressed2.len());
    CompressedSimilarities {
        shape,
        weights: new_weights,
        transl_row1: compressed1,
        transl_row2: compressed2,
    }
}

/// solving a name-less instance of maximum-weight alignment merging problem
pub fn solve_matching_problem(
    shape: (usize, usize),
    weights: AHashMap<(u32, u32), f64>,
) -> Vec<i32> {
    let compressed = compress_coordinates(weights);
    debug!(
        "compressed matching problem coordinates from {:?} to {:?}",
        shape, compressed.shape
    );
    let (n, m) = compressed.shape; // n : length of query, m : global column count
    let weights = compressed.weights;
    let mut s = Array::<f64, _>::zeros((n + 1, m + 1).f());
    let mut back = Array::<u8, _>::zeros((n + 1, m + 1).f());
    for i in 0..(n + 1) {
        for j in 0..(m + 1) {
            if i == 0 || j == 0 {
                s[[i, j]] = 0.0;
                continue;
            }
            let mut max = 0.0;
            let mut max_pt = 0u8;
            let w = weights
                .get(&((i - 1) as u32, (j - 1) as u32))
                .copied()
                .unwrap_or_default();
            let values = [s[[i - 1, j - 1]] + w, s[[i - 1, j]], s[[i, j - 1]]];
            for (i, &v) in values.iter().enumerate() {
                if i == 0 && w <= 0.0 {
                    max_pt = 1;
                    continue;
                }
                if v > max {
                    max = v;
                    max_pt = i as u8;
                }
            }
            s[[i, j]] = max;
            back[[i, j]] = max_pt;
        }
    }
    let mut res = vec![-1; shape.0];
    let (mut i, mut j) = (n, m);
    while i > 0 && j > 0 {
        let pt = back[[i, j]];
        if pt == 0 {
            i -= 1;
            j -= 1;
            res[compressed.transl_row1[i] as usize] = compressed.transl_row2[j] as i32;
        } else if pt == 1 {
            i -= 1;
        } else if pt == 2 {
            j -= 1;
        }
    }
    res
}
