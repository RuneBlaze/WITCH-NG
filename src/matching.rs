use ahash::AHashMap;
use ndarray::{Array, ShapeBuilder};

/// represents a name-less instance of maximum-weight alignment merging problem
pub fn solve_matching_problem(
    shape: (usize, usize),
    weights: &AHashMap<(u32, u32), f64>,
) -> Vec<i32> {
    let (n, m) = shape; // n : length of query, m : global column count
    let mut s = Array::<f64, _>::zeros((n + 1, m + 1).f());
    let mut back = Array::<u8, _>::zeros((n + 1, m + 1).f());
    let mut res = vec![-1; n];
    for i in 0..(n + 1) {
        for j in 0..(m + 1) {
            if i == 0 || j == 0 {
                s[[i, j]] = 0.0;
                continue;
            }
            let mut max = 0.0;
            let mut max_pt = 0u8;
            let w = weights
                .get(&(i as u32, j as u32))
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
    let (mut i, mut j) = (n, m);
    while i > 0 && j > 0 {
        let pt = back[[i, j]];
        if pt == 0 {
            i -= 1;
            j -= 1;
            res[i] = j as i32;
        } else if pt == 1 {
            i -= 1;
        } else if pt == 2 {
            j -= 1;
        }
    }
    return res;
}
