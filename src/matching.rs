use ahash::AHashMap;

/// represents a name-less instance of maximum-weight alignment merging problem
pub fn solve_matching_problem(shape : (u32, u32), weights: AHashMap<(u32, u32), f64>) -> Vec<i32> {
    let (n, m) = shape; // n : length of query, m : global column count
}