use ahash::AHashSet;
use fixedbitset::FixedBitSet;
use itertools::Itertools;
use ogcat::ogtree::*;
use seq_io::fasta::Reader;
use std::{path::PathBuf, collections::BinaryHeap};

pub struct TaxaHierarchy {
    pub reordered_taxa: Vec<usize>,
    pub decomposition_ranges : Vec<(usize, usize)>,
}

pub fn hierarchical_decomp(tree: &Tree, max_size: usize) -> TaxaHierarchy {
    let n = tree.ntaxa;
    let mut reordered_taxa = (0..n).collect::<Vec<_>>();
    let mut taxa_label = FixedBitSet::with_capacity(n); // Taxa ID -> is on the left
    let mut pq = BinaryHeap::new();
    let mut cuts = AHashSet::new();
    let mut decomposition_ranges : Vec<(usize, usize)> = Vec::new();
    cuts.insert(0usize);
    pq.push(((0usize, tree.ntaxa), 0usize));
    let mut tree_sizes = vec![0u64; tree.taxa.len()];
    for i in tree.postorder() {
        if tree.is_leaf(i) {
            tree_sizes[i] = 1;
        } else {
            tree.children(i).for_each(|c| {
                tree_sizes[i] += tree_sizes[c];
            });
        }
    }
    while let Some(((lb, ub), root)) = pq.pop() {
        let size = ub - lb;
        if size < max_size {
            break;
        }
        decomposition_ranges.push((lb, ub));
        let it = PostorderIterator::from_node_excluding(tree, root, &cuts);
        let mut best_inbalance = u64::MAX;
        let mut best_cut = 0usize;
        let mut non_leaf = false;
        for i in it {
            if i == root {
                continue;
            }
            if tree.is_leaf(i) {
            } else {
                non_leaf = true;
                let inbalance = (size as u64 - tree_sizes[i]).abs_diff(tree_sizes[i]);
                if inbalance < best_inbalance {
                    best_inbalance = inbalance;
                    best_cut = i;
                }
            }
        } // finding the best cut
        if non_leaf {
            assert_ne!(best_inbalance, u64::MAX, "No cut found");
        } else {
            break;
        }
        for a in tree.ancestors(best_cut) {
            if a == root {
                break;
            }
            tree_sizes[a] -= tree_sizes[best_cut];
        }
        cuts.insert(best_cut);
        for u in tree.postorder_from(best_cut) {
            if tree.is_leaf(u) {
                let tid = tree.taxa[u] as usize;
                taxa_label.set(tid, true);
            }
        }
        let view = &mut reordered_taxa[lb..ub];
        view.sort_unstable_by_key(|e| !taxa_label[*e]);
        taxa_label.clear();
        pq.push(((lb, lb + tree_sizes[best_cut] as usize), best_cut));
        pq.push(((lb + tree_sizes[best_cut] as usize, ub), root));
    }
    TaxaHierarchy {
        reordered_taxa,
        decomposition_ranges,
    }
}

pub fn oneshot_melt(
    input: &PathBuf,
    tree: &PathBuf,
    max_size: usize,
    outdir: &PathBuf,
) -> anyhow::Result<()> {
    let collection = TreeCollection::from_newick(tree).expect("Failed to read tree");
    let decomp = hierarchical_decomp(&collection.trees[0], max_size);
    let mut reader = Reader::from_path(input)?;
    let mut records_failable : Result<Vec<_>, _> = reader.records().into_iter().into_iter().collect();
    let records = records_failable.as_mut().unwrap();
    let ts = &collection.taxon_set;
    records.sort_unstable_by_key(|r| {
        let taxon_name = String::from_utf8(r.head.clone()).unwrap();
        let id = ts.to_id[&taxon_name];
        decomp.reordered_taxa[id]
    });
    for (i, &(lb, ub)) in decomp.decomposition_ranges.iter().enumerate() {
        let to_write = &records[lb..ub];
    }
    Ok(())
}