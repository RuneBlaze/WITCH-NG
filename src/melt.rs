use crate::{external::hmmbuild, structures::*};
use ahash::AHashSet;
use fixedbitset::FixedBitSet;
use itertools::Itertools;
use ndarray::{Array, ShapeBuilder};
use ogcat::ogtree::*;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};
use seq_io::fasta::{Reader, Record};

use std::{
    collections::BinaryHeap,
    fs::{create_dir_all, File},
    io::BufWriter,
    path::PathBuf,
};
use tracing::info;

pub fn hierarchical_decomp(tree: &Tree, max_size: usize) -> TaxaHierarchy {
    let n = tree.ntaxa;
    let mut reordered_taxa = (0..n).collect::<Vec<_>>();
    let mut taxa_label = FixedBitSet::with_capacity(n); // Taxa ID -> is on the left
    let mut pq = BinaryHeap::new();
    let mut cuts = AHashSet::new();
    let mut decomposition_ranges: Vec<(usize, usize)> = Vec::new();
    cuts.insert(0usize);
    pq.push((tree.ntaxa, (0usize, tree.ntaxa), 0usize));
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
    decomposition_ranges.push((0usize, tree.ntaxa));
    while let Some((size, (lb, ub), root)) = pq.pop() {
        assert_eq!(size, ub - lb);
        if size <= max_size {
            break;
        }
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
        if tree_sizes[best_cut] >= 2 {
            decomposition_ranges.push((lb, lb + tree_sizes[best_cut] as usize));
        }
        if size - tree_sizes[best_cut] as usize > 2 {
            decomposition_ranges.push((lb + tree_sizes[best_cut] as usize, ub));
        }
        pq.push((
            tree_sizes[best_cut] as usize,
            (lb, lb + tree_sizes[best_cut] as usize),
            best_cut,
        ));
        pq.push((
            size - tree_sizes[best_cut] as usize,
            (lb + tree_sizes[best_cut] as usize, ub),
            root,
        ));
    }
    let mut taxa_positions: Vec<usize> = vec![0; n];
    for (p, t) in reordered_taxa.iter().enumerate() {
        taxa_positions[*t] = p;
    }

    TaxaHierarchy {
        reordered_taxa,
        taxa_positions,
        decomposition_ranges,
    }
}

pub fn oneshot_melt(
    input: &PathBuf,
    tree: &PathBuf,
    max_size: usize,
    outdir: &PathBuf,
) -> anyhow::Result<CrucibleCtxt> {
    let collection = TreeCollection::from_newick(tree).expect("Failed to read tree");
    let decomp = hierarchical_decomp(&collection.trees[0], max_size);
    info!(
        num_subsets = decomp.decomposition_ranges.len(),
        "decomposed input tree"
    );
    let mut reader = Reader::from_path(input)?;
    let mut records_failable: Result<Vec<_>, _> = reader.records().into_iter().collect();
    let records = records_failable.as_mut().unwrap();
    let ts = &collection.taxon_set;
    records.sort_unstable_by_key(|r| {
        let taxon_name = String::from_utf8(r.head.clone()).unwrap();
        let id = ts.to_id[&taxon_name];
        decomp.taxa_positions[id]
    });
    for (i, &t) in decomp.reordered_taxa.iter().enumerate() {
        assert_eq!(&String::from_utf8(records[i].head.clone())?, &ts.names[t]);
    }
    let n = records.len(); // # of seqs
    let k = records[0].seq.len(); // # of columns
    let mut nchars_prefix = Array::<u32, _>::zeros((n + 1, k).f());
    for i in 1..n + 1 {
        for j in 0..k {
            if i == 1 {
                nchars_prefix[[i, j]] = if records[i - 1].seq[j] == b'-' { 0 } else { 1 };
            } else {
                nchars_prefix[[i, j]] =
                    nchars_prefix[[i - 1, j]] + if records[i - 1].seq[j] == b'-' { 0 } else { 1 };
            }
        }
    }
    let subsets_root = outdir.join("subsets");
    let metadata_path = outdir.join("melt.json");
    create_dir_all(&subsets_root)?;
    for (i, &(lb, ub)) in decomp.decomposition_ranges.iter().enumerate() {
        let to_write = &records[lb..ub];
        let mut writer = BufWriter::new(File::create(subsets_root.join(format!("{}.afa", i)))?);
        for r in to_write {
            r.write_wrap(&mut writer, 60)?;
        }
    }

    decomp
        .decomposition_ranges
        .par_iter()
        .enumerate()
        .for_each(|(i, &(lb, ub))| {
            let to_write = &records[lb..ub];
            hmmbuild(
                to_write.iter(),
                format!("{}", i).as_str(),
                &subsets_root.join(format!("{}.hmm", i)),
            )
            .expect("Failed to build HMM");
        });

    let mut writer = BufWriter::new(File::create(metadata_path)?);
    let mut metadata: Vec<HmmMeta> = vec![];
    let mut buf = vec![0u32; k];
    for &decomp_range in &decomp.decomposition_ranges {
        CrucibleCtxt::retrieve_nchars_noalloc(&nchars_prefix, decomp_range, &mut buf);
        let mut nonzero_counts: Vec<u32> = vec![];
        let mut column_positions: Vec<usize> = vec![];
        for (i, &c) in buf.iter().enumerate() {
            if c > 0 {
                nonzero_counts.push(c);
                column_positions.push(i);
            }
        }
        let hmm = HmmMeta::new(decomp_range, nonzero_counts, column_positions);
        metadata.push(hmm);
    }
    let ctxt = CrucibleCtxt::new(metadata);
    serde_json::to_writer(&mut writer, &ctxt)?;
    Ok(ctxt)
}
