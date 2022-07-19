
use std::{io::Write, path::PathBuf};

use fixedbitset::FixedBitSet;
use itertools::Itertools;
use seq_io::{
    fasta::{OwnedRecord, Reader},
    BaseRecord,
};

/// data structure for keeping track of singleton columns efficiently
pub struct CompactHomologies {
    pub num_columns: usize,
    /// each sequence (by query seq id) and its positional homologies (-1 if not homologous to anything)
    pub homology_hits: Vec<Vec<i32>>,
}

pub struct FormattedHomologies {
    pub num_columns: usize,
    pub singleton_sets: Vec<FixedBitSet>,
    pub positions: Vec<Vec<u32>>,
}

impl CompactHomologies {
    pub fn new(num_columns: usize, homology_hits: Vec<Vec<i32>>) -> Self {
        Self {
            num_columns,
            homology_hits,
        }
    }

    pub fn append_consensus_column_hits(&mut self) {
        let new_positions = (0..(self.num_columns as i32)).collect_vec();
        self.homology_hits.push(new_positions);
    }

    pub fn transl(self) -> FormattedHomologies {
        let k = self.num_columns;
        let _n = self.homology_hits.len();
        let mut paddings = vec![0u32; k + 1];
        let seq_lengths = self.homology_hits.iter().map(|v| v.len()).collect_vec();
        let mut is_singletons = self
            .homology_hits
            .iter()
            .map(|h| FixedBitSet::with_capacity(h.len()))
            .collect_vec();
        for (i, hits) in self.homology_hits.iter().enumerate() {
            let mut padding_needed = 0u32;
            let mut last_hit = -1i32;
            let is_singleton = is_singletons.get_mut(i).unwrap();
            for (j, &h) in hits.iter().enumerate() {
                if h < 0 {
                    padding_needed += 1;
                    is_singleton.set(j as usize, true);
                } else {
                    paddings[h as usize] = padding_needed.max(paddings[h as usize]);
                    last_hit = h;
                    padding_needed = 0;
                }
            }
            if padding_needed > 0 {
                paddings[(last_hit + 1) as usize] =
                    padding_needed.max(paddings[(last_hit + 1) as usize]);
            }
        }
        let expanded_num_cols = paddings.iter().sum::<u32>() as usize + k;
        let mut shifted_columns = vec![0u32; k];
        for c in 0..k {
            if c == 0 {
                shifted_columns[c] = paddings[c];
            } else {
                shifted_columns[c] = shifted_columns[c - 1] + paddings[c] + 1;
            }
        }
        let positions: Vec<Vec<u32>> = self
            .homology_hits
            .into_iter()
            .enumerate()
            .map(|(i, mut hits)| {
                let seq_len = seq_lengths[i];
                let mut local_positions: Vec<u32> = vec![0u32; seq_len];
                let mut tracked_pos = -1i32;
                let mut num_backpositions = 0;
                let mut cursor = seq_len;
                while let Some(og_pos) = hits.pop() {
                    cursor -= 1;
                    if og_pos < 0 {
                        if tracked_pos < 0 {
                            num_backpositions += 1;
                            continue;
                        } else {
                            tracked_pos -= 1;
                            local_positions[cursor] = tracked_pos as u32;
                        }
                    } else {
                        let new_pos = shifted_columns[og_pos as usize];
                        tracked_pos = new_pos as i32;
                        local_positions[cursor] = new_pos;
                    }
                }
                if num_backpositions > 0 {
                    if tracked_pos < 0 {
                        for i in 0..num_backpositions {
                            local_positions[i] = i as u32;
                        }
                    } else {
                        let starting_pos = local_positions[seq_len - num_backpositions - 1];
                        for i in 0..num_backpositions {
                            local_positions[seq_len - num_backpositions + i] =
                                starting_pos + (i + 1) as u32;
                        }
                    }
                }
                local_positions
            })
            .collect();
        FormattedHomologies {
            num_columns: expanded_num_cols,
            singleton_sets: is_singletons,
            positions,
        }
    }
}

impl FormattedHomologies {
    pub fn write_all_sequences<W>(
        &self,
        queries: &[OwnedRecord],
        base_alignment_path: &PathBuf,
        w: &mut W,
    ) -> anyhow::Result<()>
    where
        W: Write,
    {
        let mut reader = Reader::from_path(base_alignment_path)?;
        let base_pos_map = self.positions.last().unwrap();
        while let Some(record_iffy) = reader.next() {
            let record = record_iffy?;
            w.write_all(b">")?;
            w.write_all(record.head())?;
            w.write_all(b"\n")?;
            let mut buf = vec![b'-'; self.num_columns];
            for (i, &c) in record.seq_lines().flatten().enumerate() {
                let target_pos = base_pos_map[i];
                buf[target_pos as usize] = c;
            }
            w.write_all(&buf)?;
            w.write_all(b"\n")?;
        }
        for (i, q) in queries.iter().enumerate() {
            w.write_all(b">")?;
            w.write_all(&q.head)?;
            w.write_all(b"\n")?;
            let mut buf: Vec<u8> = vec![b'-'; self.num_columns];
            for (j, &c) in q.seq.iter().enumerate() {
                let target_pos = self.positions[i][j];
                buf[target_pos as usize] = if self.singleton_sets[i][j] {
                    c.to_ascii_lowercase()
                } else {
                    c.to_ascii_uppercase()
                };
            }
            w.write_all(&buf)?;
            w.write_all(b"\n")?;
        }
        Ok(())
    }
}
