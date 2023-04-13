use std::{io::Write, path::PathBuf};

use fixedbitset::FixedBitSet;
use itertools::Itertools;
use seq_io::{
    fasta::{OwnedRecord, Reader},
    BaseRecord,
};

/// data structure for keeping track of singleton columns efficiently
pub struct LettersWithColors {
    /// number of homology equivalence classes in the backbone
    pub num_colors: usize,
    /// each sequence (by query seq id) and its positional homologies (-1 if not homologous to anything)
    pub letter_colors: Vec<Vec<i32>>,
}

/// the target data structure derived from CompactHomologies, see CompactHomologies::transl()
pub struct FormattedHomologies {
    /// number of actual columns in the output printtable MSA
    pub num_visual_columns: usize,
    /// which letters (ix by query seq id) are lower case letters (singleton columns)
    pub singleton_letters: Vec<FixedBitSet>,
    /// where each letter should go (ix by query seq id)
    pub letter_positions: Vec<Vec<u32>>,
}

impl LettersWithColors {
    pub fn new(num_columns: usize, homology_hits: Vec<Vec<i32>>) -> Self {
        Self {
            num_colors: num_columns,
            letter_colors: homology_hits,
        }
    }

    pub fn append_backbone_column_colors(&mut self) {
        let new_positions = (0..(self.num_colors as i32)).collect_vec();
        self.letter_colors.push(new_positions);
    }

    // main logic in the module: convert positional homologies to "global" printtable homologies
    pub fn transl(self) -> FormattedHomologies {
        let k = self.num_colors;
        let mut front_paddings = vec![0u32; k + 1];
        let seq_lengths = self.letter_colors.iter().map(|v| v.len()).collect_vec();
        let mut is_singletons = self
            .letter_colors
            .iter()
            .map(|h| FixedBitSet::with_capacity(h.len()))
            .collect_vec();
        for (i, hits) in self.letter_colors.iter().enumerate() {
            let mut num_singletons_in_front = 0u32;
            let mut first_hit = true; // true when we have only seen one non-negative value
            let is_singleton = is_singletons.get_mut(i).unwrap();
            for (j, &h) in hits.iter().enumerate() {
                if h < 0 {
                    num_singletons_in_front += 1;
                    is_singleton.set(j, true);
                } else {
                    if first_hit {
                        front_paddings[0] = num_singletons_in_front.max(front_paddings[0]);
                        first_hit = false;
                    } else {
                        // padding_needed = padding_needed.min((h as u32) - last_hit);
                        front_paddings[h as usize] =
                            num_singletons_in_front.max(front_paddings[h as usize]);
                    }
                    num_singletons_in_front = 0;
                }
            }
            if num_singletons_in_front > 0 {
                front_paddings[k] = num_singletons_in_front.max(front_paddings[k]);
            }
        }
        let expanded_num_cols = front_paddings.iter().sum::<u32>() as usize + k;
        let mut shifted_columns = vec![0u32; k];
        for c in 0..k {
            if c == 0 {
                shifted_columns[c] = front_paddings[c];
            } else {
                shifted_columns[c] = shifted_columns[c - 1] + front_paddings[c] + 1;
            }
        }
        let positions: Vec<Vec<u32>> = self
            .letter_colors
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
                        let starting_pos = (expanded_num_cols - num_backpositions) as u32;
                        for i in 0..num_backpositions {
                            local_positions[seq_len - num_backpositions + i] =
                                starting_pos + (i) as u32;
                        }
                    }
                }
                // push all the front singleton positions to the "true" front
                local_positions
                    .iter_mut()
                    .enumerate()
                    .take_while(|(j, _p)| is_singletons[i][*j])
                    .for_each(|(j, p)| {
                        *p = j as u32;
                    });
                local_positions
            })
            .collect();
        FormattedHomologies {
            num_visual_columns: expanded_num_cols,
            singleton_letters: is_singletons,
            letter_positions: positions,
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
        let base_pos_map = self.letter_positions.last().unwrap();
        while let Some(record_iffy) = reader.next() {
            let record = record_iffy?;
            w.write_all(b">")?;
            w.write_all(record.head())?;
            w.write_all(b"\n")?;
            let mut buf = vec![b'-'; self.num_visual_columns];
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
            let mut buf: Vec<u8> = vec![b'-'; self.num_visual_columns];
            for (j, &c) in q.seq.iter().enumerate() {
                let target_pos = self.letter_positions[i][j];
                buf[target_pos as usize] = if self.singleton_letters[i][j] {
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
