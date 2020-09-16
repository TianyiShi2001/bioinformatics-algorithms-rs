// Copyright (C) 2020 Tianyi Shi
//
// This file is part of rust-bio-edu.
//
// rust-bio-edu is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// rust-bio-edu is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with rust-bio-edu.  If not, see <http://www.gnu.org/licenses/>.

//! # Needleman-Wunsch Algorithm (1)
//!
//! An introduction to the Needleman-Wunsch algorithm and a straightforward implementation
//!
//! # Learning Outcomes
//!
//! - The principle of the Needleman-Wunsch algorithm, its time and space complexity
//! - Familiarize with Rust Programming: lifetimes, traits, borrowing, etc.
//!
//! Original Literature: [Needleman SB, Wunsch CD. A general method applicable to the search for similarities in the amino acid sequence of two proteins. J Mol Biol. 1970;48(3):443-453.](https://doi.org/10.1016/0022-2836(70)90057-4)
//!
#![allow(non_snake_case)]
use super::GlobalAlign;
use super::{Alignment, AlignmentMode, AlignmentOperation, Score, Seq};
type Matrix<T> = Vec<Vec<T>>;
type MatchFn = fn(u8, u8) -> Score;

/// Needleman-Wunsch Aligner.
///
/// # Fields
///
/// - `match_fn`: a function pointer to a function, which takes two characters as `u8` as arguments and computes
///    the substitution score between them.
/// - `gap_penalty`: the penalty for opening a gap; should be a negative integer
pub struct Aligner {
    pub match_fn: MatchFn,
    pub gap_penalty: Score,
}

impl GlobalAlign for Aligner {
    /// `S`: a matrix containing alignment scores. `S[i][j]` is the best alignment score between `x[0..i]` and `y[..j]`.
    /// `T`: a matrix containing alignment operations (or "directions", in the context of a matrix).
    fn global<'a>(&self, x: Seq<'a>, y: Seq<'a>) -> Alignment<'a> {
        let (m, n) = (x.len(), y.len());
        let (mut S, mut T) = self.init_matrices(m, n);
        self.fill_matrices(&mut S, &mut T, x, y);
        let operations = self.traceback(&T);
        Alignment {
            x,
            y,
            score: S[m][n],
            xstart: 0,
            ystart: 0,
            xend: m,
            yend: n,
            operations,
            mode: AlignmentMode::Global,
        }
    }
}

impl<'a> Aligner {
    pub fn new(match_fn: MatchFn, gap_penalty: Score) -> Self {
        Aligner {
            match_fn,
            gap_penalty,
        }
    }
    pub fn init_matrices(&self, m: usize, n: usize) -> (Matrix<Score>, Matrix<AlignmentOperation>) {
        let mut S = vec![vec![0 as Score; n + 1]; m + 1];
        let mut T = vec![vec![AlignmentOperation::Origin; n + 1]; m + 1];
        for j in 1..=n {
            S[0][j] = self.gap_penalty * j as Score;
            T[0][j] = AlignmentOperation::Insert;
        }
        for i in 1..=m {
            S[i][0] = self.gap_penalty * i as Score;
            T[i][0] = AlignmentOperation::Delete;
        }
        (S, T)
    }
    pub fn fill_matrices(
        &self,
        S: &mut Matrix<Score>,
        T: &mut Matrix<AlignmentOperation>,
        x: Seq,
        y: Seq,
    ) {
        for i in 1..=x.len() {
            for j in 1..=y.len() {
                let (xi, yj) = (x[i - 1], y[j - 1]);
                let diag = S[i - 1][j - 1] + (self.match_fn)(xi, yj);
                let up = S[i - 1][j] + self.gap_penalty;
                let left = S[i][j - 1] + self.gap_penalty;
                let mut max_score = diag;
                let mut operation = if xi == yj {
                    AlignmentOperation::Match
                } else {
                    AlignmentOperation::Mismatch
                };
                if up > max_score {
                    max_score = up;
                    operation = AlignmentOperation::Delete;
                }
                if left > max_score {
                    max_score = left;
                    operation = AlignmentOperation::Insert;
                }
                S[i][j] = max_score;
                T[i][j] = operation;
            }
        }
    }
    pub fn traceback(&self, T: &Matrix<AlignmentOperation>) -> Vec<AlignmentOperation> {
        let (mut i, mut j) = (T.len() - 1, T[0].len() - 1);
        let mut operations = Vec::<AlignmentOperation>::new();
        loop {
            let o = T[i][j];
            match o {
                AlignmentOperation::Match | AlignmentOperation::Mismatch => {
                    operations.push(o);
                    i -= 1; // moving diagonal
                    j -= 1;
                }
                AlignmentOperation::Delete => {
                    operations.push(o);
                    i -= 1; // moving up
                }
                AlignmentOperation::Insert => {
                    operations.push(o);
                    j -= 1; // moving left
                }
                AlignmentOperation::Origin => break, // reaching origin (0,0)
            }
        }
        operations
    }
}
