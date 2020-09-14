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

//! # Needleman-Wunsch Algorithm (2)
//!
//! Re-implement the Needleman-Wunsch algorithm to improve performance
//!
//! # Learning Outcomes
//!
//! - Benchmarking in Rust
//! - Strategies to improve the performace of a program
//!
#![allow(non_snake_case)]
use super::GlobalAlign;
use super::{Alignment, AlignmentMode, AlignmentOperation, MatchFn, Score, Seq};
use crate::utils::matrix::Matrix;

/// Needleman-Wunsch Aligner.
///
pub struct Aligner {
    pub match_fn: MatchFn,
    pub gap_penalty: Score,
}

impl GlobalAlign for Aligner {
    fn global<'a>(&self, x: Seq<'a>, y: Seq<'a>) -> Alignment<'a> {
        let (m, n) = (x.len(), y.len());
        let (mut S, mut T) = self.init_matrices(m, n);
        self.fill_matrices(&mut S, &mut T, x, y);
        let operations = self.traceback(&T);
        Alignment {
            x,
            y,
            xstart: 0,
            ystart: 0,
            xend: m,
            yend: n,
            operations: operations,
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
        let mut S = Matrix::fill(m + 1, n + 1, 0 as Score);
        let mut T = Matrix::fill(m + 1, n + 1, AlignmentOperation::Origin);
        for j in 1..=n {
            S.set(0, j, self.gap_penalty * j as Score);
            T.set(0, j, AlignmentOperation::Insert);
        }
        for i in 1..=m {
            S.set(i, 0, self.gap_penalty * i as Score);
            T.set(i, 0, AlignmentOperation::Delete);
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
                let diag = S.get(i - 1, j - 1) + (self.match_fn)(xi, yj);
                let up = S.get(i - 1, j) + self.gap_penalty;
                let left = S.get(i, j - 1) + self.gap_penalty;
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
                S.set(i, j, max_score);
                T.set(i, j, operation);
            }
        }
    }
    pub fn traceback(&self, T: &Matrix<AlignmentOperation>) -> Vec<AlignmentOperation> {
        let (mut i, mut j) = (T.nrow, T.ncol);
        let mut operations = Vec::<AlignmentOperation>::new();
        loop {
            let o = T.get(i, j);
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
