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

//! # Needleman-Wunsch Algorithm (3)
//!
//! Re-implement the Needleman-Wunsch algorithm to improve space efficiency
//!
//! # Learning Outcomes
//!
//! - Strategies to improve space efficiency
//!
//! # Improving Space Efficiency
#![allow(non_snake_case)]
#![allow(clippy::many_single_char_names)]
use super::traceback_naive::TracebackDirection;
use super::GlobalAlign;
use super::{Alignment, AlignmentMode, AlignmentOperation, Score, Seq};
use crate::utils::matrix::Matrix;
type MatchFn = fn(u8, u8) -> Score;

/// Needleman-Wunsch Aligner.
pub struct Aligner {
    pub match_fn: MatchFn,
    pub gap_penalty: Score,
}

impl GlobalAlign for Aligner {
    fn global<'a>(&self, x: Seq<'a>, y: Seq<'a>) -> Alignment<'a> {
        let (m, n) = (x.len(), y.len());
        let (S, T) = self.fill_matrices(x, y);
        let operations = self.traceback(&T, x, y);
        Alignment {
            x,
            y,
            score: S[n],
            xstart: 0,
            ystart: 0,
            xend: m,
            yend: n,
            operations,
            mode: AlignmentMode::Global,
        }
    }
}

impl Aligner {
    pub fn new(match_fn: MatchFn, gap_penalty: Score) -> Self {
        Aligner {
            match_fn,
            gap_penalty,
        }
    }

    #[allow(clippy::needless_range_loop)]
    pub fn fill_matrices(&self, x: Seq, y: Seq) -> (Vec<Score>, Matrix<TracebackDirection>) {
        let (m, n) = (x.len(), y.len());
        // init matrices
        let mut S = vec![0; n + 1];
        let mut T = Matrix::fill(m + 1, n + 1, TracebackDirection::Origin);
        let mut s: Score; // S[i - 1][j - 1]
        let mut c: Score; // S[i][j - 1]
                          // fill the first row (row 0)
        for j in 1..=n {
            S[j] = self.gap_penalty * j as Score;
            T.set(0, j, TracebackDirection::Left);
        }
        // for each row i from row 1
        for i in 1..=m {
            // setting s to previous S[0] i.e. S[i - 1][0]
            s = S[0];
            // current S[0] (S[i][0])
            c = self.gap_penalty + i as Score;
            // update S[i][0] and T[i][0] (the first element of each row)
            S[0] = c; // note that S[0] really is the 0-th element of the i-th row, i.e. S[i][0]; before updating, it was S[i - 1][0]
            T.set(i, 0, TracebackDirection::Up);
            let xi = x[i - 1];
            // xi only needs to be calculated once for each row; so instead of writing `let (xi, yj) = (x[i - 1], y[j - 1])`, separating them is better
            for j in 1..=n {
                let yj = y[j - 1];

                let diag = s + (self.match_fn)(xi, yj);
                // previously: `let diag = S.get(i - 1, j - 1) + (self.match_fn)(xi, yj);`

                let up = S[j] + self.gap_penalty;
                // previously: let up = S.get(i - 1, j) + self.gap_penalty;

                let left = c + self.gap_penalty;
                // previously: let left = S.get(i, j - 1) + self.gap_penalty;

                let mut max_score = diag;
                let mut direction = TracebackDirection::Diag;
                if up > max_score {
                    max_score = up;
                    direction = TracebackDirection::Up;
                }
                if left > max_score {
                    max_score = left;
                    direction = TracebackDirection::Left;
                }
                s = S[j]; // the value of S[j] before updating will become S[i - 1][j - 1] in the next round of calculation
                S[j] = max_score;
                c = max_score; // the value of S[j] after updating will become S[i][j - 1] in the next round of calculation
                T.set(i, j, direction);
            }
        }
        (S, T)
    }
    pub fn traceback(
        &self,
        T: &Matrix<TracebackDirection>,
        x: Seq,
        y: Seq,
    ) -> Vec<AlignmentOperation> {
        let (mut i, mut j) = (T.nrow, T.ncol);
        let mut operations = Vec::<AlignmentOperation>::new();
        loop {
            match T.get(i, j) {
                TracebackDirection::Diag => {
                    operations.push(if x[i - 1] == y[j - 1] {
                        AlignmentOperation::Match
                    } else {
                        AlignmentOperation::Mismatch
                    });
                    i -= 1; // moving diagonal
                    j -= 1;
                }
                TracebackDirection::Up => {
                    operations.push(AlignmentOperation::Delete);
                    i -= 1; // moving up
                }
                TracebackDirection::Left => {
                    operations.push(AlignmentOperation::Insert);
                    j -= 1; // moving left
                }
                TracebackDirection::Origin => break, // reaching origin (0,0)
            }
        }
        operations
    }
}
