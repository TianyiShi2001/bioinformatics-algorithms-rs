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

//! # Smith-Waterman Local Alignment (1)
//!
//! An introduction to the Needleman-Wunsch algorithm and a straightforward implementation
//!
//! # Learning Outcomes
//!
//! - The principle of the Smith-Waterman local alignment algorithm, its time and space complexity
#![allow(non_snake_case)]
use super::LocalAlign;
use super::{Alignment, AlignmentMode, AlignmentOperation, MatchFn, Score, Seq};
use crate::utils::matrix::Matrix;

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

impl LocalAlign for Aligner {
    /// `S`: a matrix containing alignment scores. `S[i][j]` is the best alignment score between `x[0..i]` and `y[..j]`.
    /// `T`: a matrix containing alignment operations (or "directions", in the context of a matrix).
    fn local<'a>(&self, x: Seq<'a>, y: Seq<'a>) -> Alignment<'a> {
        let (m, n) = (x.len(), y.len());
        let (mut S, mut T) = self.init_matrices(m, n);
        self.fill_matrices(&mut S, &mut T, x, y);
        let (score, (xend, yend)) = Self::find_max_score_and_coords(&S);
        let (operations, (xstart, ystart)) = self.traceback(&T, xend, yend);
        Alignment {
            x,
            y,
            score,
            xstart,
            ystart,
            xend,
            yend,
            operations,
            mode: AlignmentMode::Local,
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
        let S = Matrix::fill(m + 1, n + 1, 0 as Score);
        let T = Matrix::fill(m + 1, n + 1, AlignmentOperation::Origin);
        // we don't need to fill the first column & row
        //
        // for j in 1..=n {
        //     S.set(0, j, self.gap_penalty * j as Score);
        //     T.set(0, j, AlignmentOperation::Insert);
        // }
        // for i in 1..=m {
        //     S.set(i, 0, self.gap_penalty * i as Score);
        //     T.set(i, 0, AlignmentOperation::Delete);
        // }
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
                // if the max score is less than 0, set this point to origin
                if max_score < 0 {
                    max_score = 0;
                    operation = AlignmentOperation::Origin;
                }
                S.set(i, j, max_score);
                T.set(i, j, operation);
            }
        }
    }
    pub fn traceback(
        &self,
        T: &Matrix<AlignmentOperation>,
        mut i: usize,
        mut j: usize,
    ) -> (Vec<AlignmentOperation>, (usize, usize)) {
        // trackback from the max score
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
        (operations, (i, j))
    }
    pub fn find_max_score_and_coords(S: &Matrix<Score>) -> (Score, (usize, usize)) {
        let mut max_score = Score::MIN;
        let mut max_coords = (0usize, 0usize);
        for i in 0..S.nrow {
            for j in 0..S.ncol {
                let v = S.get(i, j);
                if v >= max_score {
                    max_score = v;
                    max_coords = (i, j);
                }
            }
        }
        (max_score, max_coords)
    }
}
