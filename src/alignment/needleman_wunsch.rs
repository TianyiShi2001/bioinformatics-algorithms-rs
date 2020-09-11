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
#![allow(non_snake_case)]
use super::GlobalAlign;
use super::{Alignment, AlignmentMode, AlignmentOperation, MatchFn, Matrix, Score, Seq};

/// Needleman-Wunsch

pub struct Aligner<'a> {
    x: Seq<'a>,
    y: Seq<'a>,
    S: Matrix<Score>,
    T: Matrix<AlignmentOperation>,
    match_fn: MatchFn,
    gap_penalty: Score,
}

impl<'a> GlobalAlign<'a> for Aligner<'a> {
    fn global(mut self) -> Alignment<'a> {
        let (m, n) = (self.x.len(), self.y.len());
        self.fill_matrices();
        let operations = self.traceback();
        Alignment {
            x: self.x,
            y: self.y,
            xstart: 0,
            ystart: 0,
            xend: m,
            yend: n,
            operations: operations,
            mode: AlignmentMode::Global,
        }
    }
}

impl<'a> Aligner<'a> {
    pub fn new(x: Seq<'a>, y: Seq<'a>, match_fn: MatchFn, gap_penalty: Score) -> Self {
        let (m, n) = (x.len(), y.len());
        let (S, T) = Self::init_matrices(m, n, gap_penalty);
        Aligner {
            x,
            y,
            S,
            T,
            match_fn,
            gap_penalty,
        }
    }
    fn init_matrices(
        m: usize,
        n: usize,
        gap_penalty: Score,
    ) -> (Matrix<Score>, Matrix<AlignmentOperation>) {
        let mut S = vec![vec![0 as Score; n + 1]; m + 1];
        let mut T = vec![vec![AlignmentOperation::Origin; n + 1]; m + 1];
        for j in 1..=n {
            S[0][j] = gap_penalty * j as Score;
            T[0][j] = AlignmentOperation::Insert;
        }
        for i in 1..=m {
            S[i][0] = gap_penalty * i as Score;
            T[i][0] = AlignmentOperation::Delete;
        }
        (S, T)
    }
    fn fill_matrices(&mut self) {
        for i in 1..=self.x.len() {
            for j in 1..=self.y.len() {
                let (xi, yj) = (self.x[i - 1], self.y[j - 1]);
                let diag = self.S[i - 1][j - 1] + (self.match_fn)(xi, yj);
                let up = self.S[i - 1][j] + self.gap_penalty;
                let left = self.S[i][j - 1] + self.gap_penalty;
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
                self.S[i][j] = max_score;
                self.T[i][j] = operation;
            }
        }
    }
    fn traceback(&self) -> Vec<AlignmentOperation> {
        let (mut i, mut j) = (self.x.len(), self.y.len());
        let mut operations = Vec::<AlignmentOperation>::new();
        loop {
            let o = self.T[i][j];
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
