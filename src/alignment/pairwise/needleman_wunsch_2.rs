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
//! # Improving the performace of matrix operations
//!
//! In our first implementation of the Needleman-Wunsch algorithm, we used a `Vec` of `Vec` of `Score`,
//! i.e. `Vec<Vec<Score>>`, to represent a matrix of `Score`. It was straightforward to use, but how about
//! its performance?
//!
//! ## Introducing `cargo bench`
//!
//! `cargo` is shipped along with a benchmark utility. To use it, you write a set of benchmark programs in
//! the `benches/` directory, and `cargo bench` them, without any additional requirements.
//!
//! ### Write your first benchmark program
//!
//! In `benches/matrix.rs`, write:
//!
//! ```ignore
//! #![feature(test)]
//! extern crate test;
//! use test::Bencher;
//! #[bench]
//! fn bench_vec_of_vec(b: &mut Bencher) {
//!     let (m, n) = (1000, 10000);
//!     let mut matrix = vec![vec![0; n]; m];
//!     b.iter(|| {
//!         for i in 0..m {
//!             for j in 0..n {
//!                 matrix[i][j] = i * j;
//!             }
//!         }
//!     });
//! }
//! ```
//!
//! The first three lines are boilerplate code needed for each file containing benchmark functions.
//!
//! Benchmark functions are similar to test functions. While test functions are marked with the attribute
//! `test` and are named `test_*`, benchmark function are marked with the attribute`bench`, and it is
//! idiomatic to begin the name of function with `bench_`.
//!
//! ## How can a matrix be represented in another way?
//!
//! Representing a matrix as a rectangle is intuitive to humans, but for computers, this adds complexity.
//!
//! In the Needleman-Wunsch algorithm, we only need to interact with matrices by getting at setting values
//! at given coordinates `(i, j)`, and in fact we can implement these operations in a 1D vector.
//!
//! Imagine the elements of a `3x7` matrix `M` is numbered from left to right and then from top to bottom row by row,
//! effectively collapsing the matrix into a vector `V`:
//!
//! ```ignore
//!   |   0  1  2  3  4  5  6
//! –––––––––––––––––––––––––
//! 0 | [ 0  1  2  3  4  5  6
//! 1 |   7  8  9 10 11 12 13
//! 2 |  14 15 16 17 18 19 20 ]
//! ```
//!
//! Clealy, getting `M[i][j]` can be translated to getting `V[i * n + j]`, where `n` is the number of columns.
//!
//! ## Comparing the performace of the two approaches
//!
//! # Conclusion
#![allow(non_snake_case)]
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
        let (mut S, mut T) = self.init_matrices(m, n);
        self.fill_matrices(&mut S, &mut T, x, y);
        let operations = self.traceback(&T);
        Alignment {
            x,
            y,
            score: S.get(m, n),
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
