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
use super::LocalAlign;
use super::{Alignment, AlignmentMode, AlignmentOperation, MatchFunc, Score, Seq, MatchParams};
use crate::utils::matrix::Matrix;

/// Gotoh Aligner.
///
/// # Fields
///
/// - `match_fn`: a function pointer to a function, which takes two characters as `u8` as arguments and computes
///    the substitution score between them.
/// - `gap_penalty`: the penalty for opening a gap; should be a negative integer
pub struct Aligner<F: MatchFunc> {
    pub match_fn: F,
    pub gap_open: Score,
    pub gap_extend: Score,
    S: Matrix<Score>,
    D: Matrix<Score>,
    I: Matrix<Score>,
    TS: Matrix<GotohTraceback>,
    TD: Matrix<GotohTraceback>,
    TI: Matrix<GotohTraceback>,
}

enum GotohTraceback {
    Match,
    Mismatch,
    DeleteFromS,
    DeleteFromD,
    InsertFromS,
    InsertFromI,
}
