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

//! The Pairwise Alignment Problem
//pub mod myers_miller;
pub mod myers_miller;
pub mod needleman_wunsch_1;
pub mod needleman_wunsch_2;
pub mod needleman_wunsch_3;
pub mod smith_waterman_1;

use std::fmt;

type Score = i32;

/// Two dimensional matrix
pub type Matrix<T> = Vec<Vec<T>>;
pub type Seq<'a> = &'a [u8];

#[derive(Clone, Copy, Eq, PartialEq, Debug)]
pub enum AlignmentOperation {
    Insert,
    Delete,
    Mismatch,
    Match,
}

#[derive(Debug, PartialEq, Eq)]
pub enum AlignmentMode {
    Global,
    Semiglobal,
    Local,
}

#[derive(Debug)]
pub struct Alignment<'a> {
    pub x: Seq<'a>,
    pub y: Seq<'a>,
    pub score: Score,
    pub xstart: usize,
    pub ystart: usize,
    pub xend: usize,
    pub yend: usize,
    pub operations: Vec<AlignmentOperation>,
    pub mode: AlignmentMode,
}

impl<'a> Alignment<'a> {
    pub fn new(x: Seq<'a>, y: Seq<'a>) -> Self {
        Self {
            x,
            y,
            score: 0,
            xstart: 0,
            ystart: 0,
            xend: x.len(),
            yend: y.len(),
            operations: Vec::new(),
            mode: AlignmentMode::Global,
        }
    }

    pub fn global(x: Seq<'a>, y: Seq<'a>) -> Self {
        Alignment::new(x, y).mode(AlignmentMode::Global)
    }

    pub fn local(x: Seq<'a>, y: Seq<'a>) -> Self {
        Alignment::new(x, y).mode(AlignmentMode::Local)
    }

    pub fn semiglobal(x: Seq<'a>, y: Seq<'a>) -> Self {
        Alignment::new(x, y).mode(AlignmentMode::Semiglobal)
    }

    pub fn mode(mut self, mode: AlignmentMode) -> Self {
        self.mode = mode;
        self
    }

    pub fn score(mut self, score: Score) -> Self {
        self.score = score;
        self
    }

    pub fn operations(mut self, ops: Vec<AlignmentOperation>) -> Self {
        self.operations = ops;
        self
    }

    pub fn xstart(mut self, p: usize) -> Self {
        self.xstart = p;
        self
    }

    pub fn xend(mut self, p: usize) -> Self {
        self.xend = p;
        self
    }

    pub fn ystart(mut self, p: usize) -> Self {
        self.ystart = p;
        self
    }

    pub fn yend(mut self, p: usize) -> Self {
        self.xend = p;
        self
    }

    pub fn start(self, x: usize, y: usize) -> Self {
        self.xstart(x).ystart(y)
    }

    pub fn end(self, x: usize, y: usize) -> Self {
        self.xend(x).yend(y)
    }

    pub fn termini(self, start: (usize, usize), end: (usize, usize)) -> Self {
        self.start(start.0, start.1).end(end.0, end.1)
    }
}

impl PartialEq for Alignment<'_> {
    /// Two alignments are equivalent if their
    ///
    /// - sequences
    /// - score
    /// - start/end positions on both sequences
    /// - alignment mode
    ///
    /// are equal, and the number of inserts, deletes, matches, mismatches each are equal
    fn eq(&self, other: &Self) -> bool {
        if !(
            // *self.x == *other.x
            // && *self.y == *other.y &&
            self.score == other.score
                && self.xstart == other.xstart
                && self.xend == other.xend
                && self.ystart == other.ystart
                && self.yend == other.yend
                && self.mode == other.mode
        ) {
            return false;
        }
        let (o1, o2) = (&self.operations, &other.operations);
        if o1.len() != o2.len() {
            return false;
        }
        let (mut i, mut d, mut s, mut m) = (0usize, 0usize, 0usize, 0usize);
        for o in o1 {
            match o {
                AlignmentOperation::Insert => i += 1,
                AlignmentOperation::Delete => d += 1,
                AlignmentOperation::Mismatch => s += 1,
                AlignmentOperation::Match => m += 1,
            }
        }
        for o in o2 {
            match o {
                AlignmentOperation::Insert => i -= 1,
                AlignmentOperation::Delete => d -= 1,
                AlignmentOperation::Mismatch => s -= 1,
                AlignmentOperation::Match => m -= 1,
            }
        }
        i == 0 && d == 0 && s == 0 && m == 0
    }
}

impl Eq for Alignment<'_> {}

pub trait GlobalAlign {
    fn global<'a>(&self, x: Seq<'a>, y: Seq<'a>) -> Alignment<'a>;
}

pub trait SemiglobalAlign {
    fn semiglobal<'a>(&self, x: Seq<'a>, y: Seq<'a>) -> Alignment<'a>;
}

pub trait LocalAlign {
    fn local<'a>(&self, x: Seq<'a>, y: Seq<'a>) -> Alignment<'a>;
}

impl<'a> fmt::Display for Alignment<'a> {
    /// Adapted from bio_types::alignment according to the MIT License
    /// Copyright 2014-2015 Johannes Köster, Vadim Nazarov, Patrick Marks
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let x = self.x;
        let y = self.y;
        let xlen = self.xend - self.xstart;
        let ylen = self.yend - self.ystart;
        let mut x_pretty = String::new();
        let mut y_pretty = String::new();
        let mut inb_pretty = String::new();

        if !self.operations.is_empty() {
            let mut x_i = self.xstart;
            let mut y_i = self.ystart;

            for k in x.iter().take(self.xstart) {
                x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[*k])));
                inb_pretty.push(' ');
                y_pretty.push(' ')
            }
            for k in y.iter().take(self.ystart) {
                y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[*k])));
                inb_pretty.push(' ');
                x_pretty.push(' ')
            }

            // Process the alignment.
            for i in 0..self.operations.len() {
                match self.operations[i] {
                    AlignmentOperation::Match => {
                        x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[x[x_i]])));
                        x_i += 1;

                        inb_pretty.push('|');

                        y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[y[y_i]])));
                        y_i += 1;
                    }
                    AlignmentOperation::Mismatch => {
                        x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[x[x_i]])));
                        x_i += 1;

                        inb_pretty.push('\\');

                        y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[y[y_i]])));
                        y_i += 1;
                    }
                    AlignmentOperation::Insert => {
                        x_pretty.push('-');

                        inb_pretty.push('x');

                        y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[y[y_i]])));
                        y_i += 1;
                    }
                    AlignmentOperation::Delete => {
                        x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[x[x_i]])));
                        x_i += 1;

                        inb_pretty.push('+');

                        y_pretty.push('-');
                    }
                }
            }

            for k in x.iter().take(xlen).skip(x_i) {
                x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[*k])));
                inb_pretty.push(' ');
                y_pretty.push(' ')
            }
            for k in y.iter().take(ylen).skip(y_i) {
                y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[*k])));
                inb_pretty.push(' ');
                x_pretty.push(' ')
            }
        }

        let mut s = String::new();
        let mut idx = 0;
        let step = 100; // Number of characters per line
        use std::cmp::min;

        assert_eq!(x_pretty.len(), inb_pretty.len());
        assert_eq!(y_pretty.len(), inb_pretty.len());

        let ml = x_pretty.len();

        while idx < ml {
            let rng = idx..min(idx + step, ml);
            s.push_str(&x_pretty[rng.clone()]);
            s.push('\n');

            s.push_str(&inb_pretty[rng.clone()]);
            s.push('\n');

            s.push_str(&y_pretty[rng]);
            s.push('\n');

            s.push_str("\n\n");
            idx += step;
        }

        write!(f, "{}", s)
    }
}

// ! from bio-rust

use std::i32;

/// Value to use as a 'negative infinity' score. Should be close to `i32::MIN`,
/// but avoid underflow when used with reasonable scoring parameters or even
/// adding two negative infinities. Use ~ `0.4 * i32::MIN`
pub const MIN_SCORE: i32 = -858_993_459;

/// Trait required to instantiate a MatchFunc instance
pub trait MatchFunc {
    fn score(&self, a: u8, b: u8) -> i32;
}

/// A concrete data structure which implements trait MatchFunc with constant
/// match and mismatch scores
#[derive(Debug, Clone)]
pub struct MatchParams {
    pub match_score: i32,
    pub mismatch_score: i32,
}

impl MatchParams {
    /// Create new MatchParams instance with given match and mismatch scores
    ///
    /// # Arguments
    ///
    /// * `match_score` - the score for a match (should not be negative)
    /// * `mismatch_score` - the score for a mismatch (should not be positive)
    pub fn new(match_score: i32, mismatch_score: i32) -> Self {
        assert!(match_score >= 0, "match_score can't be negative");
        assert!(mismatch_score <= 0, "mismatch_score can't be positive");
        MatchParams {
            match_score,
            mismatch_score,
        }
    }
}

impl MatchFunc for MatchParams {
    #[inline]
    fn score(&self, a: u8, b: u8) -> i32 {
        if a == b {
            self.match_score
        } else {
            self.mismatch_score
        }
    }
}

/// The trait Matchfunc is also implemented for Fn(u8, u8) -> i32 so that MatchFunc
/// can be instantiated using closures and custom user defined functions
impl<F> MatchFunc for F
where
    F: Fn(u8, u8) -> i32,
{
    fn score(&self, a: u8, b: u8) -> i32 {
        (self)(a, b)
    }
}

/// Details of scoring are encapsulated in this structure.
///
/// An [affine gap score model](https://en.wikipedia.org/wiki/Gap_penalty#Affine)
/// is used so that the gap score for a length `k` is:
/// `GapScore(k) = gap_open + gap_extend * k`
#[derive(Debug, Clone)]
pub struct Scoring<F: MatchFunc> {
    pub gap_open: i32,
    pub gap_extend: i32,
    pub match_fn: F,
}

impl Scoring<MatchParams> {
    /// Create new Scoring instance with given gap open, gap extend penalties
    /// match and mismatch scores. The clip penalties are set to `MIN_SCORE` by default
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should not be positive)
    /// * `gap_extend` - the score for extending a gap (should not be positive)
    /// * `match_score` - the score for a match
    /// * `mismatch_score` - the score for a mismatch
    pub fn from_scores(
        gap_open: i32,
        gap_extend: i32,
        match_score: i32,
        mismatch_score: i32,
    ) -> Self {
        assert!(gap_open <= 0, "gap_open can't be positive");
        assert!(gap_extend <= 0, "gap_extend can't be positive");

        Scoring {
            gap_open,
            gap_extend,
            match_fn: MatchParams::new(match_score, mismatch_score),
        }
    }
}

impl<F: MatchFunc> Scoring<F> {
    /// Create new Scoring instance with given gap open, gap extend penalties
    /// and the score function. The clip penalties are set to [`MIN_SCORE`](constant.MIN_SCORE.html) by default
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should not be positive)
    /// * `gap_extend` - the score for extending a gap (should not be positive)
    /// * `match_fn` - function that returns the score for substitutions
    ///    (see also [`bio::alignment::pairwise::Scoring`](struct.Scoring.html))
    pub fn new(gap_open: i32, gap_extend: i32, match_fn: F) -> Self {
        assert!(gap_open <= 0, "gap_open can't be positive");
        assert!(gap_extend <= 0, "gap_extend can't be positive");

        Scoring {
            gap_open,
            gap_extend,
            match_fn,
        }
    }
}

pub mod traceback_naive {
    use std::fmt;
    #[derive(Copy, Clone)]
    pub enum TracebackDirection {
        Up,
        Left,
        Diag,
        Origin,
    }

    impl fmt::Display for TracebackDirection {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            match self {
                TracebackDirection::Up => write!(f, "↑"),
                TracebackDirection::Left => write!(f, "←"),
                TracebackDirection::Diag => write!(f, "↖"),
                TracebackDirection::Origin => write!(f, "·"),
            }
        }
    }
}

pub mod traceback {
    #[derive(Copy, Clone)]
    pub struct TracebackCell(u8);

    pub const TB_ORIGIN: u8 = 0b00;
    pub const TB_LEFT: u8 = 0b01;
    pub const TB_UP: u8 = 0b10;
    pub const TB_DIAG: u8 = 0b11;

    // Traceback bit positions (LSB)
    const I_POS: u8 = 0; // Meaning bits 0,1 corresponds to I and so on
    const D_POS: u8 = 2;
    const S_POS: u8 = 4;

    impl std::default::Default for TracebackCell {
        #[inline(always)]
        fn default() -> Self {
            Self(0u8)
        }
    }

    impl TracebackCell {
        /// Initialize a blank traceback cell
        pub fn new() -> Self {
            Self::default()
        }

        /// Sets 2 bits [pos, pos+2) with the 2 LSBs of value
        #[inline(always)]
        fn set_bits(&mut self, pos: u8, value: u8) {
            let bits: u8 = (0b11) << pos;
            self.0 = (self.0 & !bits) // First clear the bits
            | (value << pos) // And set the bits
        }

        #[inline(always)]
        pub fn set_i_bits(&mut self, value: u8) {
            // Traceback corresponding to matrix I
            self.set_bits(I_POS, value);
        }

        #[inline(always)]
        pub fn set_d_bits(&mut self, value: u8) {
            // Traceback corresponding to matrix D
            self.set_bits(D_POS, value);
        }

        #[inline(always)]
        pub fn set_s_bits(&mut self, value: u8) {
            // Traceback corresponding to matrix S
            self.set_bits(S_POS, value);
        }

        // Gets 4 bits [pos, pos+4) of v
        #[inline(always)]
        fn get_bits(self, pos: u8) -> u8 {
            (self.0 >> pos) & (0b11)
        }

        #[inline(always)]
        pub fn get_i_bits(self) -> u8 {
            self.get_bits(I_POS)
        }

        #[inline(always)]
        pub fn get_d_bits(self) -> u8 {
            self.get_bits(D_POS)
        }

        #[inline(always)]
        pub fn get_s_bits(self) -> u8 {
            self.get_bits(S_POS)
        }

        /// Set all matrices to the same value.
        pub fn set_all(&mut self, value: u8) {
            self.set_i_bits(value);
            self.set_d_bits(value);
            self.set_s_bits(value);
        }
    }

    impl std::fmt::Debug for TracebackCell {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            f.debug_struct(&format!("{:06b}", self.0)).finish()
        }
    }
}

// pub mod tests {
//     // TODO: macros that help running multiple aligners?
//     use super::myers_miller::Aligner;
//     use super::*;
//     use crate::alignment::pairwise::{Alignment, AlignmentOperation::*};
//     use std::iter::repeat;

//     macro_rules! test_global_alignment {
//         ($aligner:item) => {
//             #[test]
//             fn test_global_affine_ins() {
//                 let x = b"ACGAGAACA";
//                 let y = b"ACGACA";
//                 let score = |a: u8, b: u8| if a == b { 1i32 } else { -3i32 };
//                 let aligner = $aligner::new(-5, -1, score);
//                 let alignment = aligner.global(x, y);

//                 println!("{}", alignment);
//                 assert_eq!(
//                     alignment,
//                     Alignment::global(x, y).score(-2).operations(vec![
//                         Match, Match, Match, Delete, Delete, Delete, Match, Match, Match
//                     ])
//                 );
//             }

//             #[test]
//             fn test_global_affine_ins2() {
//                 let x = b"AGATAGATAGATAGGGAGTTGTGTAGATGATCCACAGT";
//                 let y = b"AGATAGATAGATGTAGATGATCCACAGT";
//                 let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
//                 let aligner = $aligner::new(-5, -1, score);
//                 let alignment = aligner.global(x, y);

//                 println!("{:?}", alignment);
//                 println!("{}", alignment);

//                 let mut correct = Vec::new();
//                 correct.extend(repeat(Match).take(11));
//                 correct.extend(repeat(Delete).take(10));
//                 correct.extend(repeat(Match).take(17));

//                 assert_eq!(
//                     alignment,
//                     Alignment::global(x, y).score(13).operations(correct)
//                 );
//             }

//             #[test]
//             fn test_global() {
//                 let x = b"ACCGTGGAT";
//                 let y = b"AAAAACCGTTGAT";
//                 let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
//                 let aligner = $aligner::new(-5, -1, score);
//                 let alignment = aligner.global(x, y);

//                 println!("{}", alignment);

//                 assert_eq!(
//                     alignment,
//                     Alignment::global(x, y).score(-2).operations(vec![
//                         Insert, Insert, Insert, Insert, Match, Match, Match, Match, Match,
//                         Mismatch, Match, Match, Match,
//                     ])
//                 );
//             }

//             // #[test]
//             // fn test_blosum62() {
//             //     let x = b"AAAA";
//             //     let y = b"AAAA";
//             //     let score = &blosum62;
//             //     let aligner = $aligner::new(-5, -1, score);
//             //     let alignment = aligner.global(x, y);
//             //     assert_eq!(alignment.ystart, 0);
//             //     assert_eq!(alignment.xstart, 0);
//             //     assert_eq!(alignment.score, 16);
//             //     assert_eq!(alignment, [Match, Match, Match, Match]);
//             // }

//             #[test]
//             fn test_issue11() {
//                 let y = b"TACC"; //GTGGAC";
//                 let x = b"AAAAACC"; //GTTGACGCAA";
//                 let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
//                 let aligner = $aligner::new(-5, -1, score);
//                 let alignment = aligner.global(x, y);
//                 assert_eq!(alignment.ystart, 0);
//                 assert_eq!(alignment.xstart, 0);
//                 assert_eq!(
//                     alignment,
//                     Alignment::global(x, y)
//                         .score(-6)
//                         .operations(vec![Delete, Delete, Delete, Mismatch, Match, Match, Match])
//                 );
//             }

//             #[test]
//             fn test_left_aligned_del() {
//                 let x = b"GTGCATCATGTG";
//                 let y = b"GTGCATCATCATGTG";
//                 let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
//                 let aligner = $aligner::new(-5, -1, score);
//                 let alignment = aligner.global(x, y);
//                 println!("{}", alignment);

//                 assert_eq!(
//                     alignment,
//                     Alignment::global(x, y).score(4).operations(vec![
//                         Match, Match, Match, Insert, Insert, Insert, Match, Match, Match, Match,
//                         Match, Match, Match, Match, Match,
//                     ])
//                 )
//             }

//             // Test that trailing deletions are correctly handled
//             // in global mode
//             #[test]
//             fn test_global_right_del() {
//                 let x = b"AACCACGTACGTGGGGGGA";
//                 let y = b"CCACGTACGT";

//                 let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
//                 let aligner = $aligner::new(-5, -1, score);
//                 let alignment = aligner.global(x, y);

//                 println!("{}", alignment);

//                 assert_eq!(
//                     alignment,
//                     Alignment::global(x, y).score(-9).operations(vec![
//                         Delete, Delete, Match, Match, Match, Match, Match, Match, Match, Match,
//                         Match, Match, Delete, Delete, Delete, Delete, Delete, Delete, Delete,
//                     ])
//                 );
//             }

//             #[test]
//             fn test_left_aligned_ins() {
//                 let x = b"GTGCATCATCATGTG";
//                 let y = b"GTGCATCATGTG";
//                 let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
//                 let aligner = $aligner::new(-5, -1, score);
//                 let alignment = aligner.global(x, y);
//                 println!("{}", alignment);

//                 assert_eq!(alignment.ystart, 0);
//                 assert_eq!(alignment.xstart, 0);
//                 assert_eq!(
//                     alignment,
//                     Alignment::global(x, y).score(4).operations(vec![
//                         Match, Match, Match, Delete, Delete, Delete, Match, Match, Match, Match,
//                         Match, Match, Match, Match, Match,
//                     ])
//                 );
//             }
//         };
//     }
//     test_global_alignment!(myers_miller::Aligner);

//     macro_rules! test_semiglobal_alignment {
//         ($aligner:item) => {
//             #[test]
//             fn test_semiglobal() {
//                 let x = b"AAAAACCGTTGAT";
//                 let y = b"ACCGTGGAT";
//                 let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
//                 let aligner = $aligner::new(-5, -1, score);
//                 let alignment = aligner.semiglobal(x, y);
//                 println!("{}", alignment);
//                 assert_eq!(
//                     alignment,
//                     Alignment::semiglobal(x, y)
//                         .score(7)
//                         .xstart(4)
//                         .operations(vec![
//                             Match, Match, Match, Match, Match, Mismatch, Match, Match, Match,
//                         ])
//                 );
//             }
//         };
//     }

//     test_semiglobal_alignment!(myers_miller::Aligner);

//     macro_rules! test_local_alignment {
//         ($aligner:item) => {
//             #[test]
//             fn test_local_affine_ins2() {
//                 let x = b"ACGTATCATAGATAGATAGGGTTGTGTAGATGATCCACAG";
//                 let y = b"CGTATCATAGATAGATGTAGATGATCCACAGT";
//                 let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
//                 let aligner = $aligner::new(-5, -1, score);
//                 let alignment = aligner.local(x, y);
//                 assert_eq!(alignment.xstart, 1);
//                 assert_eq!(alignment.ystart, 0);
//             }

//             #[test]
//             fn test_local() {
//                 let x = b"ACCGTGGAT";
//                 let y = b"AAAAACCGTTGAT";
//                 let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
//                 let aligner = $aligner::new(-5, -1, score);
//                 let alignment = aligner.local(x, y);
//                 println!("{}", alignment);
//                 assert_eq!(
//                     alignment,
//                     Alignment::local(x, y)
//                         .xstart(0)
//                         .ystart(4)
//                         .score(7)
//                         .operations(vec![
//                             Match, Match, Match, Match, Match, Mismatch, Match, Match, Match,
//                         ])
//                 );
//             }
//         };
//     }

//     test_local_alignment!(myers_miller::Aligner);
// }
