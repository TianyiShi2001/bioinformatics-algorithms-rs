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

pub mod needleman_wunsch;

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
    Origin,
}

#[derive(Debug)]
pub enum AlignmentMode {
    Global,
    Semiglobal,
    Local,
}

#[derive(Debug)]
pub struct Alignment<'a> {
    pub x: Seq<'a>,
    pub y: Seq<'a>,
    pub xstart: usize,
    pub ystart: usize,
    pub xend: usize,
    pub yend: usize,
    pub operations: Vec<AlignmentOperation>,
    pub mode: AlignmentMode,
}

pub trait GlobalAlign<'a> {
    fn global(self) -> Alignment<'a>;
}

pub trait SemiglobalAlign<'a> {
    fn semiglobal(self) -> Alignment<'a>;
}

pub trait LocalAlign<'a> {
    fn local(self) -> Alignment<'a>;
}

type MatchFn = fn(u8, u8) -> Score;

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
            let mut x_i: usize;
            let mut y_i: usize;

            // If the alignment mode is one of the standard ones, the prefix clipping is
            // implicit so we need to process it here
            match self.mode {
                // AlignmentMode::Custom => {
                //     x_i = 0;
                //     y_i = 0;
                // }
                _ => {
                    x_i = self.xstart;
                    y_i = self.ystart;
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
                }
            }

            // Process the alignment.
            for i in 0..self.operations.len() {
                match self.operations[i] {
                    AlignmentOperation::Match => {
                        x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[x[x_i]])));
                        x_i += 1;

                        inb_pretty.push_str("|");

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
                    AlignmentOperation::Delete => {
                        x_pretty.push('-');

                        inb_pretty.push('x');

                        y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[y[y_i]])));
                        y_i += 1;
                    }
                    AlignmentOperation::Insert => {
                        x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[x[x_i]])));
                        x_i += 1;

                        inb_pretty.push('+');

                        y_pretty.push('-');
                    }
                    _ => {}
                }
            }

            // If the alignment mode is one of the standard ones, the suffix clipping is
            // implicit so we need to process it here
            match self.mode {
                // AlignmentMode::Custom => {}
                _ => {
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
            s.push_str("\n");

            s.push_str(&inb_pretty[rng.clone()]);
            s.push_str("\n");

            s.push_str(&y_pretty[rng]);
            s.push_str("\n");

            s.push_str("\n\n");
            idx += step;
        }

        write!(f, "{}", s)
    }
}
