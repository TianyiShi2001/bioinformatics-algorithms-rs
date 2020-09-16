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

//! Substitution matrices taken from [SeqAn](https://github.com/seqan/seqan/blob/master/include%2Fseqan%2Fscore%2Fscore_matrix_data.h)
//!
//! Note these special characters in the alphabet:
//!
//! | character | 3-letter code |               Definition                |
//! | :-------: | :-----------: | :-------------------------------------: |
//! |     B     |      Asx      | Asparagine or Aspartic acid (Aspartate) |
//! |     Z     |      Glx      | Glutamine or Glutamic acid (Glutamate)  |
//! |     X     |      Xaa      |        Any amino acid	All codons        |
//! |     *     |      END      |  Termination codon (translation stop)   |
//!
//! # References
//!
//! - https://www.mathworks.com/help/bioinfo/ref/aminolookup.html

pub mod blosum62;

pub use blosum62::blosum62;

#[inline]
fn lookup(a: u8) -> usize {
    if a == b'*' {
        26 as usize
    } else {
        (a - 65) as usize
    }
}
