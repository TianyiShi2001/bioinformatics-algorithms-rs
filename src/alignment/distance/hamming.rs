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

use crate::utils::Seq;

pub fn hamming_naive(x: Seq, y: Seq) -> u32 {
    assert_eq!(
        x.len(),
        y.len(),
        "The lengths of the two sequences must equal."
    );
    let mut dist = 0;
    for (p, q) in x.iter().zip(y) {
        if p != q {
            dist += 1;
        }
    }
    dist
}

pub fn hamming_unsafe(x: Seq, y: Seq) -> u32 {
    let len = x.len();
    assert_eq!(len, y.len(), "The lengths of the two sequences must equal.");
    let mut dist = 0;
    for i in 0..len {
        if unsafe { x.get_unchecked(i) != y.get_unchecked(i) } {
            dist += 1
        }
    }
    dist
}
