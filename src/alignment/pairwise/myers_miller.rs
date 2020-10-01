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

// Copyright 2020 Tianyi Shi
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Alignment with affine gap penalty in linear space, by combining Gotoh's (1982) and
//! Hirschberg's (1975) ideas, which was first implemented in C (Myers & Miller 1988).
//!
//! Myers & Miller originally used their technique to implement global alignment only,
//! but alignments of other modes can be achieved by first finding the termini of
//! the non-global alignment and then global-aligning the corresponding substrings.
//!
//! # Time Complexity
//!
//! $O(nm)$ for strings of length $m$ and $n$.
//!
//! # Space Complexity
//!
//! $O(n)$. For exact number of bits, see `global`, `semiglobal`, `local`, and `custom` methods
//!
//! # References
//!
//! - [Eugene W. Myers and Webb Miller (1988) Optimal alignments in linear space. _Bioinformatics_ **4**: 11-17.](https://doi.org/10.1093/bioinformatics/4.1.11)
//! - [Hirschberg, D. S. (1975) A linear space algorithm for computing maximal common subsequences. _Commun. Assoc. Comput. Mach._ **18**: 341-343.](https://doi.org/10.1145/360825.360861)
//! - [Gotoh, O. (1982) An improved algorithm for matching biological sequences. _J. Molec. Biol._ **162**: 705-708.](https://doi.org/10.1016/0022-2836(82)90398-9)

#![allow(non_snake_case)]
#![allow(clippy::many_single_char_names)]
use crate::alignment::pairwise::MIN_SCORE;
use crate::alignment::pairwise::{Alignment, AlignmentMode, AlignmentOperation};
use crate::alignment::pairwise::{GlobalAlign, LocalAlign, SemiglobalAlign};
use crate::alignment::pairwise::{MatchFunc, Scoring};
use crate::utils::Seq;
use std::cmp::max;

pub struct Aligner<F: MatchFunc + Sync> {
    scoring: Scoring<F>,
}

impl<F: MatchFunc + Sync> GlobalAlign for Aligner<F> {
    /// Calculate global alignment of `x` against `y`.
    ///
    /// - Time complexity: $O(mn)$, where $m$ and $n$ are the lengths of the first and the second
    ///   sequence
    /// - Space complexity: $O(n)$; specifically, about $64n$ bits. Note that `compute_recursive()`
    ///   uses less and less space as recursion proceeds.
    fn global<'a>(&self, x: Seq<'a>, y: Seq<'a>) -> Alignment<'a> {
        let (m, n) = (x.len(), y.len());
        let operations =
            self.compute_recursive(x, y, m, n, self.scoring.gap_open, self.scoring.gap_open);
        let score = self.cost_only(x, y, false, self.scoring.gap_open).0[y.len()];
        Alignment {
            x,
            y,
            score,
            xstart: 0,
            ystart: 0,
            xend: m,
            yend: n,
            operations,
            mode: AlignmentMode::Global,
        }
    }
}
impl<F: MatchFunc + Sync> SemiglobalAlign for Aligner<F> {
    /// Calculate semiglobal alignment of x against y (x is global, y is local).
    ///
    /// `xstart` is always 0 and `xend` is always `x.len()`. So this algorithm first finds `ystart`
    /// and `yend`, thus determining the coordinates of the termini of the optimal alignment. Then,
    /// a global alignment on the corresponding substrings calculates the operations.
    ///
    /// - Space complexity: $O(n)$
    /// - Time complexity: $(mn + mn')$ where $m$ and $n$ are lengths of the input sequences `x` and
    ///   `y`; $n'$ is the substrings of `y` that correspond to the optimal alignment.
    /// TODO: y against x, not the reverse
    fn semiglobal<'a>(&self, x: Seq<'a>, y: Seq<'a>) -> Alignment<'a> {
        // Compute the alignment
        let (score, xstart, xend) = self.find_semiglobal_score_and_termini(x, y);
        let xlen = xend - xstart;
        let ylen = y.len();
        let operations = self.compute_recursive(
            &x[xstart..xend],
            y,
            xlen,
            ylen,
            self.scoring.gap_open,
            self.scoring.gap_open,
        );

        Alignment {
            x,
            y,
            score,
            xstart,
            xend,
            ystart: 0,
            yend: ylen,
            operations,
            mode: AlignmentMode::Semiglobal,
        }
    }
}
impl<F: MatchFunc + Sync> LocalAlign for Aligner<F> {
    /// Calculate local alignment of x against y.
    ///
    /// This algorithm first find the pair of coordinates corresponding to the termini of one optimal
    /// local alignment path. Then, a global alignment on the corresponding substrings calculates the
    /// operations.
    ///
    /// - Space complexity: $O(n)$
    /// - Time complexity: $(mn + m'n')$ where $m$ and $n$ are lengths of the input sequences; $m'$
    ///   and $n'$ are the lengths of the substrings of the input sequences that correspond to the
    ///   optimal alignment.
    ///
    /// Termini can be determined by two approaches. Huang's method (1991) is faster, but uses about
    /// $446n$ bits. Shamir's method is slightly slower, but only uses $64n$ bits.
    ///
    /// # References
    ///
    /// - [Huang, X. and Miller, W. 1991. A time-efficient linear-space local similarity algorithm. Adv. Appl. Math. 12, 3 (Sep. 1991), 337-357](https://doi.org/10.1016/0196-8858(91)90017-D)
    fn local<'a>(&self, x: Seq<'a>, y: Seq<'a>) -> Alignment<'a> {
        let (score, xstart, ystart, xend, yend) = self.find_local_score_and_termini_huang(x, y);
        let xlen = xend - xstart;
        let ylen = yend - ystart;
        let operations = self.compute_recursive(
            &x[xstart..xend],
            &y[ystart..yend],
            xlen,
            ylen,
            self.scoring.gap_open,
            self.scoring.gap_open,
        );
        Alignment {
            x,
            y,
            score,
            xstart,
            xend,
            ystart,
            yend,
            operations,
            mode: AlignmentMode::Local,
        }
    }
}

impl<F: MatchFunc + Sync> Aligner<F> {
    /// Create new aligner instance with given gap open and gap extend penalties
    /// and the score function.
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should be negative)
    /// * `gap_extend` - the score for extending a gap (should be negative)
    /// * `match_fn` - function that returns the score for substitutions (also see bio::scores)
    pub fn new(gap_open: i32, gap_extend: i32, match_fn: F) -> Self {
        Aligner {
            scoring: Scoring::new(gap_open, gap_extend, match_fn),
        }
    }

    /// Recursively compute alignments of sub-sequences and concatenating them.
    ///
    /// - Space complexity: $O(n)$. Precisely, about $128n + log\_2{m}$, where $m =$ `x.len()` and $n =$ `y.len()`.
    ///   $128n$ is used for four `Vec<i32>` in `find_mid()` and $log\_2{m}$ is due to the activation stack (the
    ///   latter can be negligible).
    ///
    /// # Panics
    ///
    /// Rust has a [`recursion_limit` which defaults to 128](https://doc.rust-lang.org/reference/attributes/limits.html#the-recursion_limit-attribute),
    /// which means the length of the sequence `x` should not exceed $2^128 = 3.4\times10^38$ (well, you'll
    /// never encounter such gigantic biological sequences)
    fn compute_recursive(
        &self,
        x: Seq<'_>,
        y: Seq<'_>,
        m: usize,
        n: usize,
        tb: i32,
        te: i32,
    ) -> Vec<AlignmentOperation> {
        // * m = x.len(); n = y.len()
        if n == 0 {
            return vec![AlignmentOperation::Delete; m];
        }
        if m == 0 {
            return vec![AlignmentOperation::Insert; n];
        }
        if m == 1 {
            return self.nw_onerow(x[0], y, n, max(tb, te));
        }
        let (imid, jmid, join_by_deletion) = self.find_mid(x, y, m, n, tb, te); // `find_mid()` uses 128n bits
        if join_by_deletion {
            // y.len() === (&y[..jmid]).len() + (&y[jmid..]).len()
            // so the sum of the space used by the two subtasks is still 128n (n = y.len())
            let (a, b) = rayon::join(
                || self.compute_recursive(&x[..imid - 1], &y[..jmid], imid - 1, jmid, tb, 0),
                || {
                    self.compute_recursive(
                        &x[imid + 1..],
                        &y[jmid..],
                        m - imid - 1,
                        n - jmid,
                        0,
                        te,
                    )
                },
            );
            [a, vec![AlignmentOperation::Delete; 2], b].concat()
        } else {
            let (a, b) = rayon::join(
                || {
                    self.compute_recursive(
                        &x[..imid],
                        &y[..jmid],
                        imid,
                        jmid,
                        tb,
                        self.scoring.gap_open,
                    )
                },
                || {
                    self.compute_recursive(
                        &x[imid..],
                        &y[jmid..],
                        m - imid,
                        n - jmid,
                        self.scoring.gap_open,
                        te,
                    )
                },
            );
            [a, b].concat()
        }
    }

    /// Find the "midpoint" (see module-level documentation)
    ///
    /// - Space complexity: $O(n)$. Specifically, about $128n$ bits (each of `cc_upper`,
    ///   `dd_upper`, `cc_lower` and `dd_lower` are `Vec<i32>` of length $n$), where $n =$ `y.len()`.
    /// - Time complexity: $O(nm)$
    fn find_mid(
        &self,
        x: Seq<'_>,
        y: Seq<'_>,
        m: usize,
        n: usize,
        tb: i32,
        te: i32,
    ) -> (usize, usize, bool) {
        let imid = m / 2;
        let ((cc_upper, dd_upper), (cc_lower, dd_lower)) = rayon::join(
            || self.cost_only(&x[..imid], y, false, tb),
            || self.cost_only(&x[imid..], y, true, te),
        );
        let mut max = MIN_SCORE;
        let mut jmid = 0;
        let mut join_by_deletion = false;
        for j in 0..=n {
            let c = cc_upper[j] + cc_lower[n - j];
            if c > max {
                max = c;
                jmid = j;
                join_by_deletion = false;
            }
            let d = dd_upper[j] + dd_lower[n - j] - self.scoring.gap_open; // subtract duplicating open!
            if d > max {
                max = d;
                jmid = j;
                join_by_deletion = true;
            }
        }
        (imid, jmid, join_by_deletion)
    }

    /// Cost-only (score-only) Gotoh's algorithm in linear space
    ///
    /// - Space Complexity: $O(n)$; specifically, about $64n$ bits, where $n =$ `y.len() + 1`
    /// Use six scalars and two vectors of length (N + 1), where N is the length
    /// of the shorter sequence.
    /// -Time complexity: $O(nm)$
    fn cost_only(&self, x: Seq, y: Seq, rev: bool, tx: i32) -> (Vec<i32>, Vec<i32>) {
        let (m, n) = (x.len() + 1, y.len() + 1);
        let mut cc: Vec<i32> = vec![0; n]; // match/mismatch    32 * n bits
        let mut dd: Vec<i32> = vec![0; n]; // deletion          32 * n bits
        let mut e: i32; // I(i, j-1)
        let mut c: i32; // C(i, j-1)
        let mut s: i32; // C(i-1, j-1)
        let mut t: i32;
        let mut p: u8;
        let mut q: u8;
        t = self.scoring.gap_open;
        for j in 1..n {
            t += self.scoring.gap_extend;
            cc[j] = t;
            dd[j] = MIN_SCORE;
        }
        t = tx; // originally self.scoring.gap_open;
        for i in 1..m {
            s = cc[0];
            t += self.scoring.gap_extend;
            c = t;
            cc[0] = c;
            // dd[0] = c;
            e = MIN_SCORE;
            p = if rev { x[m - i - 1] } else { x[i - 1] };
            for j in 1..n {
                q = if rev { y[n - j - 1] } else { y[j - 1] };
                e = max(e, c + self.scoring.gap_open) + self.scoring.gap_extend; // update e to I[i,j]
                dd[j] = max(dd[j], cc[j] + self.scoring.gap_open) + self.scoring.gap_extend; // cc[j] = C[i-1, j]
                c = max(max(dd[j], e), s + self.scoring.match_fn.score(p, q));
                s = cc[j];
                cc[j] = c;
            }
        }
        dd[0] = cc[0]; // otherwise deletions at start/end will be free
        (cc, dd)
    }

    /// Find the maximum score in a local alignment between `x` and `y` and the coordinates
    /// of the alignment path termini
    ///
    /// - Space complexity: $O(n)$; specifically, about $448n$ bits (x64 architecture)
    ///   or $256n$ bits (x32 architecture), where $n=$ `y.len() + 1)
    /// - Time complexity: $O(nm)$
    ///
    /// # References
    ///
    /// - [Huang, X. and Miller, W. 1991. A time-efficient linear-space local similarity algorithm. Adv. Appl. Math. 12, 3 (Sep. 1991), 337-357](https://doi.org/10.1016/0196-8858(91)90017-D)
    fn find_local_score_and_termini_huang(
        &self,
        x: Seq,
        y: Seq,
    ) -> (i32, usize, usize, usize, usize) {
        let (m, n) = (x.len() + 1, y.len() + 1);
        let mut cc: Vec<i32> = vec![0; n]; // match/mismatch           32 * n bits
        let mut dd: Vec<i32> = vec![MIN_SCORE; n]; // deletion          32 * n bits
        let mut origin_cc: Vec<[usize; 2]> = Vec::with_capacity(n); // usize * 2 * n bits
        let mut origin_dd: Vec<[usize; 2]> = vec![[0, 0]; n]; //       usize * 2 * n bits
        let mut e_origin: [usize; 2] = [0, 0]; // this value will not be used
        for j in 0..n {
            origin_cc.push([0, j]);
            // origin_dd.push([0, j]); values won't be used;
            // dd[0] = [MIN, MIN, ...], dd[0][j] garanteed to be less than c[j] + self.scoring.gap_open
            // origin_ii.push([0, j]); values won't be used
        }
        let mut e: i32; // I(i, j-1)
        let mut c: i32; // C(i, j-1)
        let mut c_origin;
        let mut max_c = MIN_SCORE;
        let mut c_start = [0usize, 0usize];
        let mut c_end = [m, n];
        let mut s: i32; // C(i-1, j-1)
        let mut s_origin: [usize; 2]; // origin of C(i-1, j-1)
        let mut p: u8;
        let mut q: u8;
        for i in 1..m {
            // s and cc[0] = 0; cc[0] always equals to 0
            s = 0;
            s_origin = [i - 1, 0];
            c = 0;
            c_origin = [i, 0];
            origin_cc[0] = [i, 0];
            e = MIN_SCORE; // I[i, 0] garanteed to be less than c + self.scoring.gap_open
                           // origin_dd[0] = [i, 0]; // ! value won't be used
                           // e_origin = [i, 0]; // ! value won't be used
            p = x[i - 1];
            for j in 1..n {
                q = y[j - 1];
                e = if e > c + self.scoring.gap_open {
                    // e_origin = e_origin;
                    e
                } else {
                    e_origin = c_origin;
                    c + self.scoring.gap_open
                } + self.scoring.gap_extend; // update e to I[i,j]
                dd[j] = if dd[j] > cc[j] + self.scoring.gap_open {
                    // origin_dd[j] = origin_dd[j];
                    dd[j]
                } else {
                    origin_dd[j] = origin_cc[j];
                    cc[j]
                } + self.scoring.gap_extend; // cc[j] = C[i-1, j]
                c = s + self.scoring.match_fn.score(p, q); // substitution score
                c_origin = s_origin;
                s = cc[j];
                s_origin = origin_cc[j];
                if dd[j] > c {
                    c = dd[j];
                    c_origin = origin_dd[j];
                }
                if e > c {
                    c = e;
                    c_origin = e_origin;
                }
                if c < 0 {
                    // the critical step in local alignment
                    c = 0;
                    c_origin = [i, j]
                }
                if c >= max_c {
                    max_c = c;
                    c_start = c_origin;
                    c_end = [i, j];
                }
                origin_cc[j] = c_origin;
                cc[j] = c;
            }
        }
        (max_c, c_start[0], c_start[1], c_end[0], c_end[1])
    }

    /// Find the maximum score in a local alignment between `x` and `y` and the coordinates
    /// of the alignment path termini, based on Shamir's approach
    ///
    /// - Space complexity: $O(n)$; specifically, about $64n$ bits, where $n=$ `y.len() + 1)
    /// - Time complexity: $O(nm)$ (slightly slower than the original `find_local_score_and_termini`)
    #[allow(dead_code)]
    fn find_local_score_and_termini_shamir(
        &self,
        x: Seq,
        y: Seq,
    ) -> (i32, usize, usize, usize, usize) {
        let (m, n) = (x.len() + 1, y.len() + 1);
        let mut cc: Vec<i32> = vec![0; n]; // match/mismatch           32 * n bits
        let mut dd: Vec<i32> = vec![MIN_SCORE; n]; // deletion          32 * n bits
        let mut e: i32; // I(i, j-1)
        let mut c: i32; // C(i, j-1)
        let mut max_c = MIN_SCORE;
        let mut xstart = 0usize;
        let mut ystart = 0usize;
        let mut xend = m;
        let mut yend = n;
        let mut s: i32; // C(i-1, j-1)
        let mut p: u8;
        let mut q: u8;
        for i in 1..m {
            s = 0;
            c = 0;
            e = MIN_SCORE;
            p = x[i - 1];
            for j in 1..n {
                q = y[j - 1];
                e = max(e, c + self.scoring.gap_open) + self.scoring.gap_extend;
                dd[j] = max(dd[j], cc[j] + self.scoring.gap_open) + self.scoring.gap_extend;
                c = max(s + self.scoring.match_fn.score(p, q), max(dd[j], max(e, 0)));
                s = cc[j];
                if c >= max_c {
                    max_c = c;
                    xend = i;
                    yend = j;
                }
                cc[j] = c;
            }
        }
        // run in reverse
        let x = &x[..xend];
        let y = &y[..yend];
        cc = vec![0; yend + 1];
        dd = vec![MIN_SCORE; yend + 1];
        max_c = MIN_SCORE;
        for i in 1..=xend {
            s = 0;
            c = 0;
            e = MIN_SCORE;
            p = x[xend - i];
            for j in 1..=yend {
                q = y[yend - j];
                e = max(e, c + self.scoring.gap_open) + self.scoring.gap_extend;
                dd[j] = max(dd[j], cc[j] + self.scoring.gap_open) + self.scoring.gap_extend;
                c = s + self.scoring.match_fn.score(p, q);
                s = cc[j];
                if dd[j] > c {
                    c = dd[j];
                }
                if e > c {
                    c = e;
                }
                if c < 0 {
                    c = 0;
                }
                if c >= max_c {
                    max_c = c;
                    xstart = i; // subtract from xend later
                    ystart = j;
                }
                cc[j] = c;
            }
        }
        xstart = xend - xstart;
        ystart = yend - ystart;
        (max_c, xstart, ystart, xend, yend)
    }

    /// Find the maximum score in a semiglobal alignment of `x` against `y` (x is global, y is local), // TODO: change: y agianst x; x local, y global
    /// and the coordinates of the alignment path termini
    ///
    /// In semiglobal mode, `xstart === 0` and `xend === m`, so only `ystart` and `yend` are computed
    /// and returned
    fn find_semiglobal_score_and_termini(&self, x: Seq, y: Seq) -> (i32, usize, usize) {
        let (m, n) = (x.len(), y.len());
        let mut cc: Vec<i32> = Vec::with_capacity(n + 1); //                32 * n bits
        let mut dd: Vec<i32> = vec![MIN_SCORE; n + 1]; //                   32 * n bits
        let mut x_origin_cc: Vec<usize> = vec![0; n + 1]; //                usize * n bits
        let mut x_origin_dd: Vec<usize> = vec![0; n + 1]; //                usize * n bits
        let mut e: i32; // I(i, j-1)
        let mut e_x_origin: usize;
        let mut c: i32; // C(i, j-1)
        let mut c_x_origin;
        let mut s: i32; // C(i-1, j-1)
        let mut s_x_origin: usize; // x_origin of C(i-1, j-1)
        let mut p: u8;
        let mut q: u8;
        let mut xstart = 0;
        let mut xend = m;
        let mut max_last_col = MIN_SCORE;
        let mut t = self.scoring.gap_open;
        cc.push(0);
        for _j in 1..=n {
            t += self.scoring.gap_extend;
            cc.push(t);
        }
        for i in 1..=m {
            s = 0; //cc[0];
                   // t += self.scoring.gap_extend;
            c = 0; // t;
            cc[0] = 0; // c;
            e = MIN_SCORE;
            s_x_origin = i - 1;
            c_x_origin = i;
            e_x_origin = i;
            // x_origin_cc[0] = i;
            p = x[i - 1];
            for j in 1..=n {
                q = y[j - 1];
                e = if e > c + self.scoring.gap_open {
                    // e_x_origin = e_x_origin;
                    e
                } else {
                    e_x_origin = c_x_origin;
                    c + self.scoring.gap_open
                } + self.scoring.gap_extend; // update e to I[i,j]
                dd[j] = if dd[j] > cc[j] + self.scoring.gap_open {
                    // x_origin_dd[j] = x_origin_dd[j];
                    dd[j]
                } else {
                    x_origin_dd[j] = x_origin_cc[j];
                    cc[j]
                } + self.scoring.gap_extend; // cc[j] = C[i-1, j]
                c = s + self.scoring.match_fn.score(p, q); // substitution score
                c_x_origin = s_x_origin;
                s = cc[j];
                s_x_origin = x_origin_cc[j];
                if dd[j] > c {
                    c = dd[j];
                    c_x_origin = x_origin_dd[j];
                }
                if e > c {
                    c = e;
                    c_x_origin = e_x_origin;
                }
                cc[j] = c;
                x_origin_cc[j] = c_x_origin;

                if j == n && c > max_last_col {
                    // c + self.scoring.xclip_suffix
                    max_last_col = c;
                    xstart = c_x_origin;
                    xend = i;
                }
            }
        }
        (max_last_col, xstart, xend)
    }

    // ! Different from the original custom aligner, this one does not allow a gap to align to another gap!
    /// if all 4 clip scores are set to zero, then in all these alignments, terminal gaps are not penalized
    /// ```ignore
    /// ---TTGGCC  AAATTGG--  --TTGG---  AATTGGCCC
    /// AAATTGG--  ---TTGGCC  AATTGGCCC  --TTGG---
    /// ```
    /// But a gap is not allowed to align to another gap:
    /// ```ignore
    /// --AATTG  ATGAT--
    /// ---ATTG  ATGAT---
    /// ```
    #[allow(dead_code)]
    fn find_dual_semiglobal_score_and_termini(
        &self,
        x: Seq,
        y: Seq,
    ) -> (i32, usize, usize, usize, usize) {
        let (m, n) = (x.len() + 1, y.len() + 1);
        let mut cc: Vec<i32> = vec![0i32; n]; //                  32 * n bits
        let mut dd: Vec<i32> = vec![MIN_SCORE; n]; //             32 * n bits
        let mut origin_cc: Vec<usize> = vec![0; n]; //            usize * n bits
        let mut clip_x_cc: Vec<bool> = vec![false; n]; //         8 * n bits
        let mut origin_dd: Vec<usize> = vec![0; n]; //            usize * n bits
        let mut clip_x_dd: Vec<bool> = vec![false; n]; //         8 * n bits
        let mut origin_ii: Vec<usize> = vec![0; n]; //            usize * n bits
        let mut clip_x_ii: Vec<bool> = vec![false; n]; //         8 * n bits

        // About clip_x: false: clip y (start = [0, j]); true: clip x (start = [i, 0]), or start = [0,0]
        let mut e: i32; // I(i, j-1)
        let mut c: i32; // C(i, j-1)
        let mut max_last_column_or_row = MIN_SCORE; // tracks the maximum of the last column
        let mut c_origin;
        let mut c_clip_x: bool;
        let mut xstart = 0usize;
        let mut ystart = 0usize;
        let mut xend = m;
        let mut yend = n;
        let mut s: i32; // C(i-1, j-1)
        let mut s_origin: usize; // origin of C(i-1, j-1)
        let mut s_clip_x: bool;
        for i in 1..m {
            s = 0;
            c = 0;
            // cc[0] = 0; unchanged
            e = MIN_SCORE;

            s_origin = i - 1; // origin_cc[0] (prev)
            s_clip_x = true; //  clip_x_cc[0] (prev)
            c_origin = i;
            c_clip_x = true;
            origin_cc[0] = i;
            clip_x_cc[0] = true;
            for j in 1..n {
                e = if e > c + self.scoring.gap_open {
                    origin_ii[j] = origin_ii[j - 1];
                    clip_x_ii[j] = clip_x_ii[j - 1];
                    e
                } else {
                    origin_ii[j] = c_origin;
                    clip_x_ii[j] = c_clip_x;
                    c + self.scoring.gap_open
                } + self.scoring.gap_extend; // update e to I[i,j]
                dd[j] = if dd[j] > cc[j] + self.scoring.gap_open {
                    // origin_dd[j] = origin_dd[j];
                    // clip_x_dd[j] = clip_x_dd[j];
                    dd[j]
                } else {
                    origin_dd[j] = origin_cc[j];
                    clip_x_dd[j] = clip_x_cc[j];
                    cc[j]
                } + self.scoring.gap_extend; // cc[j] = C[i-1, j]
                c = s + self.scoring.match_fn.score(x[i - 1], y[j - 1]); // substitution score
                c_origin = s_origin;
                c_clip_x = s_clip_x;
                s = cc[j];
                s_origin = origin_cc[j];
                s_clip_x = clip_x_cc[j];
                if dd[j] > c {
                    c = dd[j];
                    c_origin = origin_dd[j];
                    c_clip_x = clip_x_dd[j];
                }
                if e > c {
                    c = e;
                    c_origin = origin_ii[j];
                    c_clip_x = clip_x_ii[j];
                }
                cc[j] = c;
                origin_cc[j] = c_origin;
                clip_x_cc[j] = c_clip_x;

                if j == n && c > max_last_column_or_row {
                    max_last_column_or_row = c;
                    xend = i;
                    // yend = n unchanged;
                    if c_clip_x {
                        xstart = c_origin;
                        ystart = 0;
                    } else {
                        xstart = 0;
                        ystart = c_origin;
                    }
                }
            }
        }
        // last (m-th) row
        max_last_column_or_row = if max_last_column_or_row > cc[n] {
            max_last_column_or_row
        } else {
            if clip_x_cc[n] {
                xstart = origin_cc[n];
                ystart = 0;
            } else {
                xstart = 0;
                ystart = origin_cc[n];
            }
            cc[n]
        };
        for j in 0..n {
            if cc[j] > max_last_column_or_row {
                max_last_column_or_row = cc[j];
                xend = m;
                yend = j;
                if clip_x_cc[j] {
                    xstart = origin_cc[j];
                    ystart = 0;
                } else {
                    xstart = 0;
                    ystart = origin_cc[j];
                }
            }
        }
        (max_last_column_or_row, xstart, ystart, xend, yend)
    }

    /// Compute the (global) alignment operations between a single letter sequence `x` and a
    /// second sequence `y`. The second sequence can be empty, i.e. `b""`
    fn nw_onerow(&self, x: u8, y: Seq, n: usize, tx: i32) -> Vec<AlignmentOperation> {
        let score_by_indels_only =
            tx + self.scoring.gap_extend * (n as i32 + 1) + self.scoring.gap_open;
        let mut max = score_by_indels_only;
        let score_with_one_substitution_base =
            (n as i32 - 1) * self.scoring.gap_extend + self.scoring.gap_open; // plus substitution score and possibly one more gap_open
        let mut maxj_ = 0usize;
        // for j_ in 0..n // TODO: speed?
        for (j_, q) in y.iter().enumerate() {
            // index of sequence instead of matrix; y[j] instead of j[j-1] is the jth character
            let score = score_with_one_substitution_base
                + self.scoring.match_fn.score(x, *q)
                + if j_ == 0 || j_ == n - 1 {
                    0
                } else {
                    self.scoring.gap_open
                };
            if score > max {
                max = score;
                maxj_ = j_;
            }
        }
        if max == score_by_indels_only {
            let mut res = Vec::with_capacity(n + 1);
            res.push(AlignmentOperation::Delete);
            for _j in 0..n {
                res.push(AlignmentOperation::Insert)
            }
            res
        } else {
            let mut res = Vec::with_capacity(n);
            for _j in 0..maxj_ {
                res.push(AlignmentOperation::Insert)
            }
            if x == y[maxj_] {
                res.push(AlignmentOperation::Match);
            } else {
                res.push(AlignmentOperation::Mismatch);
            }
            for _j in 0..(n - maxj_ - 1) {
                res.push(AlignmentOperation::Insert)
            }
            res
        }
    }
}
