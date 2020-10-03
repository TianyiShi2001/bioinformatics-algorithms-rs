// Copyright 2020 Tianyi Shi
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

#![allow(non_snake_case)]
use crate::alignment::pairwise::MIN_SCORE;
use crate::alignment::pairwise::{traceback::*, MatchFunc, Scoring};
use crate::alignment::pairwise::{Alignment, AlignmentMode, AlignmentOperation, GlobalAlign};
use crate::utils::Seq;

pub struct Aligner<F: MatchFunc> {
    scoring: Scoring<F>,
}

impl<F: MatchFunc> GlobalAlign for Aligner<F> {
    /// Fast global alignment with $O(nm)$ space. Used when `y.len()` is small.
    ///
    /// # Implementation Details
    ///
    /// ## Traceback Matrix `T: Vec<u8>`
    ///
    /// ```ignore
    /// 0b00   0b01    0b10    0b11
    /// start  insert  delete  match_or_subst
    /// ```
    ///
    /// ```ignore
    /// 0b00101101
    ///     / |  \
    ///    S  D   I
    /// ```
    fn global<'a>(&self, x: Seq<'a>, y: Seq<'a>) -> Alignment<'a> {
        let (m, n) = (x.len() + 1, y.len() + 1);
        let mut S = vec![0; m]; //                                            ! 32 * m bits
        let mut D = vec![MIN_SCORE; m]; //                                    ! 32 * m bits
        let mut T: Vec<TracebackCell> = vec![TracebackCell::new(); n * m]; // ! 8 * n * m bits
        let mut s: i32; // S[i - 1][j - 1] or S[i - 1][j]
        let mut c: i32; // S[i][j - 1] or S[i][j]
        let mut e: i32; // I[i - 1] or I [i]
        let mut idx: usize = 0; // i * n + j
        let mut x_i: &u8; // x[i] (x[i - 1]) // TODO: is not using p/q faster?
        let mut y_j: &u8; // y[j] (y[j - 1])
        let mut S_j: &mut i32; // S[j]
        let mut D_j: &mut i32; // D[j]
        let mut score_1: i32; // used when determining D[j] and I[j]
        let mut score_2: i32; // used when determining D[j] and I[j]

        // SAFETY: unchecked indexing is used here. x, y, S, D, T all have a fixed size related to n and/or m;
        //         it should have been implied by the for loops that all indexing operations are in-bound but
        //         the compiler wasn't smart enough to notice this as of October 2020.
        unsafe {
            // T[0] = TracebackCell::new() // origin at T[0 * n + 0]
            let mut t = self.scoring.gap_open;
            for j in 1..n {
                t += self.scoring.gap_extend;
                // I[0][j] = t will not be read
                *S.get_unchecked_mut(j) = t;
                let mut tb = TracebackCell::new();
                tb.set_s_bits(TB_LEFT);
                idx += 1;
                *T.get_unchecked_mut(idx) = tb; // T[0 * n + j]
            }

            t = self.scoring.gap_open;
            for i in 1..m {
                s = *S.get_unchecked(0);
                t += self.scoring.gap_extend;
                c = t;
                *S.get_unchecked_mut(0) = c;
                e = MIN_SCORE;
                // D[0] = t will not be read
                let mut tb = TracebackCell::new();
                tb.set_s_bits(TB_UP);
                idx += 1;
                *T.get_unchecked_mut(idx) = tb; // T[j * m + 0]

                x_i = x.get_unchecked(i - 1);
                for j in 1..n {
                    y_j = y.get_unchecked(j - 1);
                    S_j = S.get_unchecked_mut(j);
                    D_j = D.get_unchecked_mut(j);
                    let mut tb = TracebackCell::new();

                    score_1 = e + self.scoring.gap_extend;
                    score_2 = c + self.scoring.gap_open + self.scoring.gap_extend;
                    e = if score_1 > score_2 {
                        tb.set_i_bits(TB_LEFT);
                        score_1
                    } else {
                        tb.set_i_bits(T.get_unchecked(idx).get_s_bits()); // T[i-1][j]
                        score_2
                    };

                    idx += 1; // ! update idx

                    score_1 = *D_j + self.scoring.gap_extend;
                    score_2 = *S_j + self.scoring.gap_open + self.scoring.gap_extend;
                    *D_j = if score_1 > score_2 {
                        tb.set_d_bits(TB_UP);
                        score_1
                    } else {
                        tb.set_d_bits(T.get_unchecked(idx - n).get_s_bits()); //T[i][j-1]
                        score_2
                    };

                    c = s + self.scoring.match_fn.score(*x_i, *y_j);
                    tb.set_s_bits(TB_DIAG); // no need to be exact at this stage

                    if e > c {
                        c = e;
                        tb.set_s_bits(TB_LEFT);
                    }

                    if *D_j > c {
                        c = *D_j;
                        tb.set_s_bits(TB_UP);
                    }

                    s = *S_j;
                    *S_j = c;
                    *T.get_unchecked_mut(idx) = tb;
                }
            }

            let (mut operations, i, j) = Self::traceback(&T, x, y, m - 1, n - 1, n);
            operations.resize(operations.len() + i, AlignmentOperation::Delete); // reaching at (i, 0)
            operations.resize(operations.len() + j, AlignmentOperation::Insert); // reaching at (0, j)

            operations.reverse();
            Alignment {
                x,
                y,
                score: *S.get_unchecked(n - 1),
                xstart: 0,
                ystart: 0,
                xend: m - 1,
                yend: n - 1,
                operations,
                mode: AlignmentMode::Global,
            }
        }
    }
}

impl<F: MatchFunc> Aligner<F> {
    pub fn new(gap_open: i32, gap_extend: i32, match_fn: F) -> Self {
        Aligner {
            scoring: Scoring::new(gap_open, gap_extend, match_fn),
        }
    }

    pub unsafe fn traceback(
        T: &[TracebackCell],
        x: Seq,
        y: Seq,
        mut i: usize,
        mut j: usize,
        n: usize,
    ) -> (Vec<AlignmentOperation>, usize, usize) {
        let mut operations = Vec::with_capacity(n);
        let mut next_layer = T.get_unchecked(i * n + j).get_s_bits(); // start from the last tb cell
        loop {
            match next_layer {
                TB_ORIGIN => break,
                TB_LEFT => {
                    operations.push(AlignmentOperation::Delete);
                    next_layer = T.get_unchecked(i * n + j).get_i_bits();
                    j -= 1;
                }
                TB_UP => {
                    operations.push(AlignmentOperation::Insert);
                    next_layer = T.get_unchecked(i * n + j).get_d_bits();
                    i -= 1;
                }
                TB_DIAG => {
                    i -= 1;
                    j -= 1;
                    next_layer = T.get_unchecked(i * n + j).get_s_bits(); // T[i - 1][j - 1]
                    operations.push(if *y.get_unchecked(j) == *x.get_unchecked(i) {
                        AlignmentOperation::Match
                    } else {
                        AlignmentOperation::Mismatch
                    });
                }
                _ => unreachable!(),
            }
        }
        (operations, i, j)
    }
}
