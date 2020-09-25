use crate::utils::Seq;
use std::cmp::min;
use std::mem;

pub fn levenshtein_space_efficient(alpha: Seq<'_>, beta: Seq<'_>) -> u32 {
    let (m, n) = (alpha.len(), beta.len());
    let mut dp_matrix: Vec<u32> = Vec::with_capacity(n + 1); // the dynamic programming matrix (only 1 column stored)
    let mut s_diag: u32; // dp_matrix[i - 1][j - 1]
    let mut s_left: u32; // dp_matrix[i][j - 1]
    let mut a: u8; // alpha[i - 1]
    let mut b: u8; // beta[j - 1]

    // 0th row
    for j in 0..=(n as u32) {
        dp_matrix.push(j);
    }
    // rows 1 to m
    for i in 1..=m {
        s_diag = (i - 1) as u32;
        s_left = i as u32;
        a = alpha[i - 1];
        for j in 1..=n {
            b = beta[j - 1];
            s_left = min(
                s_diag + if a == b { 0 } else { 1 },
                min(s_left + 1, dp_matrix[j] + 1),
            );
            s_diag = dp_matrix[j];
            dp_matrix[j] = s_left;
        }
    }

    dp_matrix[n]
}

pub fn levenshtein_unsafe(alpha: Seq<'_>, beta: Seq<'_>) -> u32 {
    let (m, n) = (alpha.len(), beta.len());
    let mut dp_matrix: Vec<u32> = Vec::with_capacity(n + 1); // the dynamic programming matrix (only 1 column stored)
    let mut s_diag: u32; // dp_matrix[i - 1][j - 1]
    let mut s_left: u32; // dp_matrix[i][j - 1]
    let mut a: u8; // alpha[i - 1]
    let mut b: u8; // beta[j - 1]

    // 0th row
    for j in 0..=(n as u32) {
        dp_matrix.push(j);
    }
    // rows 1 to m
    for i in 1..=m {
        s_diag = (i - 1) as u32;
        s_left = i as u32;
        unsafe {
            a = *alpha.get_unchecked(i - 1);
        }
        for j in 1..=n {
            unsafe {
                b = *beta.get_unchecked(j - 1);
                s_left = min(
                    s_diag + if a == b { 0 } else { 1 },
                    min(s_left + 1, dp_matrix.get_unchecked(j) + 1),
                );
                s_diag = *dp_matrix.get_unchecked(j);
                *dp_matrix.get_unchecked_mut(j) = s_left;
            }
        }
    }

    dp_matrix[n]
}
