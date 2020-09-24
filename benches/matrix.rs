#![feature(test)]

extern crate test;
use na::Matrix;
use na::{Dynamic, VecStorage};
use nalgebra as na;

type DMatrixi32 = Matrix<u8, Dynamic, Dynamic, VecStorage<u8, Dynamic, Dynamic>>;

use test::Bencher;

#[bench]
fn bench_vec_of_vec(b: &mut Bencher) {
    let (m, n) = (10000, 10000);
    let mut matrix = vec![vec![0u8; n]; m];
    b.iter(|| {
        for i in 0..m {
            for j in 0..n {
                matrix[i][j] = 1u8;
            }
        }
    });
}

#[bench]
fn bench_vec(b: &mut Bencher) {
    let (m, n) = (10000, 10000);
    let mut matrix = vec![0u8; n * m];
    b.iter(|| {
        for i in 0..m {
            for j in 0..n {
                matrix[i * n + j] = 1u8;
            }
        }
    });
}

#[bench]
fn bench_vec_of_vec_unsafe(b: &mut Bencher) {
    let (m, n) = (10000, 10000);
    let mut matrix = vec![vec![0u8; n]; m];
    b.iter(|| {
        for i in 0..m {
            for j in 0..n {
                unsafe {
                    *matrix.get_unchecked_mut(i).get_unchecked_mut(j) = 1u8;
                }
            }
        }
    });
}

#[bench]
fn bench_vec_unsafe(b: &mut Bencher) {
    let (m, n) = (10000, 10000);
    let mut matrix = vec![0u8; n * m];
    b.iter(|| {
        for i in 0..m {
            for j in 0..n {
                unsafe { *matrix.get_unchecked_mut(i * n + j) = 1u8 };
            }
        }
    });
}

#[bench]
fn bench_nalgebra(b: &mut Bencher) {
    let (m, n) = (10000, 10000);
    let mut matrix = DMatrixi32::from_vec(m, n, vec![0u8; n * m]);
    b.iter(|| {
        for i in 0..m {
            for j in 0..n {
                matrix[(i, j)] = 1u8;
            }
        }
    });
}

// test bench_vec               ... bench:  84,694,933 ns/iter (+/- 7,412,836)
// test bench_vec_of_vec        ... bench:  87,083,636 ns/iter (+/- 1,171,842)
// test bench_vec_unsafe        ... bench:  41,440,947 ns/iter (+/- 752,463)
// test bench_vec_of_vec_unsafe ... bench:  44,532,595 ns/iter (+/- 629,209)
