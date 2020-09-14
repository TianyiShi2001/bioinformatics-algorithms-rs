#![feature(test)]

extern crate test;

use test::Bencher;

#[bench]
fn bench_vec_of_vec(b: &mut Bencher) {
    let (m, n) = (1000, 10000);
    let mut matrix = vec![vec![0; n]; m];
    b.iter(|| {
        for i in 0..m {
            for j in 0..n {
                matrix[i][j] = i * j;
            }
        }
    });
}

#[bench]
fn bench_vec(b: &mut Bencher) {
    let (m, n) = (1000, 10000);
    let mut matrix = vec![0; n * m];
    b.iter(|| {
        for i in 0..m {
            for j in 0..n {
                matrix[i * n + j] = i * j;
            }
        }
    });
}
