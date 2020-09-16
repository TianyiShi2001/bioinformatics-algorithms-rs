use bio_edu::alignment::pairwise::*;
use bio_edu::utils::matrix::Matrix;

fn main() {
    let x = b"ATGCTGCCGC";
    let y = b"ATGCTTGCATGC";
    // let match_fn = |a, b| {
    //     if a == b {
    //         1i32
    //     } else {
    //         -1i32
    //     }
    // };
    // let aligner = needleman_wunsch_2::Aligner::new(match_fn, -1);
    // //let res = aligner.global();
    // //println!("{:?}", res);
    // let (mut S, mut T) = aligner.init_matrices(x.len(), y.len());
    // aligner.fill_matrices(&mut S, &mut T, x, y);
    // println!("{}", T);

    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    let aligner = myers_miller::Aligner::new(-5, -1, score);
    let alignment = aligner.global(x, y);
    println!("{:?}", alignment);
    println!("{}", alignment.pretty(x, y));
}
