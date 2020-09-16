use bio_edu::alignment::pairwise::*;
use bio_edu::utils::matrix::Matrix;

fn main() {
    let x = b"CCCATGCTGCCGCAAA";
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

    let aligner = smith_waterman_1::Aligner::new(score, -2);
    let alignment = aligner.local(x, y); 
    println!("{:?}", alignment);
    let aligner2 = smith_waterman_1::Aligner::from_scores(1, -1, -2);
    let alignment = aligner2.local(x, y);
    println!("{:?}", alignment);
}
