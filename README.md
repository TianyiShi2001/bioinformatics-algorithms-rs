# Learn Bioinformatics by Writing Entirely Too Many Algorithms in Rust

Bioinformatics Algorithms implemented in Rust **for bioinformatics learners**, with a [book](https://tianyishi2001.github.io/too-many-bioinformatics-algorithms).

## TL;DR

### Should I Read on?

This repository is for those who want to seriously learn bioinformatics algorithms without learning C/C++

### What Can I Learn?

- Bioinformatics algorithms, from those "ancient" ones found only in textbooks to the state-of-the-art ones currently being used and discussed
- You'll understand those algorithms, and be able to implement them by yourself.
- **This is a very young project initiated on 2020-09-12 and maintained by an undergraduate student, so currently don't expect too much.**

### Why Rust/What is Rust?

- Simple programming languages such as Python, Ruby and Perl are popular among biologists. That's why BioPython, BioRuby and BioPerl are much popular than SeqAn (C++).
- Python, Ruby and Perl are too slow for core algorithms in bioinformatics, and these algorithms are traditionally implemented in C/C++.
- C/C++ is difficult for many biologists. Of those interested in the internals of the bioinformatics tools they use, few take a step further to learn C/C++.
- A young language called [**Rust**](https://www.rust-lang.org) turns out to be both simple and fast (as fast as C). Indeed, some algorithms implemented in [**rust-bio**](https://github.com/rust-bio/rust-bio) outperforms their counterparts in SeqAn.
- Thus I believe Rust is the perfect language for bioinformatics. People who want to learn bioinformatics algorithms seriously should learn Rust.
- This repository helps you to **learn** bioinformatics with Rust.
- It is not just a bunch of code to be downloaded and run. It helps you to **learn** by guiding you through the steps to reach the final solution. Any special implementation details are thoroughly explained, and the only prerequisite is Chapters 1-14 of [The Book](https://doc.rust-lang.org/book/).

# Ready to go?

Read the book [here](https://tianyishi2001.github.io/too-many-bioinformatics-algorithms).

# Contribution

Contributions of any form are welcome!

Currently this project is at is very early stage of development, so I don't mind any breaking changes.

In addition, English is not my first language, so I would be grateful if you can help rephrasing words that can be improved!

# Roadmap

✅: Algorithm Implemented
✅✅: Comments Added
✅✅✅: Documentation Finished

## Alignment

### Pairwise Alignment

| Implemented |       Author        | Space complexity | Global | Local | Semiglobal | Affine Gap | Log Gap |
| :---------: | :-----------------: | :--------------: | :----: | :---: | :--------: | :--------: | :-----: |
|      ✅      |  Needleman-Wunsch   |   **_O(nm)_**    |   ✅    |       |            |            |         |
|      ✅      |   Smith-Waterman    |   **_O(nm)_**    |        |   ✅   |            |            |         |
|             |        Gotoh        |   **_O(nm)_**    |   ✅    |       |            |     ✅      |         |
|             |      Hirshberg      |    **_O(n)_**    |   ✅    |       |            |            |         |
|             | Myers-Miller (1988) |    **_O(n)_**    |   ✅    |       |            |     ✅      |         |
|             |    Huang (1991)     |    **_O(n)_**    |        |   ✅   |            |     ✅      |         |

### Multiple Alignment

MUSCLE

ClustalW