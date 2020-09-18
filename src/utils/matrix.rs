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

//! A naïve implementation of a fixed-shape matrix.
//!
//! Typically, you initilize an empty $m \times n$ matrix with `Matrix::with_capacity(m, n)` and then push
//! the values row by row with `.push(v)`. You can also initialize a matrix filled with a constant value by
//! `Matrix::fill(m, n, v)`. Then, you can get and set the value at a given index $(i, j)$ by `.get(i, j)`
//! and `.set(i, j, v)`, respectively. There are also unsafe (i.e. without bound checking) versions of the
//! two methods, namely `.get_unsafe(i, j)` and `.set_unsafe(i, j, v)`, which do their job faster (about 2x
//! speed).

use std::fmt;

#[derive(Debug, PartialEq, Eq)]
pub struct Matrix<T: Clone + Copy> {
    pub nrow: usize,
    pub ncol: usize,
    pub values: Vec<T>,
}

impl<T: Clone + Copy> Matrix<T> {
    pub fn with_capacity(nrow: usize, ncol: usize) -> Self {
        Self {
            nrow,
            ncol,
            values: Vec::<T>::with_capacity(nrow * ncol),
        }
    }
    pub fn fill(nrow: usize, ncol: usize, v: T) -> Self {
        Self {
            nrow,
            ncol,
            values: vec![v; nrow * ncol],
        }
    }
    pub fn get(&self, i: usize, j: usize) -> T {
        self.values[i * self.ncol + j]
    }
    pub fn set(&mut self, i: usize, j: usize, v: T) {
        self.values[i * self.ncol + j] = v;
    }
    pub unsafe fn get_unsafe(&self, i: usize, j: usize) -> T {
        *self.values.get_unchecked(i * self.ncol + j)
    }
    pub unsafe fn set_unsafe(&mut self, i: usize, j: usize, v: T) {
        *self.values.get_unchecked_mut(i * self.ncol + j) = v;
    }
    pub fn push(&mut self, v: T) {
        self.values.push(v)
    }
}

impl<T: Clone + Copy + fmt::Display> fmt::Display for Matrix<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.nrow {
            for j in 0..self.ncol {
                write!(f, "{:>3} ", self.get(i, j))?;
            }
            write!(f, "\n")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_matrix_fill() {
        let m = Matrix::fill(3, 2, 0i32);
        assert_eq!(
            m,
            Matrix {
                nrow: 3,
                ncol: 2,
                values: vec![0i32; 6]
            }
        )
    }

    #[test]
    fn test_matrix_get() {
        let v = Matrix::fill(3, 2, 0i32).get(2, 0);
        assert_eq!(v, 0i32)
    }

    #[test]
    fn test_matrix_set() {
        let mut m = Matrix::fill(3, 2, 0i32);
        m.set(2, 0, 23i32);
        assert_eq!(m.get(2, 0), 23i32)
    }
}
