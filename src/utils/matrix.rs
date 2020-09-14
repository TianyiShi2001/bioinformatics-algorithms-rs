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

use std::default::Default;
use std::fmt;

#[derive(Debug, PartialEq, Eq)]
pub struct Matrix<T: Clone + Copy> {
    pub nrow: usize,
    pub ncol: usize,
    pub value: Vec<T>,
}

impl<T: Clone + Copy> Default for Matrix<T> {
    fn default() -> Self {
        Matrix {
            ncol: 0,
            nrow: 0,
            value: Vec::new(),
        }
    }
}

impl<T: Clone + Copy> Matrix<T> {
    pub fn new() -> Self {
        Matrix::default()
    }
    pub fn with_capacity(nrow: usize, ncol: usize) -> Self {
        Self {
            nrow,
            ncol,
            value: Vec::<T>::with_capacity(nrow * ncol),
        }
    }
    pub fn fill(nrow: usize, ncol: usize, v: T) -> Self {
        Self {
            nrow,
            ncol,
            value: vec![v; nrow * ncol],
        }
    }
    pub fn get(&self, i: usize, j: usize) -> T {
        self.value[i * self.ncol + j]
    }
    pub fn get_mut(&mut self, i: usize, j: usize) -> &mut T {
        &mut self.value[i * self.ncol + j]
    }
    pub fn set(&mut self, i: usize, j: usize, v: T) {
        self.value[i * self.ncol + j] = v;
    }
    pub fn push(&mut self, v: T) {
        self.value.push(v)
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
    fn test_matrix_new() {
        let m = Matrix::<i32>::new();
        assert_eq!(m.ncol, 0);
        assert_eq!(m.nrow, 0);
        assert_eq!(m.value, Vec::<i32>::new());
    }

    #[test]
    fn test_matrix_fill() {
        let m = Matrix::fill(3, 2, 0i32);
        assert_eq!(
            m,
            Matrix {
                nrow: 3,
                ncol: 2,
                value: vec![0i32; 6]
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
