use crate::complex::Complex;
use crate::constants::*;
use crate::{complex, matrix, QuantumVec};
use lazy_static::lazy_static;
use std::ops::Mul;
#[derive(Debug, Clone, Default, PartialEq)]
/// A square matrix representing a quantum gate.
/// Must be of dimensions 2^n x 2^n for some n.
pub struct QMatrix {
    dims: usize,
    data: Vec<Complex>,
}
impl QuantumVec for QMatrix {
    fn data_slice(&self) -> &[Complex] {
        &self.data
    }
}

impl Mul<QMatrix> for Complex {
    type Output = QMatrix;
    fn mul(self, rhs: QMatrix) -> Self::Output {
        QMatrix {
            dims: rhs.dims,
            data: rhs.data.into_iter().map(|z| self * z).collect(),
        }
    }
}
impl Mul<QMatrix> for f64 {
    type Output = QMatrix;
    fn mul(self, rhs: QMatrix) -> Self::Output {
        complex!(self, 0.0) * rhs
    }
}
impl QMatrix {
    pub fn from_data(data: Vec<Complex>) -> QMatrix {
        assert!(is_square_number(data.len()));
        let dims = (data.len() as f32).sqrt() as usize;
        assert!(dims.is_power_of_two());
        QMatrix { dims, data }
    }
    pub fn from_vecs(data: Vec<Vec<Complex>>) -> QMatrix {
        let rows = data.len();
        // QMatrix must be square
        assert!(data.iter().all(|row| row.len() == rows));
        QMatrix {
            dims: rows,
            data: data.into_iter().flatten().collect(),
        }
    }
    pub fn mul(&self, lhs: &QMatrix) -> QMatrix {
        assert_eq!(self.dims, lhs.dims);

        let dims = self.dims;
        let mut matrix_data = vec![Complex::default(); dims.pow(2)];
        for col in 0..dims {
            for row in 0..dims {
                for ele in 0..dims {
                    matrix_data[col * dims + row] +=
                        lhs.data[ele * dims + row] * self.data[col * dims + ele];
                }
            }
        }

        QMatrix {
            dims,
            data: matrix_data,
        }
    }
    pub fn kronecker(&self, lhs: &QMatrix) -> QMatrix {
        let matrix_dims = self.dims.pow(2);
        let sub_dims = self.dims;
        let mut matrix_data = vec![Complex::default(); matrix_dims.pow(2)];
        for lhs_col in 0..sub_dims {
            for lhs_row in 0..sub_dims {
                for col in 0..sub_dims {
                    for row in 0..sub_dims {
                        let matrix_data_idx = lhs_col * sub_dims.pow(3)
                            + lhs_row * sub_dims.pow(2)
                            + col * sub_dims
                            + row;
                        let lhs_idx = sub_dims * lhs_col + col;
                        let self_idx = sub_dims * lhs_row + row;

                        matrix_data[matrix_data_idx] = lhs.data[lhs_idx] * self.data[self_idx];
                        println!(
                            "matrix_data index: {}, lhs index: {}, self index: {}",
                            matrix_data_idx, lhs_idx, self_idx,
                        );
                    }
                }
            }
        }
        QMatrix {
            dims: matrix_dims,
            data: matrix_data,
        }
    }
    // self[0 * cols + 0] * lhs[0 * cols + 0],
    // self[0 * cols + 0] * lhs[0 * cols + 1]
}
fn is_square_number(num: usize) -> bool {
    let sqrt = (num as f32).sqrt() as usize;
    sqrt.pow(2) == num
}
lazy_static! {
    pub static ref X_GATE: QMatrix = QMatrix::from_data(vec![C_ZERO, C_ONE, C_ONE, C_ZERO]);
    pub static ref Y_GATE: QMatrix = QMatrix::from_data(vec![
        C_ZERO,
        complex!(0.0, 1.0),
        complex!(0.0, -1.0),
        C_ZERO
    ]);
    pub static ref Z_GATE: QMatrix = QMatrix::from_data(vec![C_ONE, C_ZERO, C_ZERO, -C_ONE]);
    pub static ref H_GATE: QMatrix = QMatrix::from_data(vec![
        complex!(IR2, 0.0),
        complex!(IR2, 0.0),
        complex!(IR2, 0.0),
        complex!(-IR2, 0.0),
    ]);
    pub static ref S_GATE: QMatrix =
        QMatrix::from_data(vec![C_ONE, C_ZERO, C_ZERO, complex!(0.0, 1.0)]);
    pub static ref T_GATE: QMatrix =
        QMatrix::from_data(vec![C_ONE, C_ZERO, C_ZERO, IR2 * complex!(1.0, 1.0)]);
    pub static ref CNOT_GATE: QMatrix = QMatrix::from_data(vec![
        C_ONE, C_ZERO, C_ZERO, C_ZERO, C_ZERO, C_ONE, C_ZERO, C_ZERO, C_ZERO, C_ZERO, C_ZERO,
        C_ONE, C_ZERO, C_ZERO, C_ONE, C_ZERO
    ]);
    pub static ref CZ: QMatrix = QMatrix::from_data(vec![
        C_ONE, C_ZERO, C_ZERO, C_ZERO, C_ZERO, C_ONE, C_ZERO, C_ZERO, C_ZERO, C_ZERO, C_ONE,
        C_ZERO, C_ZERO, C_ZERO, C_ZERO, -C_ONE
    ]);
    pub static ref SWAP_GATE: QMatrix = QMatrix::from_data(vec![
        C_ONE, C_ZERO, C_ZERO, C_ZERO, C_ZERO, C_ZERO, C_ONE, C_ZERO, C_ZERO, C_ONE, C_ZERO,
        C_ZERO, C_ZERO, C_ZERO, C_ZERO, C_ONE
    ]);
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic]
    fn creating_qmatrix_from_data_with_non_square_dims_panics() {
        let data = vec![Complex::default(); 9999];
        let _ = QMatrix::from_data(data);
    }

    #[test]
    fn muliply_two_qmatrices_gives_the_correct_result() {
        let matrix_1 = QMatrix {
            dims: 2,
            data: vec![
                complex!(1.0, 0.0),
                complex!(2.0, 0.0),
                complex!(3.0, 0.0),
                complex!(4.0, 0.0),
            ],
        };
        let matrix_2 = QMatrix {
            dims: 2,
            data: vec![
                complex!(1.0, 0.0),
                complex!(1.0, 0.0),
                complex!(1.0, 1.0),
                complex!(-1.0, -1.0),
            ],
        };
        let result = matrix_1.mul(&matrix_2);
        let expected_result = QMatrix {
            dims: 2,
            data: vec![
                complex!(3.0, 2.0),
                complex!(-1.0, -2.0),
                complex!(7.0, 4.0),
                complex!(-1.0, -4.0),
            ],
        };
        dbg!(&result);
        dbg!(&expected_result);
        assert_eq!(result, expected_result);
    }
    #[test]
    fn kronecker_product_of_two_matrices_yields_correct_output() {
        let matrix_1 = QMatrix {
            dims: 2,
            data: vec![
                complex!(1.0, 0.0),
                complex!(2.0, 0.0),
                complex!(3.0, 0.0),
                complex!(4.0, 0.0),
            ],
        };
        let matrix_2 = QMatrix {
            dims: 2,
            data: vec![
                complex!(0.0, 1.0),
                complex!(1.0, 1.0),
                complex!(0.0, 1.0),
                complex!(-1.0, 1.0),
            ],
        };
        let result = matrix_1.kronecker(&matrix_2);
        let expected_result = QMatrix {
            dims: 4,
            data: vec![
                complex!(0.0, 1.0),
                complex!(0.0, 2.0),
                complex!(1.0, 1.0),
                complex!(2.0, 2.0),
                complex!(0.0, 3.0),
                complex!(0.0, 4.0),
                complex!(3.0, 3.0),
                complex!(4.0, 4.0),
                complex!(0.0, 1.0),
                complex!(0.0, 2.0),
                complex!(-1.0, 1.0),
                complex!(-2.0, 2.0),
                complex!(0.0, 3.0),
                complex!(0.0, 4.0),
                complex!(-3.0, 3.0),
                complex!(-4.0, 4.0),
            ],
        };
        assert!(result.equals(&expected_result));
    }
}
