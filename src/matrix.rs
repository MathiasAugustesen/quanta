use crate::complex;
use crate::complex::Complex;
#[derive(Debug, Clone, Default, PartialEq)]
pub struct QMatrix {
    dims: usize,
    data: Vec<Complex>,
}
impl QMatrix {
    pub fn from_data(data: Vec<Complex>) -> QMatrix {
        let dims = data.len();
        assert!(is_square_number(dims));
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
    pub fn combine(&self, other: QMatrix) -> QMatrix {
        todo!()
    }
}
fn is_square_number(num: usize) -> bool {
    let sqrt = (num as f32).sqrt() as usize;
    sqrt.pow(2) == num
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
        let result = matrix_1.combine(matrix_2);
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
}
