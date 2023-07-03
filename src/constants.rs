use crate::complex;
use crate::matrix::QMatrix;
use crate::Complex;
use lazy_static::lazy_static;
use std::f64::consts::FRAC_1_SQRT_2;

// Inverse Root 2
pub const IR2: f64 = FRAC_1_SQRT_2;

pub const C_ZERO: Complex = complex!(0.0, 0.0);
pub const C_ONE: Complex = complex!(1.0, 0.0);

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
