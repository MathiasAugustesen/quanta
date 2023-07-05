use crate::complex;
use crate::matrix::QMatrix;
use crate::Complex;
use lazy_static::lazy_static;
use std::f64::consts::FRAC_1_SQRT_2;

// Inverse Root 2
pub const IR2: f64 = FRAC_1_SQRT_2;

pub const ZERO: Complex = complex!(0.0, 0.0);
pub const ONE: Complex = complex!(1.0, 0.0);
pub const I: Complex = complex!(0.0, 1.0);
pub const C_IR2: Complex = complex!(IR2, 0.0);

lazy_static! {
    pub static ref X_GATE: QMatrix = QMatrix::from_data(vec![ZERO, ONE, ONE, ZERO]);
    pub static ref Y_GATE: QMatrix = QMatrix::from_data(vec![ZERO, -I, I, ZERO]);
    pub static ref Z_GATE: QMatrix = QMatrix::from_data(vec![ONE, ZERO, ZERO, -ONE]);
    pub static ref H_GATE: QMatrix = QMatrix::from_data(vec![C_IR2, C_IR2, C_IR2, -C_IR2,]);
    pub static ref S_GATE: QMatrix = QMatrix::from_data(vec![ONE, ZERO, ZERO, I]);
    pub static ref T_GATE: QMatrix =
        QMatrix::from_data(vec![ONE, ZERO, ZERO, IR2 * complex!(1.0, 1.0)]);
    #[rustfmt::skip]
    pub static ref CNOT_GATE: QMatrix = QMatrix::from_data(vec![
        ONE, ZERO, ZERO, ZERO,
        ZERO, ONE, ZERO, ZERO,
        ZERO, ZERO, ZERO, ONE,
        ZERO, ZERO, ONE, ZERO
    ]);

    pub static ref CZ: QMatrix = QMatrix::from_data(vec![
        ONE, ZERO, ZERO, ZERO,
        ZERO, ONE, ZERO, ZERO,
        ZERO, ZERO, ONE, ZERO,
        ZERO, ZERO, ZERO, -ONE
    ]);
    pub static ref SWAP_GATE: QMatrix = QMatrix::from_data(vec![
        ONE, ZERO, ZERO, ZERO,
        ZERO, ZERO, ONE, ZERO,
        ZERO, ONE, ZERO, ZERO,
        ZERO, ZERO, ZERO, ONE
    ]);
}
