use crate::complex;
use crate::{Complex, QState, Qubit};
use std::f64::consts::{FRAC_1_SQRT_2, SQRT_2};

// Inverse Root 2
pub const IR2: f64 = FRAC_1_SQRT_2;

pub const C_ZERO: Complex = complex!(0.0, 0.0);
pub const C_ONE: Complex = complex!(1.0, 0.0);
mod gate {
    use super::IR2;
    use crate::Complex;
    use crate::{complex, matrix::QMatrix};
    use lazy_static::lazy_static;
}
