#[macro_export]
macro_rules! complex {
    ($re:expr, $im:expr) => {{
        Complex { re: $re, im: $im }
    }};
}

use float_cmp::approx_eq;
use std::ops::{Add, AddAssign, Div, Mul, Neg};

#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}
impl Neg for Complex {
    type Output = Complex;
    fn neg(self) -> Self::Output {
        Complex {
            re: -self.re,
            im: -self.im,
        }
    }
}
impl Add for Complex {
    type Output = Complex;
    fn add(self, rhs: Self) -> Self::Output {
        complex!(self.re + rhs.re, self.im + rhs.im)
    }
}
impl AddAssign for Complex {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}
impl Mul for Complex {
    type Output = Complex;
    fn mul(self, rhs: Self) -> Self::Output {
        let re = (self.re * rhs.re) - (self.im * rhs.im);
        let im = (self.im * rhs.re) + (self.re * rhs.im);
        complex!(re, im)
    }
}
impl Mul<f64> for Complex {
    type Output = Complex;
    fn mul(self, rhs: f64) -> Self::Output {
        complex!(self.re * rhs, self.im * rhs)
    }
}
impl Mul<Complex> for f64 {
    type Output = Complex;
    fn mul(self, rhs: Complex) -> Self::Output {
        rhs * self
    }
}
impl Div for Complex {
    type Output = Complex;
    fn div(self, rhs: Self) -> Self::Output {
        let denominator_scalar = rhs.re.powi(2) + rhs.im.powi(2);
        let re = ((self.re * rhs.re) + (self.im * rhs.im)) / denominator_scalar;
        let im = ((self.im * rhs.re) - (self.re * rhs.im)) / denominator_scalar;
        complex!(re, im)
    }
}
impl Div<f64> for Complex {
    type Output = Complex;
    fn div(self, rhs: f64) -> Self::Output {
        self / complex!(rhs, 0.0)
    }
}
impl Complex {
    #[inline]
    pub fn mag(self) -> f64 {
        (self.re.powi(2) + self.im.powi(2)).sqrt()
    }
    #[inline]
    pub fn conj(self) -> Complex {
        Complex {
            re: self.re,
            im: -self.im,
        }
    }
    pub fn equals(self, other: Complex) -> bool {
        approx_eq!(f64, self.re, other.re) && approx_eq!(f64, self.im, other.im)
    }
}
