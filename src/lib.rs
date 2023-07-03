pub mod complex;
pub mod constants;
pub mod matrix;
use complex::Complex;
use float_cmp::approx_eq;
use matrix::QMatrix;
trait QuantumVec {
    fn data_slice(&self) -> &[Complex];
    fn equals(&self, other: &impl QuantumVec) -> bool {
        self.data_slice().len() == self.data_slice().len()
            && self
                .data_slice()
                .iter()
                // Check that all elements are equal
                .zip(other.data_slice().iter())
                .all(|(a, b)| a.equals(*b))
    }
}
#[derive(Debug, Clone, Copy)]
pub struct Qubit {
    alpha: Complex,
    beta: Complex,
}
impl Qubit {
    pub fn new(alpha: Complex, beta: Complex) -> Qubit {
        let qubit = Qubit { alpha, beta };
        assert!(qubit.is_normalized());
        qubit
    }
    pub fn is_normalized(self) -> bool {
        approx_eq!(f64, self.alpha.mag().powi(2) + self.beta.mag().powi(2), 1.0)
    }
}
impl From<Qubit> for QState {
    fn from(value: Qubit) -> Self {
        QState {
            state: vec![value.alpha, value.beta],
        }
    }
}
#[derive(Debug, Clone)]
pub struct QState {
    state: Vec<Complex>,
}
impl QuantumVec for QState {
    fn data_slice(&self) -> &[Complex] {
        &self.state
    }
}
impl QState {
    pub fn from_qubits(qubits: &[Qubit]) -> Self {
        assert!(!qubits.is_empty());
        qubits
            .iter()
            .rev()
            .map(|&q| q.into())
            .reduce(|acc: QState, e| acc.tensor_product(e))
            .unwrap()
    }
    /// Calculates the tensor product between two quantum states as
    /// other ⊗ self.
    pub fn tensor_product(self, lhs: QState) -> Self {
        QState {
            state: lhs
                .state
                .iter()
                .flat_map(|&x| self.state.iter().map(move |&y| x * y))
                .collect(),
        }
    }
    pub fn apply(self, gate: &QMatrix) -> Self {
        todo!()
    }
}
#[cfg(test)]
mod tests {

    use super::*;
    use constants::*;

    #[test]
    fn tensor_product_of_two_qubit_states_yields_correct_output() {
        let qubit_1 = Qubit::new(complex!(0.5, 0.5), complex!(0.5, 0.5));
        let qubit_2 = Qubit::new(complex!(IR2, IR2), complex!(0.0, 0.0));

        let qstate_1: QState = qubit_1.into();
        let qstate_2: QState = qubit_2.into();
        let result = qstate_1.tensor_product(qstate_2);
        let expected_result = QState {
            state: vec![
                complex!(0.5, 0.5) * complex!(IR2, IR2),
                complex!(0.5, 0.5) * complex!(IR2, IR2),
                Complex::default(),
                Complex::default(),
            ],
        };
        assert!(expected_result.equals(&result));
    }
    #[test]
    fn tensor_product_of_random_states_gives_random_state() {
        let qstate_1: QState = Qubit::new(complex!(IR2, 0.0), complex!(IR2, 0.0)).into();
        let qstate_2: QState = Qubit::new(complex!(IR2, 0.0), complex!(IR2, 0.0)).into();

        let result = qstate_1.tensor_product(qstate_2);
        let expected_result = QState {
            state: vec![
                complex!(0.5, 0.0),
                complex!(0.5, 0.0),
                complex!(0.5, 0.0),
                complex!(0.5, 0.0),
            ],
        };
        dbg!(&expected_result);
        dbg!(&result);
        assert!(expected_result.equals(&result));
    }
    #[test]
    fn tensor_product_of_three_chained_qubit_states_yields_correct_output() {
        let starting_state = QState {
            state: vec![complex!(1.0, 0.0), C_ZERO, C_ZERO, C_ZERO],
        };
        let random_state = QState {
            state: vec![complex!(IR2, 0.0), complex!(IR2, 0.0)],
        };
        let result = starting_state.tensor_product(random_state);
        let expected_result = QState {
            state: vec![
                complex!(IR2, 0.0),
                C_ZERO,
                C_ZERO,
                C_ZERO,
                complex!(IR2, 0.0),
                C_ZERO,
                C_ZERO,
                C_ZERO,
            ],
        };
        assert!(expected_result.equals(&result));
    }
}
