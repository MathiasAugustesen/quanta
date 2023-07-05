pub mod complex;
pub mod constants;
pub mod matrix;
use std::ops::Mul;

use complex::Complex;
use constants::{ONE, ZERO};
use float_cmp::approx_eq;
use matrix::QMatrix;
use rand::{thread_rng, Rng};
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
pub enum ClassicalBit {
    Off = 0,
    On = 1,
}
impl ClassicalBit {
    // Maybe this is how the universe does it?
    pub fn from_probs(prob_0: f64, prob_1: f64) -> ClassicalBit {
        let mut rng = thread_rng();
        match rng.gen_bool(prob_1) {
            true => ClassicalBit::On,
            false => ClassicalBit::Off
        }
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
    pub fn from_classical(bit: ClassicalBit) -> Qubit {
        match bit {
            ClassicalBit::Off => Qubit::new(ONE, ZERO),
            ClassicalBit::On => Qubit::new(ZERO, ONE),
        }
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
impl Mul<Complex> for QState {
    type Output = QState;
    fn mul(self, rhs: Complex) -> Self::Output {
        let new_state: Vec<Complex> = self.state.into_iter().map(|z| z * rhs).collect();
        QState { state: new_state }
    }
}
impl QuantumVec for QState {
    fn data_slice(&self) -> &[Complex] {
        &self.state
    }
}
impl QState {
    pub fn from_data(state: Vec<Complex>) -> Self {
        QState { state }
    }
    pub fn from_qubits(qubits: &[Qubit]) -> Self {
        assert!(!qubits.is_empty());
        qubits
            .iter()
            .rev()
            .map(|&q| q.into())
            .reduce(|acc: QState, e| acc.state_tensor(&e))
            .unwrap()
    }
    /// Calculates the tensor product between two quantum states as
    /// other ⊗ self.
    pub fn state_tensor(&self, lhs: &QState) -> Self {
        QState {
            state: lhs
                .state
                .iter()
                .flat_map(|&x| self.state.iter().map(move |&y| x * y))
                .collect(),
        }
    }
    pub fn apply(&self, gate: &QMatrix) -> Self {
        QState {
            state: (0..gate.dims())
                .into_iter()
                .map(|row| {
                    let slice =
                        &gate.data_slice()[(gate.dims() * row)..(gate.dims() * row + gate.dims())];
                    slice
                        .iter()
                        .zip(self.state.iter())
                        .map(|(&x, &y)| x * y)
                        .sum::<Complex>()
                })
                .collect(),
        }
    }
    pub fn measure(&mut self, bit: usize) -> ClassicalBit {
        // TODO: Create a lookup table for the indices
        let mut prob_0_indices = Vec::with_capacity(self.state.len() / 2);
        let mut prob_1_indices = Vec::with_capacity(self.state.len() / 2);
        let mut adding_prob_0 = true;
        let mut step_counter = 0;
        let step = 2_u32.pow(bit as u32);

        for i in 0..self.state.len() {
            match adding_prob_0 {
                true => prob_0_indices.push(i),
                false => prob_1_indices.push(i),
            }
            step_counter += 1;
            if step_counter == step {
                adding_prob_0 = !adding_prob_0;
            }
        }
        let prob_0: f64 = prob_0_indices.iter().map(|x| self.state[*x].prob()).sum();
        let prob_1: f64 = prob_1_indices.iter().map(|x| self.state[*x].prob()).sum();

        let measurement_outcome = ClassicalBit::from_probs(prob_0, prob_1);
        let (delete_indices, normalize_indices) = match measurement_outcome {
            ClassicalBit::Off => (prob_1_indices, prob_0_indices),
            ClassicalBit::On => (prob_0_indices, prob_1_indices),
        };
        for i in delete_indices {
            self.state[i] = ZERO;
        }
        let normalization_factor = normalize_indices
            .iter()
            .map(|x| self.state[*x].prob())
            .sum::<f64>()
            .sqrt();
        for i in normalize_indices {
            self.state[i] /= normalization_factor;
        }
        measurement_outcome
    }
}
#[cfg(test)]
mod tests {

    use super::*;
    use constants::*;

    #[test]
    fn state_tensor_of_two_qubit_states_yields_correct_output() {
        let qubit_1 = Qubit::new(complex!(0.5, 0.5), complex!(0.5, 0.5));
        let qubit_2 = Qubit::new(complex!(IR2, IR2), complex!(0.0, 0.0));

        let qstate_1: QState = qubit_1.into();
        let qstate_2: QState = qubit_2.into();
        let result = qstate_1.state_tensor(&qstate_2);
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
    fn state_tensor_of_random_states_gives_random_state() {
        let qstate_1: QState = Qubit::new(complex!(IR2, 0.0), complex!(IR2, 0.0)).into();
        let qstate_2: QState = Qubit::new(complex!(IR2, 0.0), complex!(IR2, 0.0)).into();

        let result = qstate_1.state_tensor(&qstate_2);
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
    fn state_tensor_of_three_chained_qubit_states_yields_correct_output() {
        let starting_state = QState {
            state: vec![complex!(1.0, 0.0), ZERO, ZERO, ZERO],
        };
        let random_state = QState {
            state: vec![complex!(IR2, 0.0), complex!(IR2, 0.0)],
        };
        let result = starting_state.state_tensor(&random_state);
        let expected_result = QState {
            state: vec![
                complex!(IR2, 0.0),
                ZERO,
                ZERO,
                ZERO,
                complex!(IR2, 0.0),
                ZERO,
                ZERO,
                ZERO,
            ],
        };
        assert!(expected_result.equals(&result));
    }
    #[test]
    fn applying_hadamard_to_0_qubit_creates_random_state() {
        let starting_state = QState::from_qubits(&[Qubit::new(ONE, ZERO)]);
        let result = starting_state.apply(&H_GATE);
        let expected_result = QState {
            state: vec![C_IR2, C_IR2],
        };
        assert!(result.equals(&expected_result));
    }
    #[test]
    fn measuring_single_qubit_collapses_quantum_state() {
        let mut state = QState::from_qubits(&[
            Qubit::new(I, ZERO),
            Qubit::new(C_IR2, C_IR2),
            Qubit::new(C_IR2, C_IR2),
        ]);
        let measurement_result = state.measure(1);
        let expected_result = QState::from_qubits(&[
            Qubit::new(I, ZERO),
            Qubit::from_classical(measurement_result),
            Qubit::new(C_IR2, C_IR2),
        ]);
        dbg!(&state, &expected_result);
        assert!(state.equals(&expected_result));
    }
}
