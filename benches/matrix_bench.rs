use criterion::{black_box, criterion_group, criterion_main, Criterion};
use quanta::complex::Complex;
use quanta::matrix::QMatrix;
use quanta::QState;
use quanta::{constants::*, matrix};
use rand::distributions::{Distribution, Standard};
use rand::rngs::StdRng;
use rand::{random, thread_rng, Rng, SeedableRng};
impl Distribution<P> for Standard {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> P {
        P(Complex {
            re: rng.gen(),
            im: rng.gen(),
        })
    }
}
struct P(Complex);
pub fn criterion_benchmark(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(420);
    let matrix_size: usize = 256;
    let matrix_1_data = random_qdata(matrix_size.pow(2), &mut rng);
    let matrix_1 = QMatrix::from_data(matrix_1_data);

    let matrix_2_data = random_qdata(matrix_size.pow(2), &mut rng);
    let matrix_2 = QMatrix::from_data(matrix_2_data);

    c.bench_function("matrix multiplication", |b| {
        b.iter(|| black_box(matrix_1.mul(&matrix_2)))
    });
}
fn random_qdata(size: usize, rng: &mut impl Rng) -> Vec<Complex> {
    Standard
        .sample_iter(rng)
        .map(|x: P| x.0)
        .take(size)
        .collect()
}
pub fn kronecker_bench(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(69);
    let matrix_size: usize = 32;
    let matrix_1_data = random_qdata(matrix_size.pow(2), &mut rng);
    let matrix_1 = QMatrix::from_data(matrix_1_data);

    let matrix_2_data = random_qdata(matrix_size.pow(2), &mut rng);
    let matrix_2 = QMatrix::from_data(matrix_2_data);

    c.bench_function("Kronecker product with equal sized matrices", |b| {
        b.iter(|| black_box(matrix_1.kronecker(&matrix_2)))
    });
}

pub fn matrix_vector_bench(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(100);
    let dims: usize = 256;
    let matrix_data = random_qdata(dims.pow(2), &mut rng);
    let matrix = QMatrix::from_data(matrix_data);

    let state_data = random_qdata(dims, &mut rng);
    let state = QState::from_data(state_data);

    c.bench_function("matrix vector multiplication", |b| {
        b.iter(|| black_box(state.apply(&matrix)))
    });
}
criterion_group!(benches, criterion_benchmark, kronecker_bench);
criterion_main!(benches);
