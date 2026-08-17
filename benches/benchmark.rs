use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use ginger::{
    aberth, aberth_atomic, aberth_mt, initial_aberth, initial_autocorr, initial_guess,
    pbairstow_autocorr, pbairstow_autocorr_atomic, pbairstow_autocorr_mt, pbairstow_even,
    pbairstow_even_atomic, pbairstow_even_mt, Options,
};
use std::hint::black_box;

/// FIR autocorrelation coefficients (degree 48)
const FIR_COEFFS: [f64; 49] = [
    -0.00196191,
    -0.00094597,
    -0.00023823,
    0.00134667,
    0.00380494,
    0.00681596,
    0.0097864,
    0.01186197,
    0.0121238,
    0.00985211,
    0.00474894,
    -0.00281751,
    -0.01173923,
    -0.0201885,
    -0.02590168,
    -0.02658216,
    -0.02035729,
    -0.00628271,
    0.01534627,
    0.04279982,
    0.0732094,
    0.10275561,
    0.12753013,
    0.14399228,
    0.15265722,
    0.14399228,
    0.12753013,
    0.10275561,
    0.0732094,
    0.04279982,
    0.01534627,
    -0.00628271,
    -0.02035729,
    -0.02658216,
    -0.02590168,
    -0.0201885,
    -0.01173923,
    -0.00281751,
    0.00474894,
    0.00985211,
    0.0121238,
    0.01186197,
    0.0097864,
    0.00681596,
    0.00380494,
    0.00134667,
    -0.00023823,
    -0.00094597,
    -0.00196191,
];

fn bench_deg8(c: &mut Criterion) {
    let coeffs = black_box([10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]);
    let options = black_box(Options {
        max_iters: 2000,
        tolerance: 1e-12,
        tol_ind: 1e-15,
    });

    let mut group = c.benchmark_group("deg8");
    group.bench_function("pbairstow_even", |b| {
        let vrs = initial_guess(&coeffs);
        b.iter_batched_ref(
            || vrs.clone(),
            |vrs| pbairstow_even(&coeffs, vrs, &options),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("pbairstow_even_mt", |b| {
        let vrs = initial_guess(&coeffs);
        b.iter_batched_ref(
            || vrs.clone(),
            |vrs| pbairstow_even_mt(&coeffs, vrs, &options),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("pbairstow_even_atomic", |b| {
        let vrs = initial_guess(&coeffs);
        b.iter_batched_ref(
            || vrs.clone(),
            |vrs| pbairstow_even_atomic(&coeffs, vrs, &options),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("pbairstow_autocorr", |b| {
        let vrs = initial_autocorr(&coeffs);
        b.iter_batched_ref(
            || vrs.clone(),
            |vrs| pbairstow_autocorr(&coeffs, vrs, &options),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("pbairstow_autocorr_mt", |b| {
        let vrs = initial_autocorr(&coeffs);
        b.iter_batched_ref(
            || vrs.clone(),
            |vrs| pbairstow_autocorr_mt(&coeffs, vrs, &options),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("pbairstow_autocorr_atomic", |b| {
        let vrs = initial_autocorr(&coeffs);
        b.iter_batched_ref(
            || vrs.clone(),
            |vrs| pbairstow_autocorr_atomic(&coeffs, vrs, &options),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("aberth", |b| {
        let zs = initial_aberth(&coeffs);
        b.iter_batched_ref(
            || zs.clone(),
            |zs| aberth(&coeffs, zs, &options),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("aberth_mt", |b| {
        let zs = initial_aberth(&coeffs);
        b.iter_batched_ref(
            || zs.clone(),
            |zs| aberth_mt(&coeffs, zs, &options),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("aberth_atomic", |b| {
        let zs = initial_aberth(&coeffs);
        b.iter_batched_ref(
            || zs.clone(),
            |zs| aberth_atomic(&coeffs, zs, &options),
            BatchSize::SmallInput,
        )
    });
    group.finish();
}

fn bench_fir(c: &mut Criterion) {
    let coeffs = black_box(FIR_COEFFS);

    // FIR pbairstow variants use tolerance 1e-2 (like C++ BM_fir.cpp)
    let fir_options = black_box(Options {
        max_iters: 2000,
        tolerance: 1e-2,
        tol_ind: 1e-15,
    });

    let mut group = c.benchmark_group("fir_deg48");
    group.bench_function("pbairstow_even", |b| {
        let vrs = initial_guess(&coeffs);
        b.iter_batched_ref(
            || vrs.clone(),
            |vrs| pbairstow_even(&coeffs, vrs, &fir_options),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("pbairstow_even_mt", |b| {
        let vrs = initial_guess(&coeffs);
        b.iter_batched_ref(
            || vrs.clone(),
            |vrs| pbairstow_even_mt(&coeffs, vrs, &fir_options),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("pbairstow_even_atomic", |b| {
        let vrs = initial_guess(&coeffs);
        b.iter_batched_ref(
            || vrs.clone(),
            |vrs| pbairstow_even_atomic(&coeffs, vrs, &fir_options),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("pbairstow_autocorr", |b| {
        let vrs = initial_autocorr(&coeffs);
        b.iter_batched_ref(
            || vrs.clone(),
            |vrs| pbairstow_autocorr(&coeffs, vrs, &fir_options),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("pbairstow_autocorr_mt", |b| {
        let vrs = initial_autocorr(&coeffs);
        b.iter_batched_ref(
            || vrs.clone(),
            |vrs| pbairstow_autocorr_mt(&coeffs, vrs, &fir_options),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("pbairstow_autocorr_atomic", |b| {
        let vrs = initial_autocorr(&coeffs);
        b.iter_batched_ref(
            || vrs.clone(),
            |vrs| pbairstow_autocorr_atomic(&coeffs, vrs, &fir_options),
            BatchSize::SmallInput,
        )
    });
    group.finish();

    // FIR aberth uses tolerance 1e-8 (like C++ BM_aberth.cpp)
    let aberth_options = black_box(Options {
        max_iters: 2000,
        tolerance: 1e-8,
        tol_ind: 1e-15,
    });

    let mut group = c.benchmark_group("fir_aberth");
    group.bench_function("aberth", |b| {
        let zs = initial_aberth(&coeffs);
        b.iter_batched_ref(
            || zs.clone(),
            |zs| aberth(&coeffs, zs, &aberth_options),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("aberth_mt", |b| {
        let zs = initial_aberth(&coeffs);
        b.iter_batched_ref(
            || zs.clone(),
            |zs| aberth_mt(&coeffs, zs, &aberth_options),
            BatchSize::SmallInput,
        )
    });
    group.bench_function("aberth_atomic", |b| {
        let zs = initial_aberth(&coeffs);
        b.iter_batched_ref(
            || zs.clone(),
            |zs| aberth_atomic(&coeffs, zs, &aberth_options),
            BatchSize::SmallInput,
        )
    });
    group.finish();
}

criterion_group!(benches, bench_deg8, bench_fir);
criterion_main!(benches);
