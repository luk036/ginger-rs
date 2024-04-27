use ginger::{
    aberth, aberth_mt, initial_aberth, initial_autocorr, initial_guess, pbairstow_autocorr,
    pbairstow_autocorr_mt, pbairstow_even, pbairstow_even_mt, Options,
};
use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};

fn bench(c: &mut Criterion) {
    let coeffs = black_box([10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]);
    let options = black_box(Options {
        max_iters: 2000,
        tolerance: 1e-12,
        tol_ind: 1e-15,
    });

    let vrs = initial_guess(&coeffs);
    c.bench_function("initial_guess", |b| {
        b.iter_with_large_drop(|| initial_guess(&coeffs));
    });
    c.bench_function("pbairstow_even", |b| {
        b.iter_batched_ref(
            || vrs.clone(),
            |vrs| pbairstow_even(&coeffs, vrs, &options),
            BatchSize::SmallInput,
        )
    });
    c.bench_function("pbairstow_even_mt", |b| {
        b.iter_batched_ref(
            || vrs.clone(),
            |vrs| pbairstow_even_mt(&coeffs, vrs, &options),
            BatchSize::SmallInput,
        )
    });

    let vrs = initial_autocorr(&coeffs);
    c.bench_function("initial_autocorr", |b| {
        b.iter_with_large_drop(|| initial_autocorr(&coeffs));
    });
    c.bench_function("pbairstow_autocorr", |b| {
        b.iter_batched_ref(
            || vrs.clone(),
            |vrs| pbairstow_autocorr(&coeffs, vrs, &options),
            BatchSize::SmallInput,
        )
    });
    c.bench_function("pbairstow_autocorr_mt", |b| {
        b.iter_batched_ref(
            || vrs.clone(),
            |vrs| pbairstow_autocorr_mt(&coeffs, vrs, &options),
            BatchSize::SmallInput,
        )
    });

    // let options = black_box(Options {
    //     max_iters: 2000,
    //     tolerance: 1e-12,
    //     tol_ind: 1e-15,
    // });
    let zs = initial_aberth(&coeffs);
    c.bench_function("initial_aberth", |b| {
        b.iter_with_large_drop(|| initial_aberth(&coeffs));
    });
    c.bench_function("aberth", |b| {
        b.iter_batched_ref(
            || zs.clone(),
            |zs| aberth(&coeffs, zs, &options),
            BatchSize::SmallInput,
        )
    });
    c.bench_function("aberth_mt", |b| {
        b.iter_batched_ref(
            || zs.clone(),
            |zs| aberth_mt(&coeffs, zs, &options),
            BatchSize::SmallInput,
        )
    });
}

criterion_group!(benches, bench);
criterion_main!(benches);
