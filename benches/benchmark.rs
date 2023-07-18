use bairstow::{
    aberth, aberth_mt, initial_aberth, initial_guess, pbairstow_autocorr, pbairstow_autocorr_mt,
    pbairstow_even, pbairstow_even_mt, Options,
};
use criterion::{black_box, Criterion, criterion_group, criterion_main};

fn bench(c: &mut Criterion) {
    let coeffs = black_box([10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]);
    let vrs = initial_guess(&coeffs);
    let options = Options {
        max_iters: 2000,
        tol: 1e-14,
        tol_ind: 1e-15,
    };

    c.bench_function("pbairstow_even", |b| {
        b.iter(|| pbairstow_even(&coeffs, &mut vrs.clone(), &options))
    });
    c.bench_function("pbairstow_even_mt", |b| {
        b.iter(|| pbairstow_even_mt(&coeffs, &mut vrs.clone(), &options))
    });
    c.bench_function("pbairstow_autocorr", |b| {
        b.iter(|| pbairstow_autocorr(&coeffs, &mut vrs.clone(), &options))
    });
    c.bench_function("pbairstow_autocorr_mt", |b| {
        b.iter(|| pbairstow_autocorr_mt(&coeffs, &mut vrs.clone(), &options))
    });

    let options = Options {
        max_iters: 2000,
        tol: 1e-12,
        tol_ind: 1e-15,
    };
    let zs = initial_aberth(&coeffs);
    c.bench_function("aberth", |b| {
        b.iter(|| aberth(&coeffs, &mut zs.clone(), &options))
    });
    c.bench_function("aberth_mt", |b| {
        b.iter(|| aberth_mt(&coeffs, &mut zs.clone(), &options))
    });
}

criterion_group!(benches, bench);
criterion_main!(benches);
