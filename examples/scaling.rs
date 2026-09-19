//! Degree-scaling measurements for the root-finding solvers.
//!
//! Reports iteration counts, convergence and relative residuals alongside
//! wall-clock time. These solvers cost `O(niter * n^2)`, so a change that
//! silently pushes `niter` to `max_iters` shows up as an orders-of-magnitude
//! slowdown without any change to the per-iteration work.
//!
//! Run: `cargo run --release --example scaling`

use ginger::horner::horner_eval_c;
use ginger::leja_order::leja_order;
use ginger::rootfinding::{horner, initial_guess, pbairstow_even, Options, Vec2};
use ginger::{aberth, initial_aberth, poly_from_roots};
use num_complex::Complex;
use std::time::Instant;

const DEGREES: [usize; 3] = [8, 16, 32];
const LEJA_SIZES: [usize; 4] = [50, 100, 200, 400];

/// Minimal deterministic LCG so the example needs no RNG dependency.
struct Lcg(u64);

impl Lcg {
    fn next_f64(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }
}

fn real_poly(degree: usize, seed: u64) -> Vec<f64> {
    let mut rng = Lcg(seed);
    let mut roots = Vec::with_capacity(degree);
    while roots.len() < degree {
        let re = -3.0 + 6.0 * rng.next_f64();
        let im = 0.2 + 2.8 * rng.next_f64();
        roots.push(Complex::new(re, im));
        if roots.len() < degree {
            roots.push(Complex::new(re, -im));
        }
    }
    poly_from_roots(&roots)
}

fn scale_of(coeffs: &[f64]) -> f64 {
    coeffs.iter().fold(1.0_f64, |m, c| m.max(c.abs()))
}

fn scaled(coeffs: &[f64]) -> Vec<f64> {
    let s = scale_of(coeffs);
    coeffs.iter().map(|c| c / s).collect()
}

fn poly_residual(coeffs: &[f64], zs: &[Complex<f64>]) -> f64 {
    let s = scaled(coeffs);
    zs.iter()
        .map(|z| horner_eval_c(&s, z).norm())
        .fold(0.0_f64, f64::max)
}

fn factor_residual(coeffs: &[f64], vrs: &[Vec2]) -> f64 {
    let s = scaled(coeffs);
    let degree = s.len() - 1;
    vrs.iter()
        .map(|vr| {
            let mut local = s.clone();
            horner(&mut local, degree, vr).norm_inf()
        })
        .fold(0.0_f64, f64::max)
}

fn best_ms<F: FnMut()>(mut f: F, repeats: u32) -> f64 {
    let mut best = f64::INFINITY;
    for _ in 0..repeats {
        let t0 = Instant::now();
        f();
        best = best.min(t0.elapsed().as_secs_f64() * 1e3);
    }
    best
}

fn main() {
    println!(
        "{:<10} {:>6} {:>12} {:>8} {:>6} {:>12} {:>10}",
        "solver", "degree", "max|coeff|", "niter", "conv", "rel.resid", "ms"
    );

    for &degree in &DEGREES {
        let coeffs = real_poly(degree, 42);
        let zs0 = initial_aberth(&coeffs);
        let ms = best_ms(
            || {
                let mut zs = zs0.clone();
                let _ = aberth(&coeffs, &mut zs, &Options::default());
            },
            3,
        );
        let mut zs = zs0.clone();
        let (niter, ok) = aberth(&coeffs, &mut zs, &Options::default());
        println!(
            "{:<10} {:>6} {:>12.3e} {:>8} {:>6} {:>12.3e} {:>10.3}",
            "aberth",
            degree,
            scale_of(&coeffs),
            niter,
            ok,
            poly_residual(&coeffs, &zs),
            ms
        );
    }

    for &degree in &DEGREES {
        let coeffs = real_poly(degree, 42);
        let vrs0 = initial_guess(&coeffs);
        let ms = best_ms(
            || {
                let mut vrs = vrs0.clone();
                let _ = pbairstow_even(&coeffs, &mut vrs, &Options::default());
            },
            3,
        );
        let mut vrs = vrs0.clone();
        let (niter, ok) = pbairstow_even(&coeffs, &mut vrs, &Options::default());
        println!(
            "{:<10} {:>6} {:>12.3e} {:>8} {:>6} {:>12.3e} {:>10.3}",
            "pbairstow",
            degree,
            scale_of(&coeffs),
            niter,
            ok,
            factor_residual(&coeffs, &vrs),
            ms
        );
    }

    println!();
    println!("{:<10} {:>6} {:>12}", "leja", "n", "ms");
    for &n in &LEJA_SIZES {
        let mut rng = Lcg(7);
        let pts: Vec<Complex<f64>> = (0..n)
            .map(|_| Complex::new(rng.next_f64(), rng.next_f64()))
            .collect();
        let ms = best_ms(
            || {
                let _ = leja_order(pts.clone());
            },
            1,
        );
        println!("{:<10} {:>6} {:>12.3}", "leja_order", n, ms);
    }
}
