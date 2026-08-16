//! Order-independence experiments using the REAL ginger-rs library.
//!
//! Cross-validates the Python mirror in `order_exp.py` and tests the two
//! variants directly:
//!   - `pbairstow_even_mt`: true Jacobi (frozen snapshot `vrsc`)
//!   - `pbairstow_even`:    Gauss-Seidel (live `vrs` array passed to jobs)
//!
//! Run: cargo run --example order_experiment

use ginger::rootfinding::{
    initial_guess, pbairstow_even, pbairstow_even_mt, poly_from_quadratic_factors, Options,
};
use ginger::Vector2;
use num_complex::Complex;

type Vec2 = Vector2<f64>;

fn roots_from_quadratic(vr: &Vec2) -> (Complex<f64>, Complex<f64>) {
    let r = vr.x_;
    let q = vr.y_;
    let disc = r * r + 4.0 * q;
    if disc >= 0.0 {
        let s = disc.sqrt();
        (
            Complex::new((r + s) / 2.0, 0.0),
            Complex::new((r - s) / 2.0, 0.0),
        )
    } else {
        let s = (-disc).sqrt();
        (
            Complex::new(r / 2.0, s / 2.0),
            Complex::new(r / 2.0, -s / 2.0),
        )
    }
}

/// Collect all roots from a set of converged factors, sorted for set comparison.
fn sorted_roots(vrs: &[Vec2]) -> Vec<(f64, f64)> {
    let mut roots: Vec<(f64, f64)> = Vec::new();
    for vr in vrs {
        let (a, b) = roots_from_quadratic(vr);
        roots.push((a.re, a.im));
        roots.push((b.re, b.im));
    }
    roots.sort_by(|x, y| x.partial_cmp(y).unwrap());
    roots
}

fn max_factor_diff(a: &[Vec2], b: &[Vec2]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x.x_ - y.x_).abs().max((x.y_ - y.y_).abs()))
        .fold(0.0_f64, f64::max)
}

fn main() {
    let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
    let opts = Options::default();

    println!("======================================================================");
    println!("Rust experiment: REAL ginger-rs library");
    println!("polynomial: {:?}", coeffs);
    println!("======================================================================");

    let vrs0 = initial_guess(&coeffs);
    println!("\ninitial guesses: {:?}", vrs0);

    // ---- Jacobi variant (pbairstow_even_mt) with permuted initial guesses ----
    println!("\n--- pbairstow_even_mt (JACOBI / frozen snapshot) ---");
    let perms: Vec<Vec<Vec2>> = vec![
        vrs0.clone(),
        vrs0.iter().rev().cloned().collect(),
        vec![vrs0[1], vrs0[3], vrs0[0], vrs0[2]],
        vec![vrs0[2], vrs0[0], vrs0[3], vrs0[1]],
    ];
    let mut mt_iters = Vec::new();
    let mut mt_rootsets = Vec::new();
    for (k, perm) in perms.iter().enumerate() {
        let mut vrs = perm.clone();
        let (niter, found) = pbairstow_even_mt(&coeffs, &mut vrs, &opts);
        let roots = sorted_roots(&vrs);
        mt_iters.push((niter, found));
        mt_rootsets.push(roots.clone());
        println!(
            "  perm {}: niter={:4} found={} roots={:?}",
            k, niter, found, roots
        );
    }
    let iters_same = mt_iters.iter().all(|x| *x == mt_iters[0]);
    let rootsets_same = mt_rootsets.iter().all(|x| *x == mt_rootsets[0]);
    println!(
        "  => same iteration count across perms: {}  | same root SET: {}",
        iters_same, rootsets_same
    );

    // ---- Gauss-Seidel variant (pbairstow_even) with permuted initial guesses ----
    println!("\n--- pbairstow_even (GAUSS-SEIDEL / live array) ---");
    let mut gs_iters = Vec::new();
    let mut gs_rootsets = Vec::new();
    for (k, perm) in perms.iter().enumerate() {
        let mut vrs = perm.clone();
        let (niter, found) = pbairstow_even(&coeffs, &mut vrs, &opts);
        let roots = sorted_roots(&vrs);
        gs_iters.push((niter, found));
        gs_rootsets.push(roots.clone());
        println!(
            "  perm {}: niter={:4} found={} roots={:?}",
            k, niter, found, roots
        );
    }
    let iters_same = gs_iters.iter().all(|x| *x == gs_iters[0]);
    let rootsets_same = gs_rootsets.iter().all(|x| *x == gs_rootsets[0]);
    println!(
        "  => same iteration count across perms: {}  | same root SET: {}",
        iters_same, rootsets_same
    );

    // ---- Cross-check: do mt and non-mt converge to the same root set? ----
    println!("\n--- Jacobi vs Gauss-Seidel root-set agreement ---");
    let mt_set = &mt_rootsets[0];
    let gs_set = &gs_rootsets[0];
    let sets_agree = mt_set == gs_set;
    println!("  identical root sets (Jacobi vs GS): {}", sets_agree);

    // ---- Polynomial reconstruction check (are the roots correct?) ----
    println!("\n--- sanity: reconstructed monic polynomial from Jacobi factors ---");
    let mut vrs = vrs0.clone();
    let _ = pbairstow_even_mt(&coeffs, &mut vrs, &opts);
    let rec = poly_from_quadratic_factors(&vrs);
    println!("  reconstructed monic coeffs (desc): {:?}", rec);

    // ---- Second, harder test polynomial: roots 1,2,3,4 (degree 4) ----
    println!("\n======================================================================");
    let coeffs2: Vec<f64> = vec![1.0, -10.0, 35.0, -50.0, 24.0]; // (x-1)(x-2)(x-3)(x-4)
    println!("polynomial2: {:?}", coeffs2);
    println!("======================================================================");
    let vrs0b = initial_guess(&coeffs2);
    let perms2: Vec<Vec<Vec2>> = vec![vrs0b.clone(), vrs0b.iter().rev().cloned().collect()];
    for variant in ["mt", "gs"] {
        println!("\n  --- variant: {} ---", variant);
        let mut results = Vec::new();
        for perm in &perms2 {
            let mut vrs = perm.clone();
            let r = if variant == "mt" {
                pbairstow_even_mt(&coeffs2, &mut vrs, &opts)
            } else {
                pbairstow_even(&coeffs2, &mut vrs, &opts)
            };
            let roots = sorted_roots(&vrs);
            println!("    niter={:4} found={} roots={:?}", r.0, r.1, roots);
            results.push(roots);
        }
        println!(
            "    => same root SET across perms: {}",
            results.iter().all(|x| *x == results[0])
        );
    }

    // ---- Third test: near-degenerate suppression (two close factors) ----
    println!("\n======================================================================");
    println!("polynomial3: two nearly-equal factors (x^2+2x+2)^2 * (x^2-1)");
    println!("======================================================================");
    // (x^2+2x+2)^2 * (x^2-1) = x^6 + 4x^5 + 8x^4 + 8x^3 + 3x^2 - 4x - 4
    let coeffs3: Vec<f64> = vec![1.0, 4.0, 8.0, 8.0, 3.0, -4.0, -4.0];
    let vrs0c = initial_guess(&coeffs3);
    println!("  initial guesses: {:?}", vrs0c);
    let perms3: Vec<Vec<Vec2>> = vec![vrs0c.clone(), vrs0c.iter().rev().cloned().collect()];
    for variant in ["mt", "gs"] {
        let mut results = Vec::new();
        let mut niters = Vec::new();
        for perm in &perms3 {
            let mut vrs = perm.clone();
            let r = if variant == "mt" {
                pbairstow_even_mt(&coeffs3, &mut vrs, &opts)
            } else {
                pbairstow_even(&coeffs3, &mut vrs, &opts)
            };
            let roots = sorted_roots(&vrs);
            println!(
                "  [{}] niter={:4} found={} roots={:?}",
                variant, r.0, r.1, roots
            );
            results.push(roots);
            niters.push(r.0);
        }
        println!(
            "  => same root SET: {}  | same niter: {}",
            results.iter().all(|x| *x == results[0]),
            niters.iter().all(|x| *x == niters[0])
        );
    }

    // ---- Reproducibility: run the same Jacobi config twice ----
    println!("\n--- Jacobi determinism check (same input twice) ---");
    let mut a = vrs0.clone();
    let mut b = vrs0.clone();
    let (na, fa) = pbairstow_even_mt(&coeffs, &mut a, &opts);
    let (nb, fb) = pbairstow_even_mt(&coeffs, &mut b, &opts);
    println!(
        "  niter: {} vs {}, found: {} vs {}, max factor diff: {:.3e}",
        na,
        nb,
        fa,
        fb,
        max_factor_diff(&a, &b)
    );

    println!("\n--- Gauss-Seidel determinism check (same input twice) ---");
    let mut c = vrs0.clone();
    let mut d = vrs0.clone();
    let (nc, fc) = pbairstow_even(&coeffs, &mut c, &opts);
    let (nd, fd) = pbairstow_even(&coeffs, &mut d, &opts);
    println!(
        "  niter: {} vs {}, found: {} vs {}, max factor diff: {:.3e}",
        nc,
        nd,
        fc,
        fd,
        max_factor_diff(&c, &d)
    );
}
