//! Per-iteration trace of the real Rust implementation, for cross-validation
//! against the Python mirror (trace_compare.py).
//!
//! Run: cargo run --example trace_rust

use ginger::rootfinding::{initial_guess, Options};
use ginger::Vector2;

type Vec2 = Vector2<f64>;

fn job(
    coeffs: &[f64],
    i: usize,
    vri: &mut Vec2,
    converged: &mut bool,
    vrsc: &[Vec2],
) -> Option<f64> {
    let mut coeffs1 = coeffs.to_owned();
    let degree = coeffs1.len() - 1;
    let mut v_big_a = ginger::rootfinding::horner(&mut coeffs1, degree, vri);
    let tol_i = v_big_a.norm_inf();
    if tol_i < 1e-15 {
        *converged = true;
        return None;
    }
    let mut v_big_a1 = ginger::rootfinding::horner(&mut coeffs1, degree - 2, vri);
    for (_, vrj) in vrsc.iter().enumerate().filter(|t| t.0 != i) {
        ginger::rootfinding::suppress_old(&mut v_big_a, &mut v_big_a1, vri, vrj);
    }
    let dt = ginger::rootfinding::delta(&v_big_a, vri, &v_big_a1);
    *vri -= dt;
    Some(tol_i)
}

/// Same as `job` but with the CORRECT infinity norm (max of abs values),
/// isolating the effect of the signed `norm_inf` on convergence.
fn job_corrected(
    coeffs: &[f64],
    i: usize,
    vri: &mut Vec2,
    converged: &mut bool,
    vrsc: &[Vec2],
) -> Option<f64> {
    let mut coeffs1 = coeffs.to_owned();
    let degree = coeffs1.len() - 1;
    let mut v_big_a = ginger::rootfinding::horner(&mut coeffs1, degree, vri);
    let tol_i = v_big_a.x_.abs().max(v_big_a.y_.abs());
    if tol_i < 1e-15 {
        *converged = true;
        return None;
    }
    let mut v_big_a1 = ginger::rootfinding::horner(&mut coeffs1, degree - 2, vri);
    for (_, vrj) in vrsc.iter().enumerate().filter(|t| t.0 != i) {
        ginger::rootfinding::suppress_old(&mut v_big_a, &mut v_big_a1, vri, vrj);
    }
    let dt = ginger::rootfinding::delta(&v_big_a, vri, &v_big_a1);
    *vri -= dt;
    Some(tol_i)
}

fn trace_jacobi<F>(coeffs: &[f64], vrs: &mut [Vec2], max_iters: usize, job: &F) -> (usize, bool)
where
    F: Fn(&[f64], usize, &mut Vec2, &mut bool, &[Vec2]) -> Option<f64>,
{
    let m_rs = vrs.len();
    let mut converged = vec![false; m_rs];
    for niter in 1..max_iters {
        let vrsc = vrs.to_vec();
        let mut t: f64 = 0.0;
        for i in 0..m_rs {
            if converged[i] {
                continue;
            }
            let mut vri = vrs[i];
            if let Some(tol_i) = job(coeffs, i, &mut vri, &mut converged[i], &vrsc) {
                t = t.max(tol_i);
            }
            vrs[i] = vri;
        }
        println!(
            "  it {:2}: tol={:.3e}  vrs={:?}",
            niter,
            t,
            vrs.iter()
                .map(|v| (v.x_, v.y_))
                .collect::<Vec<(f64, f64)>>()
        );
        if t < 1e-12 {
            return (niter, true);
        }
    }
    (max_iters, false)
}

fn main() {
    let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
    let vrs0 = initial_guess(&coeffs);
    println!("initial guesses (full precision):");
    for v in &vrs0 {
        println!("  ({:?}, {:?})", v.x_, v.y_);
    }

    println!("\n--- STOCK library job (signed norm_inf) ---");
    let mut vrs = vrs0.clone();
    let (niter, found) = trace_jacobi(&coeffs, &mut vrs, 2000, &job);
    println!("final (stock): niter={} found={}", niter, found);
    println!(
        "  factors: {:?}",
        vrs.iter()
            .map(|v| (v.x_, v.y_))
            .collect::<Vec<(f64, f64)>>()
    );

    println!("\n--- CORRECTED job (max of abs) ---");
    let mut vrs2 = vrs0.clone();
    let (niter2, found2) = trace_jacobi(&coeffs, &mut vrs2, 2000, &job_corrected);
    println!("final (corrected): niter={} found={}", niter2, found2);
    println!(
        "  factors: {:?}",
        vrs2.iter()
            .map(|v| (v.x_, v.y_))
            .collect::<Vec<(f64, f64)>>()
    );

    let opts = Options::default();
    println!(
        "options: max_iters={} tolerance={} tol_ind={}",
        opts.max_iters, opts.tolerance, opts.tol_ind
    );
}
