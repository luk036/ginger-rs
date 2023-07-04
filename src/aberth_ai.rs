use crate::lds::Vdcorput;
use crate::robin::Robin;
use crate::rootfinding::{horner_eval, Options};
use num::complex::Complex;
use std::convert::TryInto;
use std::f64::consts::PI;

type FoC = Complex<f64>;

fn initial_aberth_orig(coeffs: &Vec<FoC>) -> Vec<Complex<f64>> {
    let degree: usize = coeffs.len() - 1;
    let centerenter: FoC = -coeffs[1] / (degree.try_into().unwrap() * coeffs[0]);
    let mut pa_copy = coeffs.clone();
    let Pc: FoC = horner_eval(&mut pa_copy, degree, centerenter);
    let re: f64 = (-Pc).sqrt().re;

    let mut z0s: Vec<Complex<f64>> = Vec::new();
    let mut vgen = Vdcorput::new(2);
    vgen.reseed(1);
    for _ in 0..degree {
        let vdc = 2.0 * PI * vgen.pop();
        z0s.push(centerenter + Complex::new(re, 0.0) * Complex::from_polar(&1.0, &vdc));
    }
    z0s
}

fn aberth(
    coeffs: &Vec<FoC>,
    zs: &mut Vec<Complex<f64>>,
    options: &Options,
) -> (Vec<Complex<f64>>, usize, bool) {
    let M = zs.len();
    let degree = coeffs.len() - 1;
    let mut converged = vec![false; M];
    let mut robin = Robin::new(M);
    for niter in 0..options.max_iters {
        let mut tol = 0.0;
        for i in (0..M).filter(|&i| !converged[i]) {
            let mut pb = coeffs.clone();
            let P = horner_eval(&mut pb, degree, zs[i]);
            let tol_i = P.norm();
            if tol_i < options.tol_ind {
                converged[i] = true;
                continue;
            }
            let mut P1 = horner_eval(&mut pb, degree - 1, zs[i]);
            tol = tol.max(tol_i);
            for j in robin.exclude(i) {
                P1 -= P / (zs[i] - zs[j]);
            }
            zs[i] -= P / P1;
        }
        if tol < options.tol {
            return (zs.clone(), niter.try_into().unwrap(), true);
        }
    }
    (zs.clone(), options.max_iters, false)
}
