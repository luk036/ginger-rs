use crate::lds::Vdcorput;
use crate::robin::Robin;
use crate::rootfinding::{horner_eval, Options};
use num::complex::Complex;
use std::convert::TryInto;
use std::f64::consts::PI;

type FoC = Complex<f64>;

fn initial_aberth_orig(pa: &Vec<FoC>) -> Vec<Complex<f64>> {
    let N: usize = pa.len() - 1;
    let c: FoC = -pa[1] / (N.try_into().unwrap() * pa[0]);
    let mut pa_copy = pa.clone();
    let Pc: FoC = horner_eval(&mut pa_copy, N, c);
    let re: f64 = (-Pc).sqrt().re;

    let mut z0s: Vec<Complex<f64>> = Vec::new();
    let mut vgen = Vdcorput::new(2);
    vgen.reseed(1);
    for _ in 0..N {
        let vdc = 2.0 * PI * vgen.pop();
        z0s.push(c + Complex::new(re, 0.0) * Complex::from_polar(&1.0, &vdc));
    }
    z0s
}

fn aberth(
    pa: &Vec<FoC>,
    zs: &mut Vec<Complex<f64>>,
    options: &Options,
) -> (Vec<Complex<f64>>, usize, bool) {
    let M = zs.len();
    let N = pa.len() - 1;
    let mut converged = vec![false; M];
    let mut robin = Robin::new(M);
    for niter in 0..options.max_iters {
        let mut tol = 0.0;
        for i in (0..M).filter(|&i| !converged[i]) {
            let mut pb = pa.clone();
            let P = horner_eval(&mut pb, N, zs[i]);
            let tol_i = P.norm();
            if tol_i < options.tol_ind {
                converged[i] = true;
                continue;
            }
            let mut P1 = horner_eval(&mut pb, N - 1, zs[i]);
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
