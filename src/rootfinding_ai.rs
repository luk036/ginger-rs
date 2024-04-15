use std::f64::consts::PI;
use crate::lds::Vdcorput;
use crate::matrix2::Matrix2;
use crate::robin::Robin;
use crate::vector2::Vector2;

const TOL: f64 = 1e-12;
const TOL_IND: f64 = 1e-15;
const MAX_ITERS: i32 = 2000;

fn delta(vA: Vector2, vr: Vector2, vp: Vector2) -> Vector2 {
    let r = vr.x;
    let q = vr.y;
    let p = vp.x;
    let s = vp.y;
    let mp = Matrix2::new(Vector2::new(s, -p), Vector2::new(-p * q, p * r + s));
    mp.mdot(vA) / mp.det()
}

fn suppress_old(vA: &mut Vector2, vA1: &mut Vector2, vri: Vector2, vrj: Vector2) {
    let A = vA.x;
    let B = vA.y;
    let A1 = vA1.x;
    let B1 = vA1.y;
    let vp = vri - vrj;
    let r = vri.x;
    let q = vri.y;
    let p = vp.x;
    let s = vp.y;
    let f = r * p + s;
    let qp = q * p;
    let e = f * s - qp * p;
    let a = A * s - B * p;
    let b = B * f - A * qp;
    let c = A1 * e - a;
    let d = B1 * e - b - a * p;
    vA._x = a * e;
    vA._y = b * e;
    vA1._x = c * s - d * p;
    vA1._y = d * f - c * qp;
}

fn suppress(vA: Vector2, vA1: Vector2, vri: Vector2, vrj: Vector2) -> (Vector2, Vector2) {
    let vp = vri - vrj;
    let r = vri.x;
    let q = vri.y;
    let p = vp.x;
    let s = vp.y;
    let m_adjoint = Matrix2::new(Vector2::new(s, -p), Vector2::new(-p * q, p * r + s));
    let e = m_adjoint.det();
    let va = m_adjoint.mdot(vA);
    let vc = vA1 * e - va;
    vc._y -= va._x * p;
    let va = va * e;
    let va1 = m_adjoint.mdot(vc);
    (va, va1)
}

fn horner_eval(coeffs: &mut [f64], degree: usize, zval: f64) -> f64 {
    for i in 0..degree {
        coeffs[i + 1] += coeffs[i] * zval;
    }
    coeffs[degree]
}

fn horner(coeffs: &[f64], degree: usize, vr: Vector2) -> Vector2 {
    let mut coeffs = coeffs.to_vec();
    for i in 0..degree - 1 {
        coeffs[i + 1] += coeffs[i] * vr.x;
        coeffs[i + 2] += coeffs[i] * vr.y;
    }
    Vector2::new(coeffs[degree - 1], coeffs[degree])
}

fn initial_guess(coeffs: &[f64]) -> Vec<Vector2> {
    let degree = coeffs.len() - 1;
    let centroid = -coeffs[1] / (degree as f64 * coeffs[0]);

    let pc = horner_eval(&mut coeffs.to_vec(), degree, centroid);
    let reff = pc.abs().powf(1.0 / degree as f64);
    let m = centroid * centroid + reff * reff;
    let mut vr0s = Vec::new();
    let mut degree = degree / 2;
    degree *= 2;

    let mut vgen = Vdcorput::new(2);
    vgen.reseed(1);
    for _ in (1..degree).step_by(2) {
        let temp = reff * (PI * vgen.pop()).cos();
        let r0 = 2.0 * (centroid + temp);
        let t0 = m + 2.0 * centroid * temp;
        vr0s.push(Vector2::new(r0, -t0));
    }
    vr0s
}

fn pbairstow_even(coeffs: &[f64], vrs: &mut [Vector2]) -> (Vec<Vector2>, i32, bool) {
    let m_rs = vrs.len();
    let degree = coeffs.len() - 1;
    let mut converged = vec![false; m_rs];
    let mut robin = Robin::new(m_rs);
    for niter in 0..MAX_ITERS {
        let mut tolerance = 0.0;

        for i in (0..m_rs).filter(|&i| !converged[i]) {
            let mut coeffs1 = coeffs.to_vec();
            let vA = horner(&coeffs1, degree, vrs[i]);
            let tol_i = vA.x.abs().max(vA.y.abs());
            if tol_i < TOL_IND {
                converged[i] = true;
                continue;
            }
            let mut vA1 = horner(&coeffs1, degree - 2, vrs[i]);
            tolerance = tolerance.max(tol_i);

            for j in robin.exclude(i) {
                let (va, va1) = suppress(vA, vA1, vrs[i], vrs[j]);
                vA1 = va1;
                coeffs1[0] = va._x;
                coeffs1[1] = va._y;
            }
            vrs[i] -= delta(vA, vrs[i], vA1);
        }
        if tolerance < TOL {
            return (vrs.to_vec(), niter, true);
        }
    }
    (vrs.to_vec(), MAX_ITERS, false)
}

fn find_rootq(vr: Vector2) -> (f64, f64) {
    let hr = vr.x / 2.0;
    let d = hr * hr + vr.y;
    if d < 0.0 {
        let x1 = hr + (-d).sqrt() * 1.0i64;
        let x2 = -vr.y / x1;
        (x1.re, x2.re)
    } else {
        let x1 = hr + d.sqrt().copysign(hr);
        let x2 = -vr.y / x1;
        (x1, x2)
    }
}

fn extract_autocorr(vr: Vector2) -> Vector2 {
    let r = vr.x;
    let q = vr.y;
    let hr = r / 2.0;
    let d = hr * hr + q;
    if d < 0.0 {
        if q < -1.0 {
            return Vector2::new(-r, 1.0) / q;
        }
    } else {
        let a1 = hr + d.sqrt().copysign(hr);
        let a2 = -q / a1;
        if a1.abs() > 1.0 {
            if a2.abs() > 1.0 {
                return Vector2::new(a1.recip() + a2.recip(), -a1.recip() * a2.recip());
            }
            return Vector2::new(a1.recip() + a2, -a1.recip() * a2);
        } else if a2.abs() > 1.0 {
            return Vector2::new(a1 + a2.recip(), -a1 * a2.recip());
        }
        return Vector2::new(a1 + a2, -a1 * a2);
    }
    vr
}

fn test_rootfind() {
    let h = vec![5.0, 2.0, 9.0, 6.0, 2.0];
    let mut vr0s = initial_guess(&h);
    let (vrs, niter, found) = pbairstow_even(&h, &mut vr0s);
    println!("{:?}", [niter, found]);
    assert!(niter <= 4);
}

