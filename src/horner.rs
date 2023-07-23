use num::Complex;

use crate::Vector2;

/// Horner evalution (float)
///
/// The `horner_eval_f` function in Rust implements the Horner's method for evaluating a polynomial with
/// given coefficients at a specific value.
///
/// Arguments:
///
/// * `coeffs`: A vector of floating-point coefficients representing a polynomial. The coefficients are
/// ordered from highest degree to lowest degree. For example, the polynomial 10x^8 + 34x^7 + 75x^6 +
/// 94x^5 + 150x^4 + 94x^
/// * `zval`: The `zval` parameter in the `horner_eval_f` function represents the value at which the
/// polynomial is evaluated. It is of type `f64`, which means it is a floating-point number.
///
/// Returns:
///
/// The function `horner_eval_f` returns a `f64` value, which is the result of evaluating the polynomial
/// with the given coefficients at the specified value `zval`.
///
/// # Examples:
///
/// ```
/// use bairstow::horner::horner_eval_f;
/// use approx_eq::assert_approx_eq;
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let px = horner_eval_f(&coeffs, 2.0);
///
/// assert_approx_eq!(px, 18250.0);
/// ```
pub fn horner_eval_f(coeffs: &[f64], zval: f64) -> f64 {
    coeffs
        .iter()
        .copied()
        .reduce(|res, coeff| res * zval + coeff)
        .unwrap()
}

/// Horner evalution (complex)
///
/// The `horner_eval_c` function in Rust implements the Horner evaluation method for complex
/// polynomials.
///
/// Arguments:
///
/// * `coeffs`: A vector of coefficients representing a polynomial. The coefficients are in descending
/// order of degree. For example, the polynomial 10x^8 + 34x^7 + 75x^6 + 94x^5 + 150x^4 + 94x^3 + 75
/// * `zval`: The `zval` parameter is a complex number that represents the value at which the polynomial
/// is evaluated.
///
/// Returns:
///
/// The function `horner_eval_c` returns a complex number of type `Complex<f64>`.
///
/// # Examples:
///
/// ```
/// use bairstow::horner::horner_eval_c;
/// use approx_eq::assert_approx_eq;
/// use num::Complex;
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let px = horner_eval_c(&coeffs, &Complex::new(1.0, 2.0));
///
/// assert_approx_eq!(px.re, 6080.0);
/// assert_approx_eq!(px.im, 9120.0);
/// ```
pub fn horner_eval_c(coeffs: &[f64], zval: &Complex<f64>) -> Complex<f64> {
    coeffs
        .iter()
        .map(|coeff| Complex::<f64>::new(*coeff, 0.0))
        .reduce(|res, coeff| res * zval + coeff)
        .unwrap()
}

/// The `horner_eval` function in Rust implements the Horner's method for polynomial evaluation.
///
/// Arguments:
///
/// * `coeffs`: A mutable slice of f64 values representing the coefficients of a polynomial. The
/// coefficients are ordered from highest degree to lowest degree.
/// * `degree`: The `degree` parameter represents the degree of the polynomial. In the given example,
/// the polynomial has a degree of 8.
/// * `zval`: The `zval` parameter in the `horner_eval` function represents the value at which the
/// polynomial is evaluated. It is the value of the independent variable in the polynomial expression.
///
/// Returns:
///
/// The function `horner_eval` returns a `f64` value, which is the result of evaluating the polynomial
/// with the given coefficients at the specified value `zval`.
///
/// # Examples:
///
/// ```
/// use bairstow::horner::horner_eval;
/// use approx_eq::assert_approx_eq;
///
/// let mut coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let px = horner_eval(&mut coeffs, 8, 2.0);
///
/// assert_approx_eq!(px, 18250.0);
/// assert_approx_eq!(coeffs[3], 460.0);
/// ```
#[inline]
pub fn horner_eval(coeffs: &mut [f64], degree: usize, zval: f64) -> f64 {
    for idx in 0..degree {
        coeffs[idx + 1] += coeffs[idx] * zval;
    }
    coeffs[degree]
}

/// The `horner` function implements Horner's evaluation for Bairstow's method in Rust.
///
/// Arguments:
///
/// * `coeffs`: A mutable slice of f64 values representing the coefficients of the polynomial. The
/// coefficients are in descending order of degree.
/// * `degree`: The `degree` parameter represents the degree of the polynomial. It is used to determine
/// the number of coefficients in the `coeffs` array.
/// * `vr`: The parameter `vr` is a `Vec2` struct that contains two values, `x_` and `y_`. In the
/// example, `vr` is initialized with the values `-1.0` and `-2.0`.
///
/// Returns:
///
/// The function `horner` returns a `Vec2` struct, which contains two `f64` values representing the
/// results of the Horner evaluation.
///
/// # Examples:
///
/// ```
/// use bairstow::horner::horner;
/// use bairstow::vector2::Vector2;
/// use approx_eq::assert_approx_eq;
///
/// let mut coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let px = horner(&mut coeffs, 8, &Vector2::new(-1.0, -2.0));
///
/// assert_approx_eq!(px.x_, 114.0);
/// assert_approx_eq!(px.y_, 134.0);
/// assert_approx_eq!(coeffs[3], 15.0);
/// ```
pub fn horner(coeffs: &mut [f64], degree: usize, vr: &Vector2<f64>) -> Vector2<f64> {
    let Vector2 { x_: r, y_: q } = vr;
    for idx in 0..(degree - 1) {
        coeffs[idx + 1] += coeffs[idx] * r;
        coeffs[idx + 2] += coeffs[idx] * q;
    }
    Vector2::<f64>::new(coeffs[degree - 1], coeffs[degree])
}

/// The result of:
/// ```
/// use bairstow::horner::horner;
/// use bairstow::vector2::Vector2;
/// let mut coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let vr = Vector2::new(-1.0, -2.0);
/// let mut pb = coeffs.to_owned();
/// let degree = coeffs.len() - 1;
/// let vA = horner(&mut pb, degree, &vr);
/// let vA1 = horner(&mut pb, degree - 2, &vr);
/// drop(pb);
/// ```
/// usefull for `pbairstow_even_job` and `pbairstow_autocorr_job` to avoid a
/// clone
pub fn horner_duble_with_cheese(coeffs: &[f64], vr: Vector2<f64>) -> (Vector2<f64>, Vector2<f64>) {
    assert!(coeffs.len() > 2);
    let Vector2 { x_: r, y_: q } = vr;
    let mut r1 = [0.0, coeffs[0], coeffs[1]];
    let mut r2 = [0.0f64; 3];
    fn rotate_rs(r: &mut [f64; 3], new: f64) -> f64 {
        let old = r[0];
        r[0] = r[1];
        r[1] = r[2];
        r[2] = new;
        old
    }
    for (idx, v) in coeffs.iter().skip(2).enumerate() {
        let old = rotate_rs(&mut r1, *v);
        rotate_rs(&mut r2, old);
        r1[1] += r1[0] * r;
        r1[2] += r1[0] * q;
        if idx > 2 {
            r2[1] += r2[0] * r;
            r2[2] += r2[0] * q;
        }
    }
    rotate_rs(&mut r2, r1[0]);
    r2[1] += r2[0] * r;
    r2[2] += r2[0] * q;
    (Vector2::new(r1[1], r1[2]), Vector2::new(r2[1], r2[2]))
}
