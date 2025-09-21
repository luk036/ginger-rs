use num_complex::Complex;

/// Horner evalution (float)
///
/// The `horner_eval_f` function in Rust implements the Horner's method for evaluating a polynomial with
/// given coefficients at a specific value.
///
/// Arguments:
///
/// * `coeffs`: A vector of floating-point coefficients representing a polynomial. The coefficients are
///             ordered from highest degree to lowest degree.
/// * `zval`: The `zval` parameter in the `horner_eval_f` function represents the value at which the
///             polynomial is evaluated. It is of type `f64`, which means it is a floating-point number.
///
/// Returns:
///
/// The function `horner_eval_f` returns a `f64` value, which is the result of evaluating the polynomial
/// with the given coefficients at the specified value `zval`.
pub fn horner_eval_f(coeffs: &[f64], zval: f64) -> f64 {
    coeffs.iter().fold(0.0, |acc, coeff| acc * zval + coeff)
}

/// Horner evalution (complex)
///
/// The `horner_eval_c` function in Rust implements the Horner evaluation method for complex
/// polynomials.
///
/// Arguments:
///
/// * `coeffs`: A vector of coefficients representing a polynomial. The coefficients are in descending
///             order of degree.
/// * `zval`: The `zval` parameter is a complex number that represents the value at which the polynomial is evaluated.
///
/// Returns:
///
/// The function `horner_eval_c` returns a complex number of type `Complex<f64>`.
pub fn horner_eval_c(coeffs: &[f64], zval: &Complex<f64>) -> Complex<f64> {
    coeffs
        .iter()
        .fold(Complex::<f64>::new(0.0, 0.0), |acc, coeff| {
            acc * zval + coeff
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx_eq::assert_approx_eq;

    #[test]
    fn test_horner_eval() {
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let px = horner_eval_f(&coeffs, 2.0);
        assert_approx_eq!(px, 18250.0);

        let px = horner_eval_c(&coeffs, &Complex::new(1.0, 2.0));
        assert_approx_eq!(px.re, 6080.0);
        assert_approx_eq!(px.im, 9120.0);
    }
}
