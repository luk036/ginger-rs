//! Coefficient scaling that makes the convergence test scale-invariant.

/// Normalisation factor making the convergence test scale-invariant.
///
/// At a root, `|P(z)|` has a floating-point floor proportional to the
/// coefficient magnitudes. A fixed *absolute* tolerance is therefore
/// unreachable for badly scaled polynomials, and the solver burns every
/// iteration up to `max_iters` even though the roots are already accurate.
/// Dividing the coefficients by this factor (the roots are invariant under the
/// scaling `P -> P/s`) makes [`Options::tolerance`](crate::Options::tolerance)
/// a relative measure.
///
/// The factor is floored at `1.0` so the historical absolute test is never
/// made *stricter* for well-scaled polynomials.
///
/// # Examples
///
/// ```
/// use ginger::scaling::residual_scale;
///
/// assert_eq!(residual_scale(&[10.0, 150.0, 10.0]), 150.0);
/// assert_eq!(residual_scale(&[0.1, 0.2]), 1.0);
/// assert_eq!(residual_scale(&[-3.0, 2.0]), 3.0);
/// ```
#[inline]
pub fn residual_scale(coeffs: &[f64]) -> f64 {
    coeffs.iter().fold(1.0_f64, |m, c| m.max(c.abs()))
}

/// Coefficients divided by [`residual_scale`]; the roots are invariant.
///
/// # Examples
///
/// ```
/// use ginger::scaling::scaled_coeffs;
///
/// assert_eq!(scaled_coeffs(&[10.0, 150.0, 10.0]), vec![1.0 / 15.0, 1.0, 1.0 / 15.0]);
/// ```
#[inline]
pub fn scaled_coeffs(coeffs: &[f64]) -> Vec<f64> {
    let scale = residual_scale(coeffs);
    coeffs.iter().map(|c| c / scale).collect()
}

/// Derivative coefficients of a polynomial (highest degree first).
///
/// Returns `d[i] = (degree - i) * coeffs[i]`, so evaluating the result with
/// Horner's method gives `P'(z)`.
///
/// # Examples
///
/// ```
/// use ginger::scaling::derivative_coeffs;
///
/// assert_eq!(derivative_coeffs(&[2.0, 3.0, 1.0]), vec![4.0, 3.0]);
/// ```
#[inline]
pub fn derivative_coeffs(coeffs: &[f64]) -> Vec<f64> {
    let degree = coeffs.len() - 1;
    (0..degree)
        .map(|i| coeffs[i] * (degree - i) as f64)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_scale_floors_at_one() {
        assert_eq!(residual_scale(&[0.1, 0.2]), 1.0);
        assert_eq!(residual_scale(&[10.0, 150.0, 10.0]), 150.0);
        assert_eq!(residual_scale(&[-3.0, 2.0]), 3.0);
    }

    #[test]
    fn scaled_coeffs_has_unit_max_magnitude() {
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let scaled = scaled_coeffs(&coeffs);
        assert_eq!(residual_scale(&scaled), 1.0);
        for (a, b) in coeffs.iter().zip(scaled.iter()) {
            assert!((a / 150.0 - b).abs() < 1e-15);
        }
    }

    #[test]
    fn derivative_of_linear() {
        assert_eq!(derivative_coeffs(&[2.0, 3.0]), vec![2.0]);
    }
}
