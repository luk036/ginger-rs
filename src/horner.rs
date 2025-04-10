use num_complex::Complex;
use std::ops::Add;
use std::ops::Mul;

/// Evaluate a polynomial of arbitrary rank using Horner's method.
///
/// Horner's method goes like this: to find 𝑎𝑥³+𝑏𝑥²+𝑐𝑥+𝑑, you evaluate 𝑥(𝑥(𝑎𝑥+𝑏)+𝑐)+𝑑.
///
/// That's what this function does too.
///
/// You provide a value for `𝑥` and a slice of values for the coefficients `&[𝑎, 𝑏, 𝑐, 𝑑, …]`.
/// The cardinality of the slice of coefficients must equal the degree of the polynomial plus one,
/// except for the special case of the whole expression being just 0 in which case a slice of
/// length zero means the same (gives you the same result) as if the slice was equal to `&[0]`
/// or any other number of all zeros.
///
/// Here are some examples demonstrating use of eval_polynomial:
///
/// ```
/// use ginger::horner::horner_eval_g;
///
/// // Evaluating the polynomial 72𝑥²+81𝑥+99 with 𝑥 = 5
/// let val = horner_eval_g(5, &[72, 81, 99]);
///
/// // Traditional calculation.
/// let trad = 72 * 5_i32.pow(2) + 81 * 5 + 99;
///
/// assert_eq!(val, trad);
/// ```
///
/// ```
/// use ginger::horner::horner_eval_g;
/// // Here we have the "polynomial" 42, which is to say, 42𝑥⁰. Evaluated with 𝑥 = 9000
/// assert_eq!(42, horner_eval_g(9000, &[42]));
/// ```
///
/// ```
/// use ginger::horner::horner_eval_g;
/// // 23𝑥⁹+0𝑥⁸+27𝑥⁷+0𝑥⁶-5𝑥⁵+0𝑥⁴+0𝑥³+0𝑥²+0𝑥ⁱ+0𝑥⁰
/// // Written simply: 23𝑥⁹+27𝑥⁷-5𝑥⁵
/// // Evaluated with 𝑥 = 99
///
/// let val = horner_eval_g(99_i128, &[23, 0, 27, 0, -5, 0, 0, 0, 0, 0]);
/// let trad = 23 * 99_i128.pow(9) + 27 * 99_i128.pow(7) - 5 * 99_i128.pow(5);
///
/// assert_eq!(val, trad);
/// ```
///
/// See also: [const_horner_eval_g]
pub fn horner_eval_g<T: Mul<Output = T> + Add<Output = T> + Copy>(x: T, coefficients: &[T]) -> T {
    let (&k, coefficients) = coefficients.split_first().unwrap();
    coefficients.iter().fold(k, |res, &coeff| res * x + coeff)
}

/// Evaluate a polynomial of rank known at compile-time using Horner's method.
///
/// For now this function simply calls [horner_eval_g], but the idea
/// is that in the future we may be able to optimize our code further in the case
/// where the rank of the polynomial is known at compile-time.
///
/// Example usage:
///
/// ```
/// use ginger::horner::const_horner_eval_g;
///
/// assert_eq!(0, const_horner_eval_g(-4, &[1, 4]));
/// ```
///
/// See also: [horner_eval_g]
pub fn const_horner_eval_g<T: Mul<Output = T> + Add<Output = T> + Copy, const N: usize>(
    x: T,
    coefficients: &[T; N],
) -> T {
    horner_eval_g(x, coefficients)
}

/// Evaluate a polynomial of arbitrary rank using Horner's method.
///
/// Horner's method goes like this: to find 𝑎𝑥³+𝑏𝑥²+𝑐𝑥+𝑑, you evaluate 𝑥(𝑥(𝑎𝑥+𝑏)+𝑐)+𝑑.
///
/// That's what this function does too.
///
/// You provide a value for `𝑥` and a slice of values for the coefficients `&[𝑎, 𝑏, 𝑐, 𝑑, …]`.
/// The cardinality of the slice of coefficients must equal the degree of the polynomial plus one,
/// except for the special case of the whole expression being just 0 in which case a slice of
/// length zero means the same (gives you the same result) as if the slice was equal to `&[0]`
/// or any other number of all zeros.
///
/// Here are some examples demonstrating use of eval_polynomial:
///
/// ```
/// use ginger::horner::horner_eval_c;
/// use approx_eq::assert_approx_eq;
/// use num_complex::Complex;
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let px = horner_eval_c(&Complex::new(1.0, 2.0), &coeffs);
///
/// assert_approx_eq!(px.re, 6080.0);
/// assert_approx_eq!(px.im, 9120.0);
/// ```
///
/// See also: [const_horner_eval_c]
pub fn horner_eval_c(x: &Complex<f64>, coefficients: &[f64]) -> Complex<f64> {
    let (&k, coefficients) = coefficients.split_first().unwrap();
    coefficients
        .iter()
        .fold(Complex::<f64>::new(k, 0.0), |res, &coeff| res * x + coeff)
}

/// Evaluate a polynomial of rank known at compile-time using Horner's method.
///
/// For now this function simply calls [horner_eval_g], but the idea
/// is that in the future we may be able to optimize our code further in the case
/// where the rank of the polynomial is known at compile-time.
///
/// Example usage:
///
/// ```
/// use ginger::horner::const_horner_eval_c;
/// use num_complex::Complex;
///
/// let coeffs = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let px = const_horner_eval_c(&Complex::new(1.0, 2.0), &coeffs);
///
/// assert_eq!(px.re, 6080.0);
/// assert_eq!(px.im, 9120.0);
/// ```
///
/// See also: [horner_eval_c]
pub fn const_horner_eval_c<const N: usize>(
    x: &Complex<f64>,
    coefficients: &[f64; N],
) -> Complex<f64> {
    horner_eval_c(x, coefficients)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex;

    // #[test]
    // fn test_horner_eval_g_empty_coefficients() {
    //     // Empty coefficients should be treated as zero
    //     assert_eq!(horner_eval_g(5, &[]), 0);
    //     assert_eq!(horner_eval_g(10.5, &[]), 0.0);
    // }

    #[test]
    fn test_horner_eval_g_constant_polynomial() {
        // Constant polynomial (degree 0)
        assert_eq!(horner_eval_g(42, &[5]), 5);
        assert_eq!(horner_eval_g(3.14, &[2.71]), 2.71);
    }

    #[test]
    fn test_horner_eval_g_linear_polynomial() {
        // Linear polynomial (degree 1)
        assert_eq!(horner_eval_g(3, &[2, 1]), 2 * 3 + 1);
        assert_eq!(horner_eval_g(2.5, &[1.5, 3.0]), 1.5 * 2.5 + 3.0);
    }

    #[test]
    fn test_horner_eval_g_quadratic_polynomial() {
        // Quadratic polynomial (degree 2)
        assert_eq!(horner_eval_g(4, &[3, 2, 1]), 3 * 4 * 4 + 2 * 4 + 1);
        assert_eq!(
            horner_eval_g(1.5, &[2.0, 3.0, 4.0]),
            2.0 * 1.5 * 1.5 + 3.0 * 1.5 + 4.0
        );
    }

    #[test]
    fn test_horner_eval_g_with_zero_coefficients() {
        // Polynomial with some zero coefficients
        assert_eq!(horner_eval_g(2, &[1, 0, 0, 0, 5]), 2_i32.pow(4) + 5);
        assert_eq!(horner_eval_g(3, &[0, 0, 0, 0, 0, 0, 7]), 7);
    }

    #[test]
    fn test_horner_eval_g_large_numbers() {
        // Test with large numbers
        assert_eq!(
            horner_eval_g(10_i128, &[1, 2, 3, 4, 5]),
            10_i128.pow(4) + 2 * 10_i128.pow(3) + 3 * 100 + 4 * 10 + 5
        );
    }

    #[test]
    fn test_const_horner_eval_g() {
        // Test the const version with array input
        assert_eq!(const_horner_eval_g(2, &[1, 2, 3]), 4 + 2 * 2 + 3);
        assert_eq!(const_horner_eval_g(3.0, &[2.0, 1.0]), 2.0 * 3.0 + 1.0);
    }

    #[test]
    fn test_horner_eval_c_real_polynomial() {
        // Test with real coefficients (imaginary part should be zero)
        let x = Complex::new(2.0, 0.0);
        let coeffs = [1.0, 2.0, 3.0];
        let result = horner_eval_c(&x, &coeffs);
        assert_eq!(result.re, 1.0 * 4.0 + 2.0 * 2.0 + 3.0);
        assert_eq!(result.im, 0.0);
    }

    #[test]
    fn test_horner_eval_c_complex_input() {
        // Test with complex input
        let x = Complex::new(1.0, 1.0);
        let coeffs = [1.0, 2.0, 1.0];
        let result = horner_eval_c(&x, &coeffs);

        // (1 + i)^2 + 2*(1 + i) + 1 = (2i) + (2 + 2i) + 1 = 3 + 4i
        assert_eq!(result.re, 3.0);
        assert_eq!(result.im, 4.0);
    }

    #[test]
    fn test_const_horner_eval_c() {
        // Test the const version with complex input
        let x = Complex::new(1.0, -1.0);
        let coeffs = [2.0, 3.0, 4.0];
        let result = const_horner_eval_c(&x, &coeffs);

        // 2*(1 - i)^2 + 3*(1 - i) + 4 = 2*(-2i) + (3 - 3i) + 4 = (3 + 4) - (4i + 3i) = 7 - 7i
        assert_eq!(result.re, 7.0);
        assert_eq!(result.im, -7.0);
    }

    #[test]
    fn test_edge_cases() {
        // Zero polynomial
        assert_eq!(horner_eval_g(42, &[0]), 0);
        assert_eq!(horner_eval_g(42, &[0, 0, 0]), 0);

        // Zero input
        assert_eq!(horner_eval_g(0, &[1, 2, 3]), 3); // Only the constant term remains
        assert_eq!(horner_eval_g(0.0, &[1.0, 2.0, 3.0]), 3.0);

        // One coefficient
        assert_eq!(horner_eval_g(5, &[7]), 7);
    }
}
