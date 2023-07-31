use num_traits::ops::mul_add::MulAddAssign;
use num::Complex;

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
/// use bairstow::horner::horner_eval_g;
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
/// use bairstow::horner::horner_eval_g;
/// // Here we have the "polynomial" 42, which is to say, 42𝑥⁰. Evaluated with 𝑥 = 9000
/// assert_eq!(42, horner_eval_g(9000, &[42]));
/// ```
///
/// ```
/// use bairstow::horner::horner_eval_g;
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
pub fn horner_eval_g<T: MulAddAssign + Copy> (x: T, coefficients: &[T]) -> T
{
    let (&k, coefficients) = coefficients.split_first().unwrap();
    let mut val = k;
    for &k in coefficients {
        val.mul_add_assign(x, k);
    }
    val
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
/// use bairstow::horner::const_horner_eval_g;
///
/// assert_eq!(0, const_horner_eval_g(-4, &[1, 4]));
/// ```
///
/// See also: [horner_eval_g]
pub fn const_horner_eval_g<T: MulAddAssign + Copy, const N: usize> (x: T, coefficients: &[T; N]) -> T
{
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
/// use bairstow::horner::horner_eval_c;
/// use approx_eq::assert_approx_eq;
/// use num::Complex;
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let px = horner_eval_c(&Complex::new(1.0, 2.0), &coeffs);
///
/// assert_approx_eq!(px.re, 6080.0);
/// assert_approx_eq!(px.im, 9120.0);
/// ```
///
/// See also: [const_horner_eval_c]
pub fn horner_eval_c(x: &Complex<f64>, coefficients: &[f64]) -> Complex<f64>
{
    let (&k, coefficients) = coefficients.split_first().unwrap();
    let mut val = Complex::<f64>::new(k, 0.0);
    for &k in coefficients {
        val.mul_add_assign(x, &Complex::<f64>::new(k, 0.0));
    }
    val
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
/// use bairstow::horner::const_horner_eval_c;
/// use num::Complex;
///
/// let coeffs = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let px = const_horner_eval_c(&Complex::new(1.0, 2.0), &coeffs);
///
/// assert_eq!(px.re, 6080.0);
/// assert_eq!(px.im, 9120.0);
/// ```
///
/// See also: [horner_eval_c]
pub fn const_horner_eval_c<const N: usize> (x: &Complex<f64>, coefficients: &[f64; N]) -> Complex<f64>
{
    horner_eval_c(x, coefficients)
}
