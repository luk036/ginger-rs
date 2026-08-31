use super::horner::horner_eval_f;
use super::{Matrix2, Vector2};
use crate::execution_policy::{atomic_decoupled_run, jacobi_mt_run, sequential_run, Step};
use crate::seqlock::AtomicVec2;
use num_complex::Complex;

pub type Vec2 = Vector2<f64>;
type Mat2 = Matrix2<f64>;

/// The below code defines a struct named Options with three fields: max_iters, tolerance, and tol_ind.
///
/// Properties:
///
/// * `max_iters`: The `max_iters` property represents the maximum number of iterations allowed for a
///   certain algorithm or process. It is of type `usize`, which means it can only hold non-negative
///   integer values.
/// * `tolerance`: The `tolerance` property is a floating-point number that represents the tolerance for convergence
///   in an algorithm. It is used to determine when the algorithm has reached a satisfactory solution.
/// * `tol_ind`: The `tol_ind` property in the `Options` struct represents the tolerance for individual
///   values. It is a floating-point number (`f64`) that determines the acceptable difference between the
///   expected value and the actual value for each element in a calculation or comparison.
#[derive(Debug)]
pub struct Options {
    pub max_iters: usize,
    pub tolerance: f64,
    pub tol_ind: f64,
}

/// The below code is implementing the `Default` trait for the `Options` struct in Rust. The `Default`
/// trait provides a default value for a type, which can be used when creating an instance of the type
/// without specifying any values. In this case, the `default` function is defined to return an instance
/// of the `Options` struct with default values for the `max_iters`, `tolerance`, and `tol_ind` fields.
impl Default for Options {
    fn default() -> Self {
        Options {
            max_iters: 2000,
            tolerance: 1e-12,
            tol_ind: 1e-15,
        }
    }
}

/// The function `make_adjoint` calculates the adjoint matrix between two vectors.
///
/// $$ \text{adj}\!\begin{bmatrix} r & q \\ p & s \end{bmatrix}
///    = \begin{bmatrix} s & -p \\ -pq & pr+s \end{bmatrix} $$
///
/// Arguments:
///
/// * `vr`: A vector representing the direction of the reference frame's x-axis.
/// * `vp`: The parameter `vp` represents a vector `vp = (p, s)`, where `p` and `s` are the components
///   of the vector.
///
/// Returns:
///
/// The function `make_adjoint` returns a `Mat2` object.
/// * `vp`: Another vector representing the row of a 2x2 matrix
///
/// Returns:
///
/// A 2x2 matrix representing the adjoint
#[inline]
pub fn make_adjoint(vr: &Vec2, vp: &Vec2) -> Mat2 {
    let (r, q) = (vr.x_, vr.y_);
    let (p, s) = (vp.x_, vp.y_);
    Mat2::new(
        Vector2::<f64>::new(s, -p),
        Vector2::<f64>::new(-p * q, p * r + s),
    )
}

/// The function `make_inverse` calculates the inverse of a 2x2 matrix.
///
/// $$ \mathbf{M}^{-1} = \frac{\text{adj}(\mathbf{M})}{\det(\mathbf{M})} $$
///
/// Arguments:
///
/// * `vr`: A vector representing the row of a 2x2 matrix. The components of the vector are vr.x_ and vr.y_.
/// * `vp`: The parameter `vp` represents a 2D vector with components `x` and `y`.
///
/// Returns:
///
/// The function `make_inverse` returns a `Mat2` object.
/// * `vp`: Another vector representing the row of a 2x2 matrix
///
/// Returns:
///
/// A 2x2 matrix representing the inverse
#[inline]
pub fn make_inverse(vr: &Vec2, vp: &Vec2) -> Mat2 {
    let (r, q) = (vr.x_, vr.y_);
    let (p, s) = (vp.x_, vp.y_);
    let m_adjoint = Mat2::new(
        Vector2::<f64>::new(s, -p),
        Vector2::<f64>::new(-p * q, p * r + s),
    );
    m_adjoint / m_adjoint.det()
}

/// The `delta` function calculates the adjustment vector for Bairstow's method.
///
/// Solves the 2x2 linear system for the optimal adjustment to current quadratic
/// factor estimates $$ (r, q) $$:
///
/// $$ \begin{bmatrix} r p + s & p \\ q p & s \end{bmatrix}
///    \begin{bmatrix} \Delta r \\ \Delta q \end{bmatrix}
///    = \begin{bmatrix} A \\ B \end{bmatrix} $$
///
/// where $$ (p, s) = (r_i - r_j,\; q_i - q_j) $$ is the difference between
/// two factor estimates, and $$ (A, B) $$ is the remainder from polynomial division.
///
/// Arguments:
///
/// * `vA`: A vector representing the coefficients of a polynomial equation.
/// * `vr`: The parameter `vr` represents the vector `[-2.0, 0.0]`.
/// * `vp`: The parameter `vp` represents the vector vr - vrj
///
/// # Examples:
///
/// ```
/// use ginger::rootfinding::delta;
/// use ginger::vector2::Vector2;
///
/// let mut vA1 = Vector2::new(1.0, 2.0);
/// let vri = Vector2::new(-2.0, 0.0);
/// let vrj = Vector2::new(4.0, 5.0);
/// let vd = delta(&vA1, &vri, &vrj);
/// assert_eq!(vd, Vector2::new(0.2, 0.4));
/// ```
#[inline]
pub fn delta(v_big_a: &Vec2, vr: &Vec2, vp: &Vec2) -> Vec2 {
    let mp = make_adjoint(vr, vp); // 2 mul's
    mp.mdot(v_big_a) / mp.det() // 6 mul's + 2 div's
}

/// delta 1 for ri - rj
///
/// Computes the Newton correction using the adjoint of $(vr, vp)$:
///
/// $$ \mathbf{M}_{\text{adj}} = \begin{bmatrix} -s & -p \\ pq & pr - s \end{bmatrix}, \qquad \Delta = \frac{\mathbf{M}_{\text{adj}} \cdot vA}{\det(\mathbf{M})} $$
///
/// # Examples:
///
/// ```
/// use ginger::rootfinding::delta1;
/// use ginger::vector2::Vector2;
///
/// let mut vA1 = Vector2::new(1.0, 2.0);
/// let vri = Vector2::new(-2.0, -0.0);
/// let vrj = Vector2::new(4.0, -5.0);
/// let vd = delta1(&vA1, &vri, &vrj);
/// assert_eq!(vd, Vector2::new(0.2, 0.4));
/// ```
#[inline]
pub fn delta1(v_big_a: &Vec2, vr: &Vec2, vp: &Vec2) -> Vec2 {
    let (r, q) = (vr.x_, vr.y_);
    let (p, s) = (vp.x_, vp.y_);
    let mp = Matrix2::new(Vec2::new(-s, -p), Vec2::new(p * q, p * r - s));
    mp.mdot(v_big_a) / mp.det() // 6 mul's + 2 div's
}

/// The `suppress_old` function performs zero suppression on a set of vectors.
///
/// Applies the 2x2 linear system solution using Cramer's rule:
///
/// $$ \begin{bmatrix} rp + s & p \\ qp & s \end{bmatrix} \begin{bmatrix} a \\ b \end{bmatrix} = \begin{bmatrix} A \\ B \end{bmatrix} $$
///
/// Arguments:
///
/// * `vA`: A mutable reference to a Vector2 object representing the coefficients of a polynomial. The
///   coefficients are stored in the x_ and y_ fields of the Vector2 object.
/// * `vA1`: vA1 is a mutable reference to a Vector2 object.
/// * `vri`: The parameter `vri` represents a vector with components `r` and `i`. It is used in the
///   `suppress_old` function to perform calculations.
/// * `vrj`: The parameter `vrj` represents a vector with components `x` and `y`.
///
/// # Examples:
///
/// ```
/// use ginger::rootfinding::delta;
/// use ginger::rootfinding::suppress_old;
/// use ginger::vector2::Vector2;
/// use approx_eq::assert_approx_eq;
///
/// let mut vA = Vector2::new(3.0, 3.0);
/// let mut vA1 = Vector2::new(1.0, 2.0);
/// let vri = Vector2::new(-2.0, 0.0);
/// let vrj = Vector2::new(4.0, 5.0);
///
/// suppress_old(&mut vA, &mut vA1, &vri, &vrj);
/// let dr = delta(&vA, &vri, &vA1);
/// assert_approx_eq!(dr.x_, -16.780821917808325);
/// assert_approx_eq!(dr.y_, 1.4383561643835612);
#[inline]
pub fn suppress_old(v_big_a: &mut Vec2, v_big_a1: &mut Vec2, vri: &Vec2, vrj: &Vec2) {
    let (big_a, big_b) = (v_big_a.x_, v_big_a.y_);
    let (big_a1, big_b1) = (v_big_a1.x_, v_big_a1.y_);
    let vp = vri - vrj;
    let (r, q) = (vri.x_, vri.y_);
    let (p, s) = (vp.x_, vp.y_);
    let f = (r * p) + s;
    let qp = q * p;
    let e = (f * s) - (qp * p);
    let new_a = ((big_a * s) - (big_b * p)) / e;
    let new_b = ((big_b * f) - (big_a * qp)) / e;
    let c = big_a1 - new_a;
    let d = (big_b1 - new_b) - (new_a * p);
    v_big_a.x_ = new_a;
    v_big_a.y_ = new_b;
    v_big_a1.x_ = ((c * s) - (d * p)) / e;
    v_big_a1.y_ = ((d * f) - (c * qp)) / e;
}

/// The `suppress` function in Rust performs zero suppression on a set of vectors.
///
/// Uses the inverse matrix to remove the contribution of root $j$ from the
/// remainder of root $i$:
///
/// $$ \mathbf{M}^{-1} = \frac{\text{adj}(\mathbf{vr}, \mathbf{vp})}{\det(\mathbf{M})}, \qquad \mathbf{a} = \mathbf{M}^{-1} \mathbf{vA} $$
///
/// Arguments:
///
/// * `vA`: A vector representing the coefficients of a polynomial function.
/// * `vA1`: The parameter `vA1` is a `Vector2` object representing a vector with two components. It is
///   used as an input parameter in the `suppress` function.
/// * `vri`: The parameter `vri` represents the vector `ri`, and `vrj` represents the vector `rj`. These
///   vectors are used in the calculation of the suppression step in the Bairstow's method for root
///   finding.
/// * `vrj`: The parameter `vrj` represents a vector with coordinates (4.0, 5.0).
///   Zero suppression
///
/// # Examples:
///
/// ```
/// use ginger::rootfinding::delta;
/// use ginger::rootfinding::suppress;
/// use ginger::vector2::Vector2;
/// use approx_eq::assert_approx_eq;
///
/// let mut vA = Vector2::new(3.0, 3.0);
/// let mut vA1 = Vector2::new(1.0, 2.0);
/// let vri = Vector2::new(-2.0, 0.0);
/// let vrj = Vector2::new(4.0, 5.0);
///
/// (vA, vA1) = suppress(&mut vA, &mut vA1, &vri, &vrj);
/// let dr = delta(&vA, &vri, &vA1);
/// assert_approx_eq!(dr.x_, -16.780821917808325);
/// assert_approx_eq!(dr.y_, 1.4383561643835612);
#[inline]
pub fn suppress(v_big_a: &Vec2, v_big_a1: &Vec2, vri: &Vec2, vrj: &Vec2) -> (Vec2, Vec2) {
    let vp = vri - vrj;
    let m_inverse = make_inverse(vri, &vp);
    let va = m_inverse.mdot(v_big_a);
    let mut vc = v_big_a1 - va;
    vc.y_ -= va.x_ * vp.x_;
    let va1 = m_inverse.mdot(&vc);
    (va, va1)
}

/// The `horner` function implements synthetic division by a quadratic factor $$ x^2 - r x - q $$.
///
/// Given polynomial $$ P(x) = \sum_{k=0}^{n} a_k x^{n-k} $$, the recurrence for the quotient
/// coefficients $$ b_k $$ is:
///
/// $$ b_0 = a_0,\quad b_1 = a_1 + r b_0,\quad b_k = a_k + r b_{k-1} + q b_{k-2} $$
///
/// with remainder $$ A = b_{n-1},\; B = b_n + q b_{n-1} $$.
///
/// Arguments:
///
/// * `coeffs`: A mutable slice of f64 values representing the coefficients of the polynomial. The
///   coefficients are in descending order of degree.
/// * `degree`: The `degree` parameter represents the degree of the polynomial. It is used to determine
///   the number of coefficients in the `coeffs` array.
/// * `vr`: The parameter `vr` is a `Vec2` struct that contains two values, `x_` and `y_`. In the
///   example, `vr` is initialized with the values `-1.0` and `-2.0`.
///
/// Returns:
///
/// The function `horner` returns a `Vec2` struct, which contains two `f64` values representing the
/// remainder $$ (A, B) $$ of the synthetic division.
///
/// # Examples:
///
/// ```
/// use ginger::rootfinding::horner;
/// use ginger::vector2::Vector2;
/// use approx_eq::assert_approx_eq;
///
/// let mut coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let px = horner(&mut coeffs, 8, &Vector2::new(-1.0, -2.0));
///
/// assert_approx_eq!(px.x_, 114.0);
/// assert_approx_eq!(px.y_, 134.0);
/// assert_approx_eq!(coeffs[3], 15.0);
/// ```
pub fn horner(coeffs: &mut [f64], degree: usize, vr: &Vec2) -> Vec2 {
    let Vec2 { x_: r, y_: q } = vr;
    for idx in 0..(degree - 1) {
        coeffs[idx + 1] += coeffs[idx] * r;
        coeffs[idx + 2] += coeffs[idx] * q;
    }
    Vector2::<f64>::new(coeffs[degree - 1], coeffs[degree])
}

/// The `initial_guess` function generates initial quadratic factor estimates for Bairstow's method.
///
/// Estimates are placed around a circle centered at $$ c $$ with radius $$ R $$:
///
/// $$ c = -\frac{a_1}{n a_0}, \qquad R = \sqrt\[n\]{|P(c)|} $$
///
/// where the angular positions come from a van der Corput low-discrepancy sequence.
///
/// Arguments:
///
/// * `coeffs`: A vector of coefficients representing a polynomial.
///
/// Returns:
///
/// The function `initial_guess` returns a vector of `Vector2` structs, which represent the initial
/// guesses for the roots of a polynomial equation.
///
/// # Examples:
///
/// ```
/// use ginger::rootfinding::initial_guess;
/// use ginger::vector2::Vector2;
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let vr0s = initial_guess(&coeffs);
/// ```
pub fn initial_guess(coeffs: &[f64]) -> Vec<Vec2> {
    let mut degree = coeffs.len() - 1;
    let center = -coeffs[1] / (coeffs[0] * degree as f64);
    let centroid = horner_eval_f(coeffs, center); // ???
    let radius = centroid.abs().powf(1.0 / (degree as f64));
    degree /= 2;
    degree *= 2; // make even
    let m = center * center + radius * radius;
    let num_points = degree / 2;
    (0..num_points)
        .map(|i| {
            let temp = radius * crate::tables::cos_pi_vdc2(i);
            let r0 = 2.0 * (center + temp);
            let t0 = m + 2.0 * center * temp;
            Vector2::<f64>::new(r0, -t0)
        })
        .collect()
}

/// One Bairstow Newton correction for a quadratic factor (even degree).
///
/// Reads the current factor via `get(idx)`, suppresses all other factors via
/// the neighbor iterable, and returns the corrected factor for the policy to
/// write. Shared by the sequential, Jacobi-MT and atomic execution policies.
///
/// This mirrors `ginger::detail::even_bairstow_step` in ginger-cpp.
pub struct EvenBairstowStep<'a> {
    /// Polynomial coefficients (highest degree first).
    pub coeffs: &'a [f64],
    /// Degree of the polynomial.
    pub degree: usize,
    /// Convergence options.
    pub options: &'a Options,
}

impl Step<Vec2> for EvenBairstowStep<'_> {
    type Cell = AtomicVec2;

    fn run<G, N>(&self, idx: usize, get: G, neighbors: N) -> (f64, Option<Vec2>)
    where
        G: Fn(usize) -> Vec2,
        N: Iterator<Item = usize>,
    {
        let vri = get(idx);
        let mut coeffs1 = self.coeffs.to_owned(); // horner corrupts the array
        let mut v_big_a = horner(&mut coeffs1, self.degree, &vri);
        let tol_i = v_big_a.norm_inf();
        if tol_i < self.options.tol_ind {
            return (0.0, None);
        }
        let mut v_big_a1 = horner(&mut coeffs1, self.degree - 2, &vri);
        for j in neighbors {
            let vrj = get(j);
            suppress_old(&mut v_big_a, &mut v_big_a1, &vri, &vrj);
        }
        let dt = delta(&v_big_a, &vri, &v_big_a1); // Gauss-Seidel fashion
        (tol_i, Some(vri - dt))
    }
}

/// Parallel Bairstow's method (even degree only)
///
/// The `pbairstow_even` function implements the parallel Bairstow's method for finding roots of
/// even-degree polynomials.
///
/// Arguments:
///
/// * `coeffs`: The `coeffs` parameter is a slice of `f64` values representing the coefficients of a polynomial.
///   It is assumed that the polynomial has an even degree.
/// * `vrs`: A vector of initial guesses for the roots of the polynomial. Each element of the vector is
///   a complex number representing a root guess.
/// * `options`: The `options` parameter is an instance of the `Options` struct, which contains the
///   following fields:
///
/// # Examples:
///
/// ```
/// use ginger::rootfinding::{initial_guess, pbairstow_even, Options};
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut vrs = initial_guess(&coeffs);
/// let (niter, found) = pbairstow_even(&coeffs, &mut vrs, &Options::default());
///
/// assert!(niter > 0);
/// assert!(found);
/// ```
pub fn pbairstow_even(coeffs: &[f64], vrs: &mut [Vec2], options: &Options) -> (usize, bool) {
    let step = EvenBairstowStep {
        coeffs,
        degree: coeffs.len() - 1,
        options,
    };
    sequential_run(vrs, options, &step)
}

/// Multi-threading Bairstow's method (even degree only)
///
/// The `pbairstow_even_mt` function implements the multi-threading parallel Bairstow's
/// method for finding roots of even-degree polynomials.
///
/// Arguments:
///
/// * `coeffs`: The `coeffs` parameter is a slice of `f64` values representing the coefficients of a polynomial.
///   It is assumed that the polynomial has an even degree.
/// * `vrs`: A vector of initial guesses for the roots of the polynomial. Each element of the vector is
///   a complex number representing a root guess.
/// * `options`: The `options` parameter is an instance of the `Options` struct, which contains the
///   following fields:
///
/// # Examples:
///
/// ```
/// use ginger::rootfinding::{initial_guess, pbairstow_even_mt, Options};
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut vrs = initial_guess(&coeffs);
/// let (niter, found) = pbairstow_even_mt(&coeffs, &mut vrs, &Options::default());
///
/// assert!(niter > 0);
/// assert!(found);
/// ```
pub fn pbairstow_even_mt(coeffs: &[f64], vrs: &mut [Vec2], options: &Options) -> (usize, bool) {
    let step = EvenBairstowStep {
        coeffs,
        degree: coeffs.len() - 1,
        options,
    };
    jacobi_mt_run(vrs, options, &step)
}

/// Atomic multi-threading Bairstow's method (even degree only, decoupled)
///
/// The `pbairstow_even_atomic` function implements the multi-threading parallel Bairstow's
/// method for finding roots of even-degree polynomials using a single atomic working buffer.
///
/// Unlike `pbairstow_even_mt` (Jacobi snapshot + per-iteration barrier), the atomic buffer is
/// built once: each thread owns a chunk of factor slots (single-writer, multi-reader via a
/// seqlock) and runs its own iteration loop independently. There is no per-iteration
/// synchronization — a thread exits as soon as its own chunk converges or the maximum number
/// of iterations is exceeded. The iteration count is therefore NON-DETERMINISTIC (the returned
/// count is the maximum across threads), and `found` is true only if every chunk converged.
///
/// Arguments:
///
/// * `coeffs`: The `coeffs` parameter is a slice of `f64` values representing the coefficients of a polynomial.
///   It is assumed that the polynomial has an even degree.
/// * `vrs`: A vector of initial guesses for the roots of the polynomial. Each element of the vector is
///   a complex number representing a root guess.
/// * `options`: The `options` parameter is an instance of the `Options` struct, which contains the
///   following fields:
///
/// # Examples:
///
/// ```
/// use ginger::rootfinding::{initial_guess, pbairstow_even_atomic, Options};
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut vrs = initial_guess(&coeffs);
/// let (niter, found) = pbairstow_even_atomic(&coeffs, &mut vrs, &Options::default());
///
/// assert!(niter > 0);
/// assert!(found);
/// ```
pub fn pbairstow_even_atomic(coeffs: &[f64], vrs: &mut [Vec2], options: &Options) -> (usize, bool) {
    let step = EvenBairstowStep {
        coeffs,
        degree: coeffs.len() - 1,
        options,
    };
    atomic_decoupled_run(vrs, options, &step)
}

/// Extract the two roots from a quadratic factor $$ x^2 - r x - q $$
///
/// $$ x = \frac{r \pm \sqrt{r^2 + 4q}}{2} $$
///
/// Given a quadratic factor represented as Vec2 where x() = r and y() = -q
/// (i.e., x^2 - r*x - q), return the two roots as complex numbers.
pub(crate) fn roots_from_quadratic(vr: &Vec2) -> (Complex<f64>, Complex<f64>) {
    let r = vr.x_;
    let q = vr.y_;
    let disc = r * r + 4.0 * q;
    if disc >= 0.0 {
        let sqrt_disc = disc.sqrt();
        (
            Complex::new((r + sqrt_disc) / 2.0, 0.0),
            Complex::new((r - sqrt_disc) / 2.0, 0.0),
        )
    } else {
        let sqrt_disc = (-disc).sqrt();
        (
            Complex::new(r / 2.0, sqrt_disc / 2.0),
            Complex::new(r / 2.0, -sqrt_disc / 2.0),
        )
    }
}

/// Reconstruct a monic polynomial from its quadratic factors
///
/// $$ P(x) = \prod_{i=1}^{m} (x^2 - r_i x - q_i) $$
///
/// Given the quadratic factors found by Bairstow's method (each representing
/// x^2 - r*x - q), multiply them together to recover the monic polynomial
/// coefficients. To get the original polynomial, multiply the result by the
/// original leading coefficient.
///
/// Arguments:
///
/// * `vrs` - Quadratic factors from pbairstow_even, each as a Vec2 with x() = r, y() = q
///
/// Returns:
///
/// Monic polynomial coefficients (highest degree first)
pub fn poly_from_quadratic_factors(vrs: &[Vec2]) -> Vec<f64> {
    if vrs.is_empty() {
        return vec![1.0];
    }
    // Extract all roots from quadratic factors and reconstruct with Leja ordering
    let mut all_roots: Vec<Complex<f64>> = Vec::with_capacity(2 * vrs.len());
    for vr in vrs {
        let (r1, r2) = roots_from_quadratic(vr);
        all_roots.push(r1);
        all_roots.push(r2);
    }
    crate::aberth::poly_from_roots(&all_roots)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx_eq::assert_approx_eq;

    #[test]
    fn test_options_default() {
        let options = Options::default();
        assert_eq!(options.max_iters, 2000);
        assert_eq!(options.tolerance, 1e-12);
        assert_eq!(options.tol_ind, 1e-15);
    }

    // #[test]
    // fn test_make_adjoint() {
    //     let vr = Vector2::new(1.0, 2.0);
    //     let vp = Vector2::new(3.0, 4.0);
    //     let adjoint = make_adjoint(&vr, &vp);

    //     assert_eq!(adjoint.x_.x_, 4.0);
    //     assert_eq!(adjoint.x_.y_, -3.0);
    //     assert_eq!(adjoint.y_.x_, -6.0);
    //     assert_eq!(adjoint.y_.y_, 11.0);
    // }

    // #[test]
    // fn test_make_inverse() {
    //     let vr = Vector2::new(1.0, 2.0);
    //     let vp = Vector2::new(3.0, 4.0);
    //     let inverse = make_inverse(&vr, &vp);

    //     // Verify inverse by multiplying with original matrix
    //     let original = Matrix2::new(vr, Vector2::new(vp.x_, 0.0));
    //     let product = original * inverse;
    //     assert_approx_eq!(product.x_.x_, 1.0);
    //     assert_approx_eq!(product.x_.y_, 0.0);
    //     assert_approx_eq!(product.y_.x_, 0.0);
    //     assert_approx_eq!(product.y_.y_, 1.0);
    // }

    #[test]
    fn test_delta() {
        let v_big_a = Vector2::new(1.0, 2.0);
        let vr = Vector2::new(-2.0, 0.0);
        let vp = Vector2::new(4.0, 5.0);
        let delta = delta(&v_big_a, &vr, &vp);

        assert_approx_eq!(delta.x_, 0.2);
        assert_approx_eq!(delta.y_, 0.4);
    }

    #[test]
    fn test_suppress_old() {
        let mut v_big_a = Vector2::new(3.0, 3.0);
        let mut v_big_a1 = Vector2::new(1.0, 2.0);
        let vri = Vector2::new(-2.0, 0.0);
        let vrj = Vector2::new(4.0, 5.0);

        suppress_old(&mut v_big_a, &mut v_big_a1, &vri, &vrj);
        let dr = delta(&v_big_a, &vri, &v_big_a1);
        assert_approx_eq!(dr.x_, -16.780821917808325);
        assert_approx_eq!(dr.y_, 1.4383561643835612);
    }

    #[test]
    fn test_suppress() {
        let v_big_a = Vector2::new(3.0, 3.0);
        let v_big_a1 = Vector2::new(1.0, 2.0);
        let vri = Vector2::new(-2.0, 0.0);
        let vrj = Vector2::new(4.0, 5.0);

        let (va, va1) = suppress(&v_big_a, &v_big_a1, &vri, &vrj);
        let dr = delta(&va, &vri, &va1);
        assert_approx_eq!(dr.x_, -16.780821917808325);
        assert_approx_eq!(dr.y_, 1.4383561643835612);
    }

    #[test]
    fn test_horner_eval() {
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let result = horner_eval_f(&coeffs, 2.0);

        assert_eq!(result, 18250.0);
    }

    #[test]
    fn test_horner() {
        let mut coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let vr = Vector2::new(-1.0, -2.0);
        let result = horner(&mut coeffs, 8, &vr);

        assert_approx_eq!(result.x_, 114.0);
        assert_approx_eq!(result.y_, 134.0);
        assert_eq!(coeffs[3], 15.0);
    }

    #[test]
    fn test_initial_guess() {
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let guesses = initial_guess(&coeffs);

        assert_eq!(guesses.len(), 4);
        // Verify the first guess is reasonable
        assert!(guesses[0].x_.abs() > 0.0);
        assert!(guesses[0].y_.abs() > 0.0);
    }

    #[test]
    fn test_pbairstow_even() {
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let mut vrs = initial_guess(&coeffs);
        let options = Options::default();

        let (niter, found) = pbairstow_even(&coeffs, &mut vrs, &options);

        assert!(niter > 0);
        assert!(found);
        // Verify at least one root is close to actual root
        let mut has_root = false;
        for vr in vrs {
            let val = horner(&mut coeffs.clone(), coeffs.len() - 1, &vr);
            if val.norm_inf() < options.tolerance {
                has_root = true;
                break;
            }
        }
        assert!(has_root);
    }

    #[test]
    fn test_pbairstow_even_atomic() {
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let mut vrs = initial_guess(&coeffs);
        let (niter, found) = pbairstow_even_atomic(&coeffs, &mut vrs, &Options::default());
        assert!(niter > 0);
        assert!(found);
    }

    const FIR_COEFFS: [f64; 49] = [
        -0.00196191,
        -0.00094597,
        -0.00023823,
        0.00134667,
        0.00380494,
        0.00681596,
        0.0097864,
        0.01186197,
        0.0121238,
        0.00985211,
        0.00474894,
        -0.00281751,
        -0.01173923,
        -0.0201885,
        -0.02590168,
        -0.02658216,
        -0.02035729,
        -0.00628271,
        0.01534627,
        0.04279982,
        0.0732094,
        0.10275561,
        0.12753013,
        0.14399228,
        0.15265722,
        0.14399228,
        0.12753013,
        0.10275561,
        0.0732094,
        0.04279982,
        0.01534627,
        -0.00628271,
        -0.02035729,
        -0.02658216,
        -0.02590168,
        -0.0201885,
        -0.01173923,
        -0.00281751,
        0.00474894,
        0.00985211,
        0.0121238,
        0.01186197,
        0.0097864,
        0.00681596,
        0.00380494,
        0.00134667,
        -0.00023823,
        -0.00094597,
        -0.00196191,
    ];

    #[test]
    fn test_pbairstow_even_atomic_fir() {
        // Decoupled threads under load (libtest runs tests in parallel) can need
        // many iterations to converge; give a generous iteration budget.
        let options = Options {
            max_iters: 20000,
            tolerance: 1e-2,
            ..Options::default()
        };
        let mut vrs = initial_guess(&FIR_COEFFS);
        let (_, found) = pbairstow_even_atomic(&FIR_COEFFS, &mut vrs, &options);
        assert!(found);
    }

    #[test]
    fn test_pbairstow_even_atomic_reconstruction() {
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let mut vrs = initial_guess(&coeffs);
        let (_, found) = pbairstow_even_atomic(&coeffs, &mut vrs, &Options::default());
        assert!(found);
        let monic = poly_from_quadratic_factors(&vrs);
        let scale = coeffs[0];
        for (i, c) in coeffs.iter().enumerate() {
            assert!(
                (monic[i] * scale - c).abs() < 1e-8,
                "coefficient {i} mismatch: {} vs {}",
                monic[i] * scale,
                c
            );
        }
    }

    #[test]
    fn test_delta1() {
        let v_big_a = Vector2::new(1.0, 2.0);
        let vr = Vector2::new(-2.0, -0.0);
        let vp = Vector2::new(4.0, -5.0);
        let delta = delta1(&v_big_a, &vr, &vp);

        assert_approx_eq!(delta.x_, 0.2);
        assert_approx_eq!(delta.y_, 0.4);
    }

    #[test]
    fn test_roots_from_quadratic_real() {
        // x^2 - 3x + 2 = (x-1)(x-2) -> r=3, q=-2
        let vr = Vec2::new(3.0, -2.0);
        let (r1, r2) = roots_from_quadratic(&vr);
        assert!((r1.re - 2.0).abs() < 1e-12);
        assert!((r2.re - 1.0).abs() < 1e-12);
        assert!(r1.im.abs() < 1e-12);
        assert!(r2.im.abs() < 1e-12);
    }

    #[test]
    fn test_roots_from_quadratic_complex() {
        // x^2 + 1 = 0 -> r=0, q=-1 (since x^2 - r*x - q = 0) -> roots: i, -i
        let vr = Vec2::new(0.0, -1.0);
        let (r1, r2) = roots_from_quadratic(&vr);
        assert!((r1.re).abs() < 1e-12);
        assert!((r1.im - 1.0).abs() < 1e-12);
        assert!((r2.re).abs() < 1e-12);
        assert!((r2.im + 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_poly_from_quadratic_factors() {
        // (x-1)(x-2) = x^2 - 3x + 2 -> r=3, q=-2
        let vrs = vec![Vec2::new(3.0, -2.0)];
        let coeffs = poly_from_quadratic_factors(&vrs);
        assert_eq!(coeffs.len(), 3);
        assert!((coeffs[0] - 1.0).abs() < 1e-12);
        assert!((coeffs[1] + 3.0).abs() < 1e-12);
        assert!((coeffs[2] - 2.0).abs() < 1e-12);
    }

    // ------------------------------------------------------------------
    // Order-independence verification
    // ------------------------------------------------------------------

    /// Collect all roots of the converged factors, sorted for set comparison.
    fn sorted_roots(vrs: &[Vec2]) -> Vec<(f64, f64)> {
        let mut roots = Vec::new();
        for vr in vrs {
            let (a, b) = roots_from_quadratic(vr);
            roots.push((a.re, a.im));
            roots.push((b.re, b.im));
        }
        roots.sort_by(|x, y| x.partial_cmp(y).unwrap());
        roots
    }

    /// Maximum coordinate-wise difference between two sorted root sets.
    fn max_root_set_diff(a: &[(f64, f64)], b: &[(f64, f64)]) -> f64 {
        a.iter()
            .zip(b.iter())
            .map(|(x, y)| (x.0 - y.0).abs().max((x.1 - y.1).abs()))
            .fold(0.0_f64, f64::max)
    }

    /// The Jacobi (multi-threaded) variant must be order-independent:
    /// permuting the initial guesses must give the same iteration count
    /// and the same converged root SET.
    ///
    /// NOTE: `pbairstow_even_mt` falls back to Gauss-Seidel when `vrs.len() <= 4`,
    /// so this test uses the degree-12 palindromic polynomial (6 factors) which
    /// exercises the true multi-threaded Jacobi snapshot path. Mirrors
    /// `test_order_mt.cpp` in ginger-cpp.
    #[test]
    fn test_jacobi_mt_order_independent() {
        // Palindromic degree-12: 6 factors, exercises the true Jacobi path.
        let coeffs = vec![
            1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.0, 3.0, 0.0, 2.0, 0.0, 1.0,
        ];
        let opts = Options::default();
        let base = initial_guess(&coeffs);
        assert_eq!(base.len(), 6);
        let perms = vec![
            base.clone(),
            base.iter().rev().cloned().collect(),
            vec![base[1], base[4], base[0], base[5], base[2], base[3]],
            vec![base[3], base[1], base[5], base[0], base[4], base[2]],
        ];

        let mut niters = Vec::new();
        let mut rootsets = Vec::new();
        for perm in &perms {
            let mut vrs = perm.clone();
            let (niter, found) = pbairstow_even_mt(&coeffs, &mut vrs, &opts);
            assert!(found, "Jacobi mt failed to converge");
            niters.push(niter);
            rootsets.push(sorted_roots(&vrs));
        }

        // Same iteration count for every permutation.
        assert!(
            niters.iter().all(|x| *x == niters[0]),
            "Jacobi mt iteration count is order-dependent: {niters:?}"
        );
        // Same converged root set (up to floating-point noise).
        for (k, rs) in rootsets.iter().enumerate() {
            let diff = max_root_set_diff(&rootsets[0], rs);
            assert!(
                diff < 1e-9,
                "Jacobi mt root set differs across perms: perm {k} diff {diff:.3e}"
            );
        }
    }

    /// The Jacobi sweep must produce a BIT-IDENTICAL next state regardless
    /// of the order in which factor steps are processed (frozen snapshot).
    #[test]
    fn test_jacobi_processing_order_bit_identical() {
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let vrs0 = initial_guess(&coeffs);
        let options = Options::default();
        let step = EvenBairstowStep {
            coeffs: &coeffs,
            degree: coeffs.len() - 1,
            options: &options,
        };
        let orders = [
            vec![0usize, 1, 2, 3],
            vec![3usize, 2, 1, 0],
            vec![1usize, 3, 0, 2],
            vec![2usize, 0, 3, 1],
        ];

        let mut next_states = Vec::new();
        for order in &orders {
            let vrsc = vrs0.clone(); // frozen snapshot, as in pbairstow_even_mt
            let mut next = vrs0.clone();
            for &i in order {
                let (_, new_value) =
                    step.run(i, |j| vrsc[j], (0..vrs0.len()).filter(move |&j| j != i));
                if let Some(value) = new_value {
                    next[i] = value;
                }
            }
            next_states.push(next);
        }

        for (k, state) in next_states.iter().enumerate().skip(1) {
            assert_eq!(
                next_states[0], *state,
                "Jacobi next state differs for processing order {:?}",
                orders[k]
            );
        }
    }

    /// The suppression order WITHIN a single job must only cause
    /// machine-epsilon-level drift (exactly commutative in exact arithmetic).
    #[test]
    fn test_suppression_order_machine_epsilon() {
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let vrs = initial_guess(&coeffs);
        let orders = [
            vec![0usize, 1, 2, 3],
            vec![3usize, 2, 1, 0],
            vec![1usize, 3, 0, 2],
            vec![2usize, 0, 3, 1],
        ];

        let mut worst = 0.0_f64;
        for i in 0..vrs.len() {
            let mut ref_vri = None;
            for order in &orders {
                let vri = vrs[i];
                let mut coeffs1 = coeffs.clone();
                let degree = coeffs1.len() - 1;
                let mut v_big_a = horner(&mut coeffs1, degree, &vri);
                if v_big_a.norm_inf() < 1e-15 {
                    continue;
                }
                let mut v_big_a1 = horner(&mut coeffs1, degree - 2, &vri);
                for &j in order {
                    if j == i {
                        continue;
                    }
                    suppress_old(&mut v_big_a, &mut v_big_a1, &vri, &vrs[j]);
                }
                let dt = delta(&v_big_a, &vri, &v_big_a1);
                let new_vri = vri - dt;
                match ref_vri {
                    None => ref_vri = Some(new_vri),
                    Some(r) => {
                        let d = (new_vri - r).norm_inf();
                        worst = worst.max(d);
                    }
                }
            }
        }
        assert!(
            worst < 1e-12,
            "suppression-order drift too large: {worst:.3e}"
        );
    }

    /// Contrast: the single-threaded Gauss-Seidel variant IS order-dependent
    /// in its iteration count when initial guesses are permuted.
    #[test]
    fn test_gs_order_dependent_iterations() {
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let opts = Options::default();
        let base = initial_guess(&coeffs);
        let perms = vec![
            base.clone(),
            base.iter().rev().cloned().collect(),
            vec![base[1], base[3], base[0], base[2]],
            vec![base[2], base[0], base[3], base[1]],
        ];

        let mut niters = Vec::new();
        for perm in &perms {
            let mut vrs = perm.clone();
            let (niter, found) = pbairstow_even(&coeffs, &mut vrs, &opts);
            assert!(found, "Gauss-Seidel failed to converge");
            niters.push(niter);
        }

        // Gauss-Seidel is NOT order-independent: iteration counts differ.
        assert!(
            !niters.iter().all(|x| *x == niters[0]),
            "expected Gauss-Seidel iteration count to be order-dependent, got {niters:?}"
        );
    }
}
