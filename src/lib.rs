#![allow(non_snake_case)]

// pub mod spectral_fact;
// pub mod robin;

//! # Ginger-rs
//! 
//! This crate provides implementations of polynomial root finding algorithms.
//! 
//! ## Modules
//! 
//! * `aberth` - Implements the Aberth method for finding the roots of a polynomial.
//! * `horner` - Implements Horner's method for polynomial evaluation.
//! * `matrix2` - Implements a simple 2x2 matrix.
//! * `rootfinding` - Implements the Bairstow's method for finding the roots of a polynomial.
//! * `leja_order` - Implements the Leja ordering.
//! * `vector2` - Implements a simple 2D vector.
//! * `vector2_ref` - Implements a simple 2D vector reference.

/// This module implements the Aberth method for finding the roots of a polynomial.
pub mod aberth;

/// This module implements Horner's method for polynomial evaluation.
pub mod horner;

/// This module implements a simple 2x2 matrix.
pub mod matrix2;

/// This module implements the Bairstow's method for finding the roots of a polynomial.
pub mod rootfinding;

/// This module implements the Leja ordering.
pub mod leja_order;

/// This module implements a simple 2D vector.
pub mod vector2;

/// This module implements a simple 2D vector reference.
pub mod vector2_ref;

pub use crate::aberth::{aberth, aberth_mt, initial_aberth};
pub use crate::horner::{horner_eval_c, horner_eval_f};
pub use crate::matrix2::Matrix2;
pub use crate::rootfinding::{
    initial_autocorr, initial_guess, pbairstow_autocorr, pbairstow_autocorr_mt, pbairstow_even,
    pbairstow_even_mt, Options,
};
pub use crate::vector2::Vector2;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let a = Vector2::<f64>::new(1.2, 2.3);
        a.scale(3.4);
        a.unscale(3.4);
        println!("{:?}", a.norm_sqr());
        println!("{:?}", a.l1_norm());

        let b = Vector2::<f64>::new(3.4, 4.5);
        println!("{:?}", a + b);
        println!("{:?}", a - b);

        let mut a = Vector2::<f64>::new(4.2, 5.3);
        a += b;
        a -= b;
        a *= 3.4;
        a /= 3.4;
        println!("{:?}", -a);
        println!("{:?}", a * 3.4);
        println!("{:?}", 3.4 * a);
        println!("{:?}", a / 3.4);

        let mm = Vector2::<Vector2<f64>>::new(a, b);
        println!("{:?}", mm);

        let mm = Matrix2::<f64>::new(a, b);
        println!("{:?}", mm);

        let b = Vector2::<i32>::new(42, 53);
        println!("{:?}", b % 3);

        let options = Options {
            max_iters: 2000,
            tolerance: 1e-14,
            tol_ind: 1e-15,
        };

        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];

        let mut vrs = initial_guess(&coeffs);
        let (niter, _found) = pbairstow_even(&coeffs, &mut vrs, &options);
        println!("{niter}");

        let mut vrs = initial_guess(&coeffs);
        let (niter, _found) = pbairstow_even_mt(&coeffs, &mut vrs, &options);
        println!("{niter}");

        let mut vrs = initial_autocorr(&coeffs);
        let (niter, _found) = pbairstow_autocorr(&coeffs, &mut vrs, &options);
        println!("{niter}");

        let mut vrs = initial_autocorr(&coeffs);
        let (niter, _found) = pbairstow_autocorr_mt(&coeffs, &mut vrs, &options);
        println!("{niter}");

        let options = Options {
            max_iters: 2000,
            tolerance: 1e-12,
            tol_ind: 1e-15,
        };

        let mut zs = initial_aberth(&coeffs);
        let (niter, _found) = aberth(&coeffs, &mut zs, &options);
        println!("{niter}");

        let mut zs = initial_aberth(&coeffs);
        let (niter, _found) = aberth_mt(&coeffs, &mut zs, &options);
        println!("{niter}");
    }
}
