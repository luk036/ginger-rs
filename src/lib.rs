#![allow(non_snake_case)]

pub mod aberth;
pub mod horner;
pub mod matrix2;
pub mod rootfinding;
pub mod vector2;
pub mod vector2_ref;
// pub mod robin;

pub use crate::aberth::{aberth, aberth_mt, initial_aberth};
pub use crate::matrix2::Matrix2;
pub use crate::rootfinding::{
    horner_eval, initial_autocorr, initial_guess, pbairstow_autocorr, pbairstow_autocorr_mt,
    pbairstow_even, pbairstow_even_mt, Options,
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
