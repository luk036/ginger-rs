#![allow(non_snake_case)]

pub mod aberth;
pub mod matrix2;
pub mod rootfinding;
pub mod vector2;

pub use crate::aberth::{aberth, aberth_th, initial_aberth};
pub use crate::matrix2::Matrix2;
pub use crate::rootfinding::{
    horner_eval, initial_autocorr, initial_guess, pbairstow_autocorr, pbairstow_autocorr_th,
    pbairstow_even, pbairstow_even_th, Options,
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
            max_iter: 2000,
            tol: 1e-14,
            tol_ind: 1e-15,
        };

        let pa = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];

        let mut vrs = initial_guess(&pa);
        let (niter, _found) = pbairstow_even(&pa, &mut vrs, &options);
        println!("{niter}");

        let mut vrs = initial_guess(&pa);
        let (niter, _found) = pbairstow_even_th(&pa, &mut vrs, &options);
        println!("{niter}");

        let mut vrs = initial_autocorr(&pa);
        let (niter, _found) = pbairstow_autocorr(&pa, &mut vrs, &options);
        println!("{niter}");

        let mut vrs = initial_autocorr(&pa);
        let (niter, _found) = pbairstow_autocorr_th(&pa, &mut vrs, &options);
        println!("{niter}");

        let options = Options {
            max_iter: 2000,
            tol: 1e-12,
            tol_ind: 1e-15,
        };

        let mut zs = initial_aberth(&pa);
        let (niter, _found) = aberth(&pa, &mut zs, &options);
        println!("{niter}");

        let mut zs = initial_aberth(&pa);
        let (niter, _found) = aberth_th(&pa, &mut zs, &options);
        println!("{niter}");
    }
}
