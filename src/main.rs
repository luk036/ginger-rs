mod vector2;
use crate::vector2::Vector2;

mod matrix2;
use crate::matrix2::Matrix2;

mod rootfinding;
use crate::rootfinding::{initial_guess, pbairstow_even, Options};

fn main() {
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
        tol: 1e-3,
    };
    let pa = vec![1.0, 3.0, 3.0, 4.0];
    let mut vrs = initial_guess(&pa);
    let _res = pbairstow_even(&pa, &mut vrs, &options);
}
