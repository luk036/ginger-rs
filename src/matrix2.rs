// #![no_std]

#[cfg(any(test, feature = "std"))]
// #[cfg_attr(test, macro_use)]
// extern crate std;

// use core::fmt;
#[cfg(test)]
use core::hash;
// use core::iter::{Product, Sum};
use core::ops::{Add, Div, Mul, Neg, Rem, Sub};
// use core::str::FromStr;
// #[cfg(feature = "std")]
// use std::error::Error;

use num_traits::{Num, Zero};

// mod vector2;
use super::Vector2;

#[derive(PartialEq, Eq, Copy, Clone, Hash, Debug, Default)]
// #[repr(C)]
pub struct Matrix2<T> {
    /// Real portion of the vector2 object
    pub x_: Vector2<T>,
    /// Imaginary portion of the vector2 object
    pub y_: Vector2<T>,
}

impl<T> Matrix2<T> {
    /// Create a new Matrix2
    #[inline]
    pub const fn new(x_: Vector2<T>, y_: Vector2<T>) -> Self {
        Matrix2 { x_, y_ }
    }
}

impl<T: Clone + Num> Matrix2<T> {
    #[inline]
    pub fn det(&self) -> T {
        self.x_.cross(&self.y_)
    }

    /// Multiplies `self` by the vector2 `v`.
    #[inline]
    pub fn mdot(&self, v: &Vector2<T>) -> Vector2<T> {
        Vector2::<T>::new(self.x_.dot(v), self.y_.dot(v))
    }

    /// Multiplies `self` by the scalar `t`.
    #[inline]
    pub fn scale(&self, t: T) -> Self {
        Self::new(self.x_.clone() * t.clone(), self.y_.clone() * t)
    }

    /// Divides `self` by the scalar `t`.
    #[inline]
    pub fn unscale(&self, t: T) -> Self {
        Self::new(self.x_.clone() / t.clone(), self.y_.clone() / t)
    }
}

macro_rules! forward_xf_xf_binop {
    (impl $imp:ident, $method:ident) => {
        impl<'a, 'b, T: Clone + Num> $imp<&'b Matrix2<T>> for &'a Matrix2<T> {
            type Output = Matrix2<T>;

            #[inline]
            fn $method(self, other: &Matrix2<T>) -> Self::Output {
                self.clone().$method(other.clone())
            }
        }
    };
}

macro_rules! forward_xf_val_binop {
    (impl $imp:ident, $method:ident) => {
        impl<'a, T: Clone + Num> $imp<Matrix2<T>> for &'a Matrix2<T> {
            type Output = Matrix2<T>;

            #[inline]
            fn $method(self, other: Matrix2<T>) -> Self::Output {
                self.clone().$method(other)
            }
        }
    };
}

macro_rules! forward_val_xf_binop {
    (impl $imp:ident, $method:ident) => {
        impl<'a, T: Clone + Num> $imp<&'a Matrix2<T>> for Matrix2<T> {
            type Output = Matrix2<T>;

            #[inline]
            fn $method(self, other: &Matrix2<T>) -> Self::Output {
                self.$method(other.clone())
            }
        }
    };
}

macro_rules! forward_all_binop {
    (impl $imp:ident, $method:ident) => {
        forward_xf_xf_binop!(impl $imp, $method);
        forward_xf_val_binop!(impl $imp, $method);
        forward_val_xf_binop!(impl $imp, $method);
    };
}

// arithmetic
forward_all_binop!(impl Add, add);

// (a, b) + (c, d) == (a + c), (b + d)
impl<T: Clone + Num> Add<Matrix2<T>> for Matrix2<T> {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self::Output::new(self.x_ + other.x_, self.y_ + other.y_)
    }
}

forward_all_binop!(impl Sub, sub);

// (a, b) - (c, d) == (a - c), (b - d)
impl<T: Clone + Num> Sub<Matrix2<T>> for Matrix2<T> {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Self::Output::new(self.x_ - other.x_, self.y_ - other.y_)
    }
}

// Op Assign

mod opassign {
    use core::ops::{AddAssign, DivAssign, MulAssign, SubAssign};

    use num_traits::NumAssign;

    use crate::Matrix2;

    impl<T: Clone + NumAssign> AddAssign for Matrix2<T> {
        fn add_assign(&mut self, other: Self) {
            self.x_ += other.x_;
            self.y_ += other.y_;
        }
    }

    impl<T: Clone + NumAssign> SubAssign for Matrix2<T> {
        fn sub_assign(&mut self, other: Self) {
            self.x_ -= other.x_;
            self.y_ -= other.y_;
        }
    }

    impl<T: Clone + NumAssign> MulAssign<T> for Matrix2<T> {
        fn mul_assign(&mut self, other: T) {
            self.x_ *= other.clone();
            self.y_ *= other;
        }
    }

    impl<T: Clone + NumAssign> DivAssign<T> for Matrix2<T> {
        fn div_assign(&mut self, other: T) {
            self.x_ /= other.clone();
            self.y_ /= other;
        }
    }

    macro_rules! forward_op_assign1 {
        (impl $imp:ident, $method:ident) => {
            impl<'a, T: Clone + NumAssign> $imp<&'a Matrix2<T>> for Matrix2<T> {
                #[inline]
                fn $method(&mut self, other: &Self) {
                    self.$method(other.clone())
                }
            }
        };
    }

    macro_rules! forward_op_assign2 {
        (impl $imp:ident, $method:ident) => {
            impl<'a, T: Clone + NumAssign> $imp<&'a T> for Matrix2<T> {
                #[inline]
                fn $method(&mut self, other: &T) {
                    self.$method(other.clone())
                }
            }
        };
    }

    forward_op_assign1!(impl AddAssign, add_assign);
    forward_op_assign1!(impl SubAssign, sub_assign);
    forward_op_assign2!(impl MulAssign, mul_assign);
    forward_op_assign2!(impl DivAssign, div_assign);
}

impl<T: Clone + Num + Neg<Output = T>> Neg for Matrix2<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self::Output::new(-self.x_, -self.y_)
    }
}

impl<'a, T: Clone + Num + Neg<Output = T>> Neg for &'a Matrix2<T> {
    type Output = Matrix2<T>;

    #[inline]
    fn neg(self) -> Self::Output {
        -self.clone()
    }
}

macro_rules! scalar_arithmetic {
    (@forward $imp:ident::$method:ident for $($scalar:ident),*) => (
        impl<'a, T: Clone + Num> $imp<&'a T> for Matrix2<T> {
            type Output = Matrix2<T>;

            #[inline]
            fn $method(self, other: &T) -> Self::Output {
                self.$method(other.clone())
            }
        }
        impl<'a, T: Clone + Num> $imp<T> for &'a Matrix2<T> {
            type Output = Matrix2<T>;

            #[inline]
            fn $method(self, other: T) -> Self::Output {
                self.clone().$method(other)
            }
        }
        impl<'a, 'b, T: Clone + Num> $imp<&'a T> for &'b Matrix2<T> {
            type Output = Matrix2<T>;

            #[inline]
            fn $method(self, other: &T) -> Self::Output {
                self.clone().$method(other.clone())
            }
        }
        $(
            impl<'a> $imp<&'a Matrix2<$scalar>> for $scalar {
                type Output = Matrix2<$scalar>;

                #[inline]
                fn $method(self, other: &Matrix2<$scalar>) -> Matrix2<$scalar> {
                    self.$method(other.clone())
                }
            }
            impl<'a> $imp<Matrix2<$scalar>> for &'a $scalar {
                type Output = Matrix2<$scalar>;

                #[inline]
                fn $method(self, other: Matrix2<$scalar>) -> Matrix2<$scalar> {
                    self.clone().$method(other)
                }
            }
            impl<'a, 'b> $imp<&'a Matrix2<$scalar>> for &'b $scalar {
                type Output = Matrix2<$scalar>;

                #[inline]
                fn $method(self, other: &Matrix2<$scalar>) -> Matrix2<$scalar> {
                    self.clone().$method(other.clone())
                }
            }
        )*
    );
    ($($scalar:ident),*) => (
        scalar_arithmetic!(@forward Mul::mul for $($scalar),*);
        // scalar_arithmetic!(@forward Div::div for $($scalar),*);
        // scalar_arithmetic!(@forward Rem::rem for $($scalar),*);

        $(
            impl Mul<Matrix2<$scalar>> for $scalar {
                type Output = Matrix2<$scalar>;

                #[inline]
                fn mul(self, other: Matrix2<$scalar>) -> Self::Output {
                    Self::Output::new(self * other.x_, self * other.y_)
                }
            }

        )*
    );
}

impl<T: Clone + Num> Mul<T> for Matrix2<T> {
    type Output = Matrix2<T>;

    #[inline]
    fn mul(self, other: T) -> Self::Output {
        Self::Output::new(self.x_ * other.clone(), self.y_ * other)
    }
}

impl<T: Clone + Num> Div<T> for Matrix2<T> {
    type Output = Self;

    #[inline]
    fn div(self, other: T) -> Self::Output {
        Self::Output::new(self.x_ / other.clone(), self.y_ / other)
    }
}

impl<T: Clone + Num> Rem<T> for Matrix2<T> {
    type Output = Matrix2<T>;

    #[inline]
    fn rem(self, other: T) -> Self::Output {
        Self::Output::new(self.x_ % other.clone(), self.y_ % other)
    }
}

scalar_arithmetic!(usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64);

// constants
impl<T: Clone + Num> Zero for Matrix2<T> {
    #[inline]
    fn zero() -> Self {
        Self::new(Zero::zero(), Zero::zero())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.x_.is_zero() && self.y_.is_zero()
    }

    #[inline]
    fn set_zero(&mut self) {
        self.x_.set_zero();
        self.y_.set_zero();
    }
}

#[cfg(test)]
fn hash<T: hash::Hash>(x: &T) -> u64 {
    use std::collections::hash_map::RandomState;
    use std::hash::{BuildHasher, Hasher};
    let mut hasher = <RandomState as BuildHasher>::Hasher::new();
    x.hash(&mut hasher);
    hasher.finish()
}

#[cfg(test)]
mod test {
    #![allow(non_upper_case_globals)]

    use super::{hash, Matrix2, Vector2};
    use core::f64;
    use num_traits::Zero;

    pub const _0_0v: Vector2<f64> = Vector2 { x_: 0.0, y_: 0.0 };
    pub const _1_0v: Vector2<f64> = Vector2 { x_: 1.0, y_: 0.0 };
    pub const _1_1v: Vector2<f64> = Vector2 { x_: 1.0, y_: 1.0 };
    pub const _0_1v: Vector2<f64> = Vector2 { x_: 0.0, y_: 1.0 };
    pub const _neg1_1v: Vector2<f64> = Vector2 { x_: -1.0, y_: 1.0 };
    pub const _05_05v: Vector2<f64> = Vector2 { x_: 0.5, y_: 0.5 };
    // pub const all_consts: [Vector2<f64>; 5] = [_0_0v, _1_0v, _1_1v, _neg1_1v, _05_05v];
    pub const _4_2v: Vector2<f64> = Vector2 { x_: 4.0, y_: 2.0 };

    pub const _0_0m: Matrix2<f64> = Matrix2 {
        x_: _0_0v,
        y_: _0_0v,
    };
    pub const _1_0m: Matrix2<f64> = Matrix2 {
        x_: _1_0v,
        y_: _0_0v,
    };
    pub const _1_1m: Matrix2<f64> = Matrix2 {
        x_: _1_1v,
        y_: _1_1v,
    };
    pub const _0_1m: Matrix2<f64> = Matrix2 {
        x_: _0_0v,
        y_: _1_0v,
    };
    pub const _neg1_1m: Matrix2<f64> = Matrix2 {
        x_: _neg1_1v,
        y_: _1_0v,
    };
    pub const _05_05m: Matrix2<f64> = Matrix2 {
        x_: _05_05v,
        y_: _05_05v,
    };
    pub const all_consts: [Matrix2<f64>; 5] = [_0_0m, _1_0m, _1_1m, _neg1_1m, _05_05m];
    pub const _4_2m: Matrix2<f64> = Matrix2 {
        x_: _4_2v,
        y_: _4_2v,
    };

    #[test]
    fn test_consts() {
        // check our constants are what Matrix2::new creates
        // fn test(c: Matrix2<f64>, r: f64, i: f64) {
        //     assert_eq!(c, Matrix2::new(r, i));
        // }
        // test(_0_0v, 0.0, 0.0);
        // test(_1_0v, 1.0, 0.0);
        // test(_1_1v, 1.0, 1.0);
        // test(_neg1_1v, -1.0, 1.0);
        // test(_05_05v, 0.5, 0.5);

        assert_eq!(_0_0m, Zero::zero());
    }

    #[test]
    fn test_scale_unscale() {
        assert_eq!(_05_05m.scale(2.0), _1_1m);
        assert_eq!(_1_1m.unscale(2.0), _05_05m);
        for &c in all_consts.iter() {
            assert_eq!(c.scale(2.0).unscale(2.0), c);
        }
    }

    #[test]
    fn test_hash() {
        let u = Vector2::new(0i32, 0i32);
        let v = Vector2::new(1i32, 0i32);
        let w = Vector2::new(0i32, 1i32);
        let a = Matrix2::new(u, v);
        let b = Matrix2::new(v, w);
        let c = Matrix2::new(w, u);
        assert!(hash(&a) != hash(&b));
        assert!(hash(&b) != hash(&c));
        assert!(hash(&c) != hash(&a));
    }
}
