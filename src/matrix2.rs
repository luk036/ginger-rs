// #![no_std]
use core::ops::{Add, Div, Mul, Neg, Rem, Sub};
use num_traits::{Num, Zero};

// mod vector2;
use super::Vector2;

/// The code defines a generic struct called Matrix2 with two fields, x_ and y_, which are both of type
/// Vector2.
///
/// Properties:
///
/// * `x_`: The `x_` property represents the first row of the `Matrix2` object. It is of type
///         `Vector2<T>`, where `T` is a generic type parameter. This means that the elements of the first row
///         are stored in a `Vector2` object.
/// * `y_`: The `y_` property is a public field of type `Vector2<T>`. It represents the second row of
///         the `Matrix2` object.
#[derive(PartialEq, Eq, Copy, Clone, Hash, Debug, Default)]
// #[repr(C)]
pub struct Matrix2<T> {
    /// The first row of the Matrix2 object
    pub x_: Vector2<T>,
    /// The second row of the Matrix2 object
    pub y_: Vector2<T>,
}

impl<T> Matrix2<T> {
    /// Creates a new [`Matrix2<T>`].
    ///
    /// The `new` function creates a new [`Matrix2`] object with the given [`Vector2`] objects.
    ///
    /// Arguments:
    ///
    /// * `x_`: A vector representing the first row of the matrix.
    /// * `y_`: The parameter `y_` is a `Vector2<T>` object representing the second row of the matrix.
    ///
    /// Returns:
    ///
    /// The `new` function returns a `Matrix2<T>` object.
    ///
    /// # Examples
    ///
    /// ```
    /// use ginger::matrix2::Matrix2;
    /// use ginger::vector2::Vector2;
    ///
    /// let x = Vector2::new(3, 4);
    /// let y = Vector2::new(5, 6);
    /// assert_eq!(Matrix2::new(x, y), Matrix2 { x_: x, y_: y });
    /// ```
    #[inline]
    pub const fn new(x_: Vector2<T>, y_: Vector2<T>) -> Self {
        Matrix2 { x_, y_ }
    }
}

impl<T: Clone + Num> Matrix2<T> {
    /// Calculate the determinant of this [`Matrix2<T>`].
    ///
    /// The `det` function calculates the determinant of a 2x2 matrix.
    ///
    /// Returns:
    ///
    /// The `det()` function returns the determinant of the `Matrix2<T>`.
    ///
    /// # Examples
    ///
    /// ```
    /// use ginger::matrix2::Matrix2;
    /// use ginger::vector2::Vector2;
    ///
    /// let x = Vector2::new(3, 4);
    /// let y = Vector2::new(5, 6);
    /// let matrix2 = Matrix2::new(x, y);
    /// assert_eq!(matrix2.det(), -2);
    /// ```
    #[inline]
    pub fn det(&self) -> T {
        self.x_.cross(&self.y_)
    }

    /// Matrix-vector multiplication
    ///
    /// The `mdot` function performs matrix-vector multiplication.
    ///
    /// Arguments:
    ///
    /// * `v`: The parameter `v` is a reference to a `Vector2<T>` object, where `T` is the type of the
    ///         elements in the vector.
    ///
    /// Returns:
    ///
    /// The `mdot` function returns a `Vector2<T>` object.
    ///
    /// # Examples
    ///
    /// ```
    /// use ginger::matrix2::Matrix2;
    /// use ginger::vector2::Vector2;
    ///
    /// let x = Vector2::new(3, 4);
    /// let y = Vector2::new(5, 6);
    /// let v = &Vector2::new(1, 1);
    /// let matrix2 = Matrix2::new(x, y);
    /// assert_eq!(matrix2.mdot(v), Vector2::new(7, 11));
    /// ```
    #[inline]
    pub fn mdot(&self, v: &Vector2<T>) -> Vector2<T> {
        Vector2::<T>::new(self.x_.dot(v), self.y_.dot(v))
    }

    /// The `scale` function multiplies a matrix by a scalar.
    ///
    /// Arguments:
    ///
    /// * `alpha`: The parameter `alpha` represents the scalar value by which the matrix is multiplied.
    ///
    /// Returns:
    ///
    /// The `scale` method returns a new `Matrix2` object.
    ///
    /// # Examples
    ///
    /// ```
    /// use ginger::matrix2::Matrix2;
    /// use ginger::vector2::Vector2;
    ///
    /// let x = Vector2::new(3, 4);
    /// let y = Vector2::new(5, 6);
    /// let matrix2 = Matrix2::new(x, y);
    /// assert_eq!(matrix2.scale(10), Matrix2 { x_: Vector2::new(30, 40), y_: Vector2::new(50, 60)});
    /// ```
    #[inline]
    pub fn scale(&self, alpha: T) -> Self {
        Self::new(self.x_.clone() * alpha.clone(), self.y_.clone() * alpha)
    }

    /// The `unscale` function divides each element of a matrix by a scalar value.
    ///
    /// Arguments:
    ///
    /// * `alpha`: The parameter `alpha` is a scalar value that is used to divide each element of the matrix
    ///             by. It is used to scale down the matrix by dividing each element by `alpha`.
    ///
    /// Returns:
    ///
    /// The `unscale` method returns a new instance of `Matrix2` with the elements of `self` divided by the
    /// scalar `alpha`.
    ///
    /// # Examples
    ///
    /// ```
    /// use ginger::matrix2::Matrix2;
    /// use ginger::vector2::Vector2;
    ///
    /// let x = Vector2::new(30, 40);
    /// let y = Vector2::new(50, 60);
    /// let matrix2 = Matrix2::new(x, y);
    /// assert_eq!(matrix2.unscale(10), Matrix2 { x_: Vector2::new(3, 4), y_: Vector2::new(5, 6)});
    /// ```
    #[inline]
    pub fn unscale(&self, alpha: T) -> Self {
        Self::new(self.x_.clone() / alpha.clone(), self.y_.clone() / alpha)
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

impl<T: Clone + Num + Neg<Output = T>> Neg for &Matrix2<T> {
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

// #[cfg(test)]
// fn hash<T: hash::Hash>(x: &T) -> u64 {
//     use std::collections::hash_map::RandomState;
//     use std::hash::{BuildHasher, Hasher};
//     let mut hasher = <RandomState as BuildHasher>::Hasher::new();
//     x.hash(&mut hasher);
//     hasher.finish()
// }

#[cfg(test)]
mod test {
    #![allow(non_upper_case_globals)]

    // use super::{hash, Matrix2, Vector2};
    use super::{Matrix2, Vector2};
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

    // #[test]
    // fn test_hash() {
    //     let u = Vector2::new(0i32, 0i32);
    //     let v = Vector2::new(1i32, 0i32);
    //     let w = Vector2::new(0i32, 1i32);
    //     let a = Matrix2::new(u, v);
    //     let b = Matrix2::new(v, w);
    //     let c = Matrix2::new(w, u);
    //     assert!(hash(&a) != hash(&b));
    //     assert!(hash(&b) != hash(&c));
    //     assert!(hash(&c) != hash(&a));
    // }

    #[test]
    fn test_new() {
        let x = Vector2::new(1, 2);
        let y = Vector2::new(3, 4);
        let m = Matrix2::new(x, y);
        assert_eq!(m.x_, x);
        assert_eq!(m.y_, y);
    }

    #[test]
    fn test_det() {
        let m = Matrix2::new(Vector2::new(3, 4), Vector2::new(5, 6));
        assert_eq!(m.det(), 3*6 - 4*5);
        
        let m2 = Matrix2::new(Vector2::new(1.0, 2.0), Vector2::new(3.0, 4.0));
        assert_eq!(m2.det(), -2.0);
    }

    #[test]
    fn test_mdot() {
        let m = Matrix2::new(Vector2::new(3, 4), Vector2::new(5, 6));
        let v = Vector2::new(1, 1);
        assert_eq!(m.mdot(&v), Vector2::new(7, 11));
        
        let m2 = Matrix2::new(Vector2::new(1.0, 2.0), Vector2::new(3.0, 4.0));
        let v2 = Vector2::new(2.0, 3.0);
        assert_eq!(m2.mdot(&v2), Vector2::new(8.0, 18.0));
    }

    #[test]
    fn test_scale() {
        let m = Matrix2::new(Vector2::new(3, 4), Vector2::new(5, 6));
        assert_eq!(m.scale(2), Matrix2::new(Vector2::new(6, 8), Vector2::new(10, 12)));
        
        let m2 = Matrix2::new(Vector2::new(1.5, 2.5), Vector2::new(3.5, 4.5));
        assert_eq!(m2.scale(2.0), Matrix2::new(Vector2::new(3.0, 5.0), Vector2::new(7.0, 9.0)));
    }

    #[test]
    fn test_unscale() {
        let m = Matrix2::new(Vector2::new(30, 40), Vector2::new(50, 60));
        assert_eq!(m.unscale(10), Matrix2::new(Vector2::new(3, 4), Vector2::new(5, 6)));
        
        let m2 = Matrix2::new(Vector2::new(3.0, 5.0), Vector2::new(7.0, 9.0));
        assert_eq!(m2.unscale(2.0), Matrix2::new(Vector2::new(1.5, 2.5), Vector2::new(3.5, 4.5)));
    }

    #[test]
    fn test_add() {
        let m1 = Matrix2::new(Vector2::new(1, 2), Vector2::new(3, 4));
        let m2 = Matrix2::new(Vector2::new(5, 6), Vector2::new(7, 8));
        assert_eq!(m1 + m2, Matrix2::new(Vector2::new(6, 8), Vector2::new(10, 12)));
        
        let m3 = &m1 + &m2;
        assert_eq!(m3, Matrix2::new(Vector2::new(6, 8), Vector2::new(10, 12)));
        
        let m4 = m1 + &m2;
        assert_eq!(m4, Matrix2::new(Vector2::new(6, 8), Vector2::new(10, 12)));
    }

    #[test]
    fn test_sub() {
        let m1 = Matrix2::new(Vector2::new(5, 6), Vector2::new(7, 8));
        let m2 = Matrix2::new(Vector2::new(1, 2), Vector2::new(3, 4));
        assert_eq!(m1 - m2, Matrix2::new(Vector2::new(4, 4), Vector2::new(4, 4)));
        
        let m3 = &m1 - &m2;
        assert_eq!(m3, Matrix2::new(Vector2::new(4, 4), Vector2::new(4, 4)));
    }

    #[test]
    fn test_neg() {
        let m = Matrix2::new(Vector2::new(1, -2), Vector2::new(-3, 4));
        assert_eq!(-m, Matrix2::new(Vector2::new(-1, 2), Vector2::new(3, -4)));
        assert_eq!(-&m, Matrix2::new(Vector2::new(-1, 2), Vector2::new(3, -4)));
    }

    #[test]
    fn test_scalar_mul() {
        let m = Matrix2::new(Vector2::new(2, 3), Vector2::new(4, 5));
        assert_eq!(m * 4, Matrix2::new(Vector2::new(8, 12), Vector2::new(16, 20)));
        assert_eq!(&m * 4, Matrix2::new(Vector2::new(8, 12), Vector2::new(16, 20)));
        assert_eq!(4 * m, Matrix2::new(Vector2::new(8, 12), Vector2::new(16, 20)));
        assert_eq!(4 * &m, Matrix2::new(Vector2::new(8, 12), Vector2::new(16, 20)));
    }

    #[test]
    fn test_scalar_div() {
        let m = Matrix2::new(Vector2::new(10, 20), Vector2::new(30, 40));
        assert_eq!(m / 5, Matrix2::new(Vector2::new(2, 4), Vector2::new(6, 8)));
    }

    #[test]
    fn test_scalar_rem() {
        let m = Matrix2::new(Vector2::new(10, 21), Vector2::new(32, 43));
        assert_eq!(m % 3, Matrix2::new(Vector2::new(1, 0), Vector2::new(2, 1)));
    }

    #[test]
    fn test_zero() {
        let zero = Matrix2::<i32>::zero();
        assert_eq!(zero, Matrix2::new(Vector2::zero(), Vector2::zero()));
        assert!(zero.is_zero());
        
        let mut m = Matrix2::new(Vector2::new(1, 2), Vector2::new(3, 4));
        assert!(!m.is_zero());
        m.set_zero();
        assert!(m.is_zero());
    }

    #[test]
    fn test_add_assign() {
        let mut m1 = Matrix2::new(Vector2::new(1, 2), Vector2::new(3, 4));
        let m2 = Matrix2::new(Vector2::new(5, 6), Vector2::new(7, 8));
        m1 += m2;
        assert_eq!(m1, Matrix2::new(Vector2::new(6, 8), Vector2::new(10, 12)));
        
        let mut m3 = Matrix2::new(Vector2::new(1, 2), Vector2::new(3, 4));
        m3 += &m2;
        assert_eq!(m3, Matrix2::new(Vector2::new(6, 8), Vector2::new(10, 12)));
    }

    #[test]
    fn test_sub_assign() {
        let mut m1 = Matrix2::new(Vector2::new(5, 6), Vector2::new(7, 8));
        let m2 = Matrix2::new(Vector2::new(1, 2), Vector2::new(3, 4));
        m1 -= m2;
        assert_eq!(m1, Matrix2::new(Vector2::new(4, 4), Vector2::new(4, 4)));
        
        let mut m3 = Matrix2::new(Vector2::new(5, 6), Vector2::new(7, 8));
        m3 -= &m2;
        assert_eq!(m3, Matrix2::new(Vector2::new(4, 4), Vector2::new(4, 4)));
    }

    #[test]
    fn test_mul_assign() {
        let mut m = Matrix2::new(Vector2::new(1, 2), Vector2::new(3, 4));
        m *= 3;
        assert_eq!(m, Matrix2::new(Vector2::new(3, 6), Vector2::new(9, 12)));
        
        let mut m2 = Matrix2::new(Vector2::new(1, 2), Vector2::new(3, 4));
        let scalar = 3;
        m2 *= &scalar;
        assert_eq!(m2, Matrix2::new(Vector2::new(3, 6), Vector2::new(9, 12)));
    }

    #[test]
    fn test_div_assign() {
        let mut m = Matrix2::new(Vector2::new(6, 9), Vector2::new(12, 15));
        m /= 3;
        assert_eq!(m, Matrix2::new(Vector2::new(2, 3), Vector2::new(4, 5)));
        
        let mut m2 = Matrix2::new(Vector2::new(6, 9), Vector2::new(12, 15));
        let scalar = 3;
        m2 /= &scalar;
        assert_eq!(m2, Matrix2::new(Vector2::new(2, 3), Vector2::new(4, 5)));
    }

    #[test]
    fn test_float_operations() {
        let m = Matrix2::new(Vector2::new(1.5, 2.5), Vector2::new(3.5, 4.5));
        assert_eq!(m.scale(2.0), Matrix2::new(Vector2::new(3.0, 5.0), Vector2::new(7.0, 9.0)));
        assert_eq!(m.unscale(0.5), Matrix2::new(Vector2::new(3.0, 5.0), Vector2::new(7.0, 9.0)));
        assert_eq!(m.det(), 1.5*4.5 - 2.5*3.5);
    }

    #[test]
    fn test_clone_and_eq() {
        let m1 = Matrix2::new(Vector2::new(1, 2), Vector2::new(3, 4));
        let m2 = m1.clone();
        assert_eq!(m1, m2);
        
        let m3 = Matrix2::new(Vector2::new(2, 1), Vector2::new(4, 3));
        assert_ne!(m1, m3);
    }

    #[test]
    fn test_debug() {
        let m = Matrix2::new(Vector2::new(1, 2), Vector2::new(3, 4));
        assert_eq!(format!("{:?}", m), "Matrix2 { x_: Vector2 { x_: 1, y_: 2 }, y_: Vector2 { x_: 3, y_: 4 } }");
    }
}