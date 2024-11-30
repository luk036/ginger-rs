// #![no_std]
use core::ops::{Add, Div, Mul, Neg, Rem, Sub};
use num_traits::{Num, Signed, Zero};

/// The `Vector2` struct represents a 2-dimensional vector with elements of type `T`.
///
/// Properties:
///
/// * `x_`: The `x_` property represents the first element of the `Vector2` object. It is of type `T`,
///         which means it can be any type that is specified when creating an instance of `Vector2`.
/// * `y_`: The `y_` property is the second element of the `Vector2` object. It represents the
///         y-coordinate of a 2D vector.
///
/// # Examples:
///
/// ```
/// use ginger::vector2::Vector2;
///
/// assert_eq!(Vector2::new(3, 4), Vector2 { x_: 3, y_: 4});
/// ```
#[derive(PartialEq, Eq, Copy, Clone, Hash, Debug, Default)]
pub struct Vector2<T> {
    /// The first element of the vector2 object
    pub x_: T,
    /// The second element of the vector2 object
    pub y_: T,
}

impl<T> Vector2<T> {
    /// Creates a new [`Vector2<T>`].
    ///
    /// The `new` function creates a new `Vector2` instance with the given `x` and `y` values.
    ///
    /// Arguments:
    ///
    /// * `x_`: The parameter `x_` represents the x-coordinate of the vector. It is of type `T`, which means
    ///         it can be any type that implements the necessary operations for vector calculations (e.g., addition,
    ///         subtraction, multiplication).
    /// * `y_`: The `y_` parameter represents the y-coordinate of the vector.
    ///
    /// Returns:
    ///
    /// The `new` function returns a new instance of the `Vector2<T>` struct.
    ///
    /// # Examples
    ///
    /// ```
    /// use ginger::vector2::Vector2;
    ///
    /// assert_eq!(Vector2::new(3, 4), Vector2 { x_: 3, y_: 4});
    /// ```
    #[inline]
    pub const fn new(x_: T, y_: T) -> Self {
        Vector2 { x_, y_ }
    }
}

impl<T: Clone + Num> Vector2<T> {
    /// The `dot` function calculates the dot product of two vectors.
    ///
    /// Arguments:
    ///
    /// * `other`: The `other` parameter is a reference to another `Vector2` object that we want to
    ///            calculate the dot product with.
    ///
    /// Returns:
    ///
    /// The dot product of two vectors is being returned.
    ///
    /// # Examples
    ///
    /// ```
    /// use ginger::vector2::Vector2;
    ///
    /// let vector2 = &Vector2::new(3, 4);
    /// let other = &Vector2::new(5, 6);
    /// assert_eq!(vector2.dot(other), 15 + 24);
    /// assert_eq!(vector2.dot(&vector2), 9 + 16);
    /// ```
    #[inline]
    pub fn dot(&self, other: &Self) -> T {
        self.x_.clone() * other.x_.clone() + self.y_.clone() * other.y_.clone()
    }

    /// The `cross` function calculates the cross product of two vectors.
    ///
    /// Arguments:
    ///
    /// * `other`: The `other` parameter is a reference to another `Vector2` object that we want to
    ///            calculate the cross product with.
    ///
    /// Returns:
    ///
    /// The cross product of two vectors is being returned.
    ///
    /// # Examples
    ///
    /// ```
    /// use ginger::vector2::Vector2;
    ///
    /// let vector2 = &Vector2::new(3, 4);
    /// let other = &Vector2::new(5, 6);
    /// assert_eq!(vector2.cross(other), 18 - 20);
    /// assert_eq!(vector2.cross(&vector2), 0);
    /// ```
    #[inline]
    pub fn cross(&self, other: &Self) -> T {
        self.x_.clone() * other.y_.clone() - self.y_.clone() * other.x_.clone()
    }

    /// Returns the norm sqr of this [`Vector2<T>`].
    ///
    /// The `norm_sqr` function calculates the squared norm of a `Vector2` object.
    ///
    /// Returns:
    ///
    /// The `norm_sqr` function returns the squared norm of the `Vector2<T>`.
    ///
    /// # Examples
    ///
    /// ```
    /// use ginger::vector2::Vector2;
    ///
    /// let vector2 = &Vector2::new(3, 4);
    /// assert_eq!(vector2.norm_sqr(), 9 + 16);
    /// ```
    #[inline]
    pub fn norm_sqr(&self) -> T {
        self.dot(self)
    }

    /// The `scale` function multiplies the x and y components of a `Vector2` object by a given scalar
    /// value.
    ///
    /// Arguments:
    ///
    /// * `alpha`: The parameter `alpha` represents the scaling factor that will be applied to the vector.
    ///
    /// Returns:
    ///
    /// The `scale` method returns a new `Vector2` object.
    ///
    /// # Examples
    ///
    /// ```
    /// use ginger::vector2::Vector2;
    ///
    /// let vector2 = &Vector2::new(3.0, 4.0);
    /// assert_eq!(vector2.scale(10.0), Vector2::new(30.0, 40.0));
    /// assert_eq!(vector2.scale(0.5), Vector2::new(1.5, 2.0));
    /// ```
    #[inline]
    pub fn scale(&self, alpha: T) -> Self {
        Self::new(self.x_.clone() * alpha.clone(), self.y_.clone() * alpha)
    }

    /// The `unscale` function divides the x and y components of a `Vector2` by a given value.
    ///
    /// Arguments:
    ///
    /// * `alpha`: The `alpha` parameter is a value of type `T` that is used to divide the `x_` and `y_`
    ///            values of the `Vector2` object.
    ///
    /// Returns:
    ///
    /// The `unscale` method returns a new `Vector2` object.
    ///
    /// # Examples
    ///
    /// ```
    /// use ginger::vector2::Vector2;
    ///
    /// let vector2 = &Vector2::new(30, 40);
    /// assert_eq!(vector2.unscale(10), Vector2::new(3, 4));
    /// ```
    #[inline]
    pub fn unscale(&self, alpha: T) -> Self {
        Self::new(self.x_.clone() / alpha.clone(), self.y_.clone() / alpha)
    }
}

impl<T: Clone + Signed> Vector2<T> {
    /// The `l1_norm` function calculates the Manhattan distance from the origin for a 2D vector.
    ///
    /// [Manhattan distance]: https://en.wikipedia.org/wiki/Taxicab_geometry
    ///
    /// Returns:
    ///
    /// The function `l1_norm` returns the L1 norm of a `Vector2` object, which is the sum of the absolute
    /// values of its `x_` and `y_` components.
    ///
    /// # Examples
    ///
    /// ```
    /// use ginger::vector2::Vector2;
    ///
    /// let vector2 = &Vector2::new(3, -4);
    /// assert_eq!(vector2.l1_norm(), 7);
    /// ```
    #[inline]
    pub fn l1_norm(&self) -> T {
        self.x_.abs() + self.y_.abs()
    }
}

impl<T: Clone + PartialOrd> Vector2<T> {
    /// The `norm_inf` function returns the maximum absolute value of the two elements in a `Vector2`
    /// object.
    ///
    /// Returns:
    ///
    /// The `norm_inf` function returns the maximum value between `self.x_` and `self.y_`.
    ///
    /// # Examples
    ///
    /// ```
    /// use ginger::vector2::Vector2;
    ///
    /// let vector2 = &Vector2::new(3, -4);
    /// assert_eq!(vector2.norm_inf(), 3);
    /// ```
    #[inline]
    pub fn norm_inf(&self) -> T {
        if self.x_ > self.y_ {
            self.x_.clone()
        } else {
            self.y_.clone()
        }
    }
}

macro_rules! forward_xf_xf_binop {
    (impl $imp:ident, $method:ident) => {
        impl<'a, 'b, T: Clone + Num> $imp<&'b Vector2<T>> for &'a Vector2<T> {
            type Output = Vector2<T>;

            #[inline]
            fn $method(self, other: &Vector2<T>) -> Self::Output {
                self.clone().$method(other.clone())
            }
        }
    };
}

macro_rules! forward_xf_val_binop {
    (impl $imp:ident, $method:ident) => {
        impl<'a, T: Clone + Num> $imp<Vector2<T>> for &'a Vector2<T> {
            type Output = Vector2<T>;

            #[inline]
            fn $method(self, other: Vector2<T>) -> Self::Output {
                self.clone().$method(other)
            }
        }
    };
}

macro_rules! forward_val_xf_binop {
    (impl $imp:ident, $method:ident) => {
        impl<'a, T: Clone + Num> $imp<&'a Vector2<T>> for Vector2<T> {
            type Output = Vector2<T>;

            #[inline]
            fn $method(self, other: &Vector2<T>) -> Self::Output {
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
impl<T: Clone + Num> Add<Vector2<T>> for Vector2<T> {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self::Output::new(self.x_ + other.x_, self.y_ + other.y_)
    }
}

forward_all_binop!(impl Sub, sub);

// (a, b) - (c, d) == (a - c), (b - d)
impl<T: Clone + Num> Sub<Vector2<T>> for Vector2<T> {
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

    use crate::Vector2;

    impl<T: Clone + NumAssign> AddAssign for Vector2<T> {
        fn add_assign(&mut self, other: Self) {
            self.x_ += other.x_;
            self.y_ += other.y_;
        }
    }

    impl<T: Clone + NumAssign> SubAssign for Vector2<T> {
        fn sub_assign(&mut self, other: Self) {
            self.x_ -= other.x_;
            self.y_ -= other.y_;
        }
    }

    impl<T: Clone + NumAssign> MulAssign<T> for Vector2<T> {
        fn mul_assign(&mut self, other: T) {
            self.x_ *= other.clone();
            self.y_ *= other;
        }
    }

    impl<T: Clone + NumAssign> DivAssign<T> for Vector2<T> {
        fn div_assign(&mut self, other: T) {
            self.x_ /= other.clone();
            self.y_ /= other;
        }
    }

    macro_rules! forward_op_assign1 {
        (impl $imp:ident, $method:ident) => {
            impl<'a, T: Clone + NumAssign> $imp<&'a Vector2<T>> for Vector2<T> {
                #[inline]
                fn $method(&mut self, other: &Self) {
                    self.$method(other.clone())
                }
            }
        };
    }

    macro_rules! forward_op_assign2 {
        (impl $imp:ident, $method:ident) => {
            impl<'a, T: Clone + NumAssign> $imp<&'a T> for Vector2<T> {
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

impl<T: Clone + Num + Neg<Output = T>> Neg for Vector2<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self::Output::new(-self.x_, -self.y_)
    }
}

impl<T: Clone + Num + Neg<Output = T>> Neg for &Vector2<T> {
    type Output = Vector2<T>;

    #[inline]
    fn neg(self) -> Self::Output {
        -self.clone()
    }
}

macro_rules! scalar_arithmetic {
    (@forward $imp:ident::$method:ident for $($scalar:ident),*) => (
        impl<'a, T: Clone + Num> $imp<&'a T> for Vector2<T> {
            type Output = Vector2<T>;

            #[inline]
            fn $method(self, other: &T) -> Self::Output {
                self.$method(other.clone())
            }
        }
        impl<'a, T: Clone + Num> $imp<T> for &'a Vector2<T> {
            type Output = Vector2<T>;

            #[inline]
            fn $method(self, other: T) -> Self::Output {
                self.clone().$method(other)
            }
        }
        impl<'a, 'b, T: Clone + Num> $imp<&'a T> for &'b Vector2<T> {
            type Output = Vector2<T>;

            #[inline]
            fn $method(self, other: &T) -> Self::Output {
                self.clone().$method(other.clone())
            }
        }
        $(
            impl<'a> $imp<&'a Vector2<$scalar>> for $scalar {
                type Output = Vector2<$scalar>;

                #[inline]
                fn $method(self, other: &Vector2<$scalar>) -> Vector2<$scalar> {
                    self.$method(other.clone())
                }
            }
            impl<'a> $imp<Vector2<$scalar>> for &'a $scalar {
                type Output = Vector2<$scalar>;

                #[inline]
                fn $method(self, other: Vector2<$scalar>) -> Vector2<$scalar> {
                    self.clone().$method(other)
                }
            }
            impl<'a, 'b> $imp<&'a Vector2<$scalar>> for &'b $scalar {
                type Output = Vector2<$scalar>;

                #[inline]
                fn $method(self, other: &Vector2<$scalar>) -> Vector2<$scalar> {
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
            impl Mul<Vector2<$scalar>> for $scalar {
                type Output = Vector2<$scalar>;

                #[inline]
                fn mul(self, other: Vector2<$scalar>) -> Self::Output {
                    Self::Output::new(self * other.x_, self * other.y_)
                }
            }

        )*
    );
}

impl<T: Clone + Num> Mul<T> for Vector2<T> {
    type Output = Vector2<T>;

    #[inline]
    fn mul(self, other: T) -> Self::Output {
        Self::Output::new(self.x_ * other.clone(), self.y_ * other)
    }
}

impl<T: Clone + Num> Div<T> for Vector2<T> {
    type Output = Self;

    #[inline]
    fn div(self, other: T) -> Self::Output {
        Self::Output::new(self.x_ / other.clone(), self.y_ / other)
    }
}

impl<T: Clone + Num> Rem<T> for Vector2<T> {
    type Output = Vector2<T>;

    #[inline]
    fn rem(self, other: T) -> Self::Output {
        Self::Output::new(self.x_ % other.clone(), self.y_ % other)
    }
}

scalar_arithmetic!(usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64);

// constants
impl<T: Clone + Num> Zero for Vector2<T> {
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

    // use super::{hash, Vector2};
    use super::Vector2;
    use core::f64;
    use num_traits::Zero;

    pub const _0_0v: Vector2<f64> = Vector2 { x_: 0.0, y_: 0.0 };
    pub const _1_0v: Vector2<f64> = Vector2 { x_: 1.0, y_: 0.0 };
    pub const _1_1v: Vector2<f64> = Vector2 { x_: 1.0, y_: 1.0 };
    pub const _0_1v: Vector2<f64> = Vector2 { x_: 0.0, y_: 1.0 };
    pub const _neg1_1v: Vector2<f64> = Vector2 { x_: -1.0, y_: 1.0 };
    pub const _05_05v: Vector2<f64> = Vector2 { x_: 0.5, y_: 0.5 };
    pub const all_consts: [Vector2<f64>; 5] = [_0_0v, _1_0v, _1_1v, _neg1_1v, _05_05v];
    pub const _4_2v: Vector2<f64> = Vector2 { x_: 4.0, y_: 2.0 };

    #[test]
    fn test_consts() {
        // check our constants are what Vector2::new creates
        fn test(c: Vector2<f64>, r: f64, i: f64) {
            assert_eq!(c, Vector2::new(r, i));
        }
        test(_0_0v, 0.0, 0.0);
        test(_1_0v, 1.0, 0.0);
        test(_1_1v, 1.0, 1.0);
        test(_neg1_1v, -1.0, 1.0);
        test(_05_05v, 0.5, 0.5);
        assert_eq!(_0_0v, Zero::zero());
    }

    #[test]
    fn test_scale_unscale() {
        assert_eq!(_05_05v.scale(2.0), _1_1v);
        assert_eq!(_1_1v.unscale(2.0), _05_05v);
        for &c in all_consts.iter() {
            assert_eq!(c.scale(2.0).unscale(2.0), c);
        }
    }

    // #[test]
    // fn test_hash() {
    //     let a = Vector2::new(0i32, 0i32);
    //     let b = Vector2::new(1i32, 0i32);
    //     let c = Vector2::new(0i32, 1i32);
    //     assert!(hash(&a) != hash(&b));
    //     assert!(hash(&b) != hash(&c));
    //     assert!(hash(&c) != hash(&a));
    // }
}
