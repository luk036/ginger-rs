#![allow(dead_code)]

/// A reference-based 2D vector implementation that holds mutable references to f64 values.
/// 
/// This structure provides mathematical operations on 2D vectors using mutable references
/// to external f64 values rather than owning the values directly.
pub struct Vector2Ref<'a> {
    /// Mutable reference to the x component
    pub x: &'a mut f64,
    /// Mutable reference to the y component
    pub y: &'a mut f64,
}

impl<'a> Vector2Ref<'a> {
    /// Creates a new Vector2Ref from mutable references to two f64 values.
    /// 
    /// # Arguments
    /// 
    /// * `x` - A mutable reference to the x component
    /// * `y` - A mutable reference to the y component
    /// 
    /// # Examples
    /// 
    /// ```
    /// use ginger::vector2_ref::Vector2Ref;
    /// 
    /// let mut x_val = 1.0;
    /// let mut y_val = 2.0;
    /// 
    /// let mut v = Vector2Ref::new(&mut x_val, &mut y_val);
    /// ```
    pub fn new(x: &'a mut f64, y: &'a mut f64) -> Self {
        Vector2Ref { x, y }
    }

    /// Computes the dot product of this vector with another vector.
    /// 
    /// # Arguments
    /// 
    /// * `other` - Another Vector2Ref to compute the dot product with
    /// 
    /// # Examples
    /// 
    /// ```
    /// use ginger::vector2_ref::Vector2Ref;
    /// 
    /// let mut x1 = 1.0;
    /// let mut y1 = 2.0;
    /// let mut x2 = 3.0;
    /// let mut y2 = 4.0;
    /// 
    /// let v1 = Vector2Ref::new(&mut x1, &mut y1);
    /// let v2 = Vector2Ref::new(&mut x2, &mut y2);
    /// let result = v1.dot(&v2);
    /// 
    /// assert_eq!(result, 11.0); // 1*3 + 2*4 = 11
    /// ```
    pub fn dot(&self, other: &Vector2Ref) -> f64 {
        *self.x * *other.x + *self.y * *other.y
    }

    /// Computes the cross product of this vector with another vector.
    /// 
    /// # Arguments
    /// 
    /// * `other` - Another Vector2Ref to compute the cross product with
    /// 
    /// # Examples
    /// 
    /// ```
    /// use ginger::vector2_ref::Vector2Ref;
    /// 
    /// let mut x1 = 1.0;
    /// let mut y1 = 2.0;
    /// let mut x2 = 3.0;
    /// let mut y2 = 4.0;
    /// 
    /// let v1 = Vector2Ref::new(&mut x1, &mut y1);
    /// let v2 = Vector2Ref::new(&mut x2, &mut y2);
    /// let result = v1.cross(&v2);
    /// 
    /// assert_eq!(result, -2.0); // 1*4 - 3*2 = -2
    /// ```
    pub fn cross(&self, other: &Vector2Ref) -> f64 {
        *self.x * *other.y - *other.x * *self.y
    }

    /// Adds another vector to this vector in-place.
    /// 
    /// # Arguments
    /// 
    /// * `other` - Another Vector2Ref to add to this vector
    /// 
    /// # Examples
    /// 
    /// ```
    /// use ginger::vector2_ref::Vector2Ref;
    /// 
    /// let mut x1 = 1.0;
    /// let mut y1 = 2.0;
    /// let mut x2 = 3.0;
    /// let mut y2 = 4.0;
    /// 
    /// let mut v1 = Vector2Ref::new(&mut x1, &mut y1);
    /// let v2 = Vector2Ref::new(&mut x2, &mut y2);
    /// v1.add_assign(&v2);
    /// 
    /// assert_eq!(*v1.x, 4.0); // 1 + 3 = 4
    /// assert_eq!(*v1.y, 6.0); // 2 + 4 = 6
    /// ```
    pub fn add_assign(&mut self, other: &Vector2Ref) {
        *self.x += *other.x;
        *self.y += *other.y;
    }

    /// Subtracts another vector from this vector in-place.
    /// 
    /// # Arguments
    /// 
    /// * `other` - Another Vector2Ref to subtract from this vector
    /// 
    /// # Examples
    /// 
    /// ```
    /// use ginger::vector2_ref::Vector2Ref;
    /// 
    /// let mut x1 = 5.0;
    /// let mut y1 = 6.0;
    /// let mut x2 = 3.0;
    /// let mut y2 = 4.0;
    /// 
    /// let mut v1 = Vector2Ref::new(&mut x1, &mut y1);
    /// let v2 = Vector2Ref::new(&mut x2, &mut y2);
    /// v1.sub_assign(&v2);
    /// 
    /// assert_eq!(*v1.x, 2.0); // 5 - 3 = 2
    /// assert_eq!(*v1.y, 2.0); // 6 - 4 = 2
    /// ```
    pub fn sub_assign(&mut self, other: &Vector2Ref) {
        *self.x -= *other.x;
        *self.y -= *other.y;
    }

    /// Scales this vector by a scalar value in-place.
    /// 
    /// # Arguments
    /// 
    /// * `alpha` - The scalar value to multiply the vector by
    /// 
    /// # Examples
    /// 
    /// ```
    /// use ginger::vector2_ref::Vector2Ref;
    /// 
    /// let mut x = 2.0;
    /// let mut y = 3.0;
    /// 
    /// let mut v = Vector2Ref::new(&mut x, &mut y);
    /// v.mul_assign(2.0);
    /// 
    /// assert_eq!(*v.x, 4.0); // 2 * 2 = 4
    /// assert_eq!(*v.y, 6.0); // 3 * 2 = 6
    /// ```
    pub fn mul_assign(&mut self, alpha: f64) {
        *self.x *= alpha;
        *self.y *= alpha;
    }

    /// Divides this vector by a scalar value in-place.
    /// 
    /// # Arguments
    /// 
    /// * `alpha` - The scalar value to divide the vector by
    /// 
    /// # Examples
    /// 
    /// ```
    /// use ginger::vector2_ref::Vector2Ref;
    /// 
    /// let mut x = 6.0;
    /// let mut y = 8.0;
    /// 
    /// let mut v = Vector2Ref::new(&mut x, &mut y);
    /// v.div_assign(2.0);
    /// 
    /// assert_eq!(*v.x, 3.0); // 6 / 2 = 3
    /// assert_eq!(*v.y, 4.0); // 8 / 2 = 4
    /// ```
    pub fn div_assign(&mut self, alpha: f64) {
        *self.x /= alpha;
        *self.y /= alpha;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vector2() {
        let mut x = 1.0;
        let mut y = 2.0;

        let mut v = Vector2Ref::new(&mut x, &mut y);
        v.mul_assign(2.0);
        assert_eq!(*v.x, 2.0);
        assert_eq!(*v.y, 4.0);

        let mut v2 = Vector2Ref::new(&mut x, &mut y);
        v2.mul_assign(2.0);
        assert_eq!(*v2.y, 8.0);
    }
}
