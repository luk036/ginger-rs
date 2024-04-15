#![allow(dead_code)]

struct Vector2Ref<'a> {
    x: &'a mut f64,
    y: &'a mut f64,
}

impl<'a> Vector2Ref<'a> {
    fn new(x: &'a mut f64, y: &'a mut f64) -> Self {
        Vector2Ref { x, y }
    }

    fn dot(&self, other: &Vector2Ref) -> f64 {
        *self.x * *other.x + *self.y * *other.y
    }

    fn cross(&self, other: &Vector2Ref) -> f64 {
        *self.x * *other.y - *other.x * *self.y
    }

    fn add_assign(&mut self, other: &Vector2Ref) {
        *self.x += *other.x;
        *self.y += *other.y;
    }

    fn sub_assign(&mut self, other: &Vector2Ref) {
        *self.x -= *other.x;
        *self.y -= *other.y;
    }

    fn mul_assign(&mut self, alpha: f64) {
        *self.x *= alpha;
        *self.y *= alpha;
    }

    fn div_assign(&mut self, alpha: f64) {
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
