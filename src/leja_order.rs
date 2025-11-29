use num_complex::Complex;
use std::cmp::Ordering;

/// Computes the Leja order of a set of complex points.
///
/// The Leja order is a sequence of points that maximizes the minimum distance between
/// consecutive points. This is useful for numerical methods that require a well-conditioned
/// set of points, such as polynomial interpolation.
///
/// # Arguments
/// * `points` - A vector of complex points to be ordered.
///
/// # Returns
/// A vector of complex points in Leja order.
///
/// # Panics
/// Panics if the input vector is empty.
///
/// # Examples
///
/// ```
/// use ginger::leja_order::leja_order;
/// use num_complex::Complex;
///
/// let points = vec![
///     Complex::new(1.0, 1.0),
///     Complex::new(2.0, -2.0),
///     Complex::new(0.5, 0.5),
///     Complex::new(-1.0, 0.0),
///     Complex::new(3.0, 3.0),
/// ];
/// let leja_ordered = leja_order(points);
/// assert_eq!(leja_ordered.len(), 5);
/// ```
pub fn leja_order(mut points: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
    // Check if input is empty and return an error if so
    if points.is_empty() {
        panic!("Input must be a non-empty vector of points.");
    }

    // Start with the point having the smallest magnitude
    let mut leja_ordered_points = vec![points.remove(
        points
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.norm().partial_cmp(&b.norm()).unwrap_or(Ordering::Equal))
            .map(|(idx, _)| idx)
            .unwrap(),
    )];

    while !points.is_empty() {
        // Compute distances from remaining points to the last point in leja_order
        let last_point = leja_ordered_points.last().unwrap();
        let distances = points.iter().map(|&p| (p - *last_point).norm());

        // Find the index of the point with the maximum minimum distance
        let next_idx = distances
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(Ordering::Equal))
            .map(|(idx, _)| idx)
            .unwrap();

        // Append this point to the leja_ordered_points
        leja_ordered_points.push(points.remove(next_idx));
    }

    leja_ordered_points
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let points = vec![
            Complex::new(1.0, 1.0),
            Complex::new(2.0, -2.0),
            Complex::new(0.5, 0.5),
            Complex::new(-1.0, 0.0),
            Complex::new(3.0, 3.0),
        ];
        let leja_ordered = leja_order(points);
        println!("{:?}", leja_ordered);
        assert_eq!(leja_ordered[4], Complex::new(1.0, 1.0));
    }
}
