use num_complex::Complex;

/// Computes the Leja order of a set of complex points.
///
/// Reorders complex points using the greedy Leja algorithm: starts with the
/// smallest-magnitude point, then iteratively selects the remaining point that
/// maximizes the minimum Euclidean distance to all already-selected points.
/// This ordering reduces numerical error when reconstructing polynomials from roots.
///
/// # Arguments
/// * `points` - Input vector of complex numbers
///
/// # Returns
/// A vector of complex points in Leja sequence. Returns an empty vector if input is empty.
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
    if points.is_empty() {
        return vec![];
    }

    // Greedy Leja ordering (O(n^2)):
    // 1. Start with the smallest-magnitude point
    // 2. Each subsequent point maximizes the MIN distance to ALL already-selected points
    points.sort_by(|a, b| a.norm().partial_cmp(&b.norm()).unwrap());

    let mut result = Vec::with_capacity(points.len());
    result.push(points[0]);
    points.remove(0);

    while !points.is_empty() {
        let mut best_idx = 0;
        let mut best_dist = -1.0f64;
        for (i, pi) in points.iter().enumerate() {
            let mut min_dist = f64::MAX;
            for p in &result {
                let d = (*pi - p).norm();
                if d < min_dist {
                    min_dist = d;
                }
            }
            if min_dist > best_dist {
                best_dist = min_dist;
                best_idx = i;
            }
        }
        result.push(points[best_idx]);
        points.remove(best_idx);
    }

    result
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

    #[test]
    fn test_empty_input() {
        let result = leja_order(vec![]);
        assert!(result.is_empty());
    }

    #[test]
    fn test_single_point() {
        let result = leja_order(vec![Complex::new(1.0, 2.0)]);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], Complex::new(1.0, 2.0));
    }
}
