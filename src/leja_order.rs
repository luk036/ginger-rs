use num_complex::Complex;

/// Computes the Leja order of a set of complex points.
///
/// Starts with $$ p_0 = \arg\min |p| $$, then iteratively selects:
///
/// $$ p_k = \arg\max_{p \in \text{remaining}} \min_{q \in \text{selected}} |p - q| $$
///
/// This greedy algorithm maximizes the minimum distance between successive points,
/// reducing numerical error when reconstructing polynomials from roots.
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

    // Greedy Leja ordering. Maintaining the running minimum distance to the
    // selected set turns the naive O(n^3) rescan into O(n^2).
    points.sort_by(|a, b| a.norm().partial_cmp(&b.norm()).unwrap());

    let mut result = Vec::with_capacity(points.len());
    result.push(points[0]);
    points.remove(0);

    let mut dists: Vec<f64> = points.iter().map(|p| (*p - result[0]).norm()).collect();

    while !points.is_empty() {
        let mut best_idx = 0;
        let mut best_dist = dists[0];
        for (i, d) in dists.iter().enumerate().skip(1) {
            if *d > best_dist {
                best_dist = *d;
                best_idx = i;
            }
        }
        let chosen = points[best_idx];
        points.remove(best_idx);
        dists.remove(best_idx);
        result.push(chosen);
        for (d, p) in dists.iter_mut().zip(points.iter()) {
            let nd = (*p - chosen).norm();
            if nd < *d {
                *d = nd;
            }
        }
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

    /// Reference O(n^3) greedy Leja ordering used to validate `leja_order`.
    fn naive_leja(mut points: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
        if points.is_empty() {
            return vec![];
        }
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

    fn lcg_points(n: usize, seed: u64) -> Vec<Complex<f64>> {
        let mut state = seed;
        let mut next = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((state >> 11) as f64) / ((1u64 << 53) as f64)
        };
        (0..n).map(|_| Complex::new(next(), next())).collect()
    }

    #[test]
    fn test_matches_naive_reference_exactly() {
        for n in [2usize, 3, 10, 37, 100] {
            let points = lcg_points(n, 7);
            assert_eq!(leja_order(points.clone()), naive_leja(points), "n={n}");
        }
    }

    #[test]
    fn test_returns_permutation_of_input() {
        let points = lcg_points(64, 11);
        let mut ordered = leja_order(points.clone());
        let mut sorted_input = points;
        let key = |z: &Complex<f64>| (z.re, z.im);
        ordered.sort_by(|a, b| key(a).partial_cmp(&key(b)).unwrap());
        sorted_input.sort_by(|a, b| key(a).partial_cmp(&key(b)).unwrap());
        assert_eq!(ordered, sorted_input);
    }
}
