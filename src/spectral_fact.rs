extern crate ndarray;
extern crate num_complex;

use ndarray::{arr1, Array1, Axis};
use num_complex::Complex;
use std::f64::consts::PI;

pub fn spectral_fact(r: &Array1<f64>) -> Array1<f64> {
    let n = r.len();
    let mult_factor = 100; 
    let m = mult_factor * n;

    let w = Array1::range(0., 2. * PI, (m as f64) / (m as f64 - 1.));
    let wn = Array1::from_iter((1..n).map(|k| w.clone().mapv(|w_i| 2. * w_i * k as f64)));
    let cos_wn = wn.mapv(f64::cos);

    let mut R = cos_wn.to_shape((m, n - 1)).unwrap().t().to_owned();
    R.insert_axis_inplace(Axis(1));
    R.slice_mut(s![.., 0]).assign(&Array1::from_elem(m, 1.));
    R *= r;

    let alpha = 0.5 * R.mapv(f64::ln).mapv(f64::abs);

    let fft_alpha = ndarrayfft::fft(&alpha.view(), false);
    let mut alphatmp = fft_alpha.to_owned();
    for i in (n/2)..m { alphatmp[i] = -alphatmp[i]; }
    alphatmp[0] = 0.0;

    let phi = ndarrayfft::ifft(&Complex::new(0., 1.) * alphatmp.view()).mapv(|c| c.re);

    let index = (0..m).step_by(mult_factor).collect::<Vec<usize>>();
    let alpha1 = alpha.select(Axis(0), &index.into());
    let phi1 = phi.select(Axis(0), &index.into());

    let h = ndarrayfft::ifft(&alpha1.mapv(Complex::from) + Complex::new(0., 1.) * phi1.mapv(Complex::from), false)
                     .mapv(|c| c.re)
                     .into_shape(n)
                     .unwrap();

    h
}

pub fn inverse_spectral_fact(h: &Array1<f64>) -> Array1<f64> {
    let mut r = Array1::zeros(h.len());
    for t in 0..h.len() {
        let slice = if t < h.len() - t { &h[t..] } else { &h[h.len() - t..] };
        r[t] = slice.dot(&h.slice(s![..h.len()-t]));
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::*;

    #[test]
    fn test_spectral_fact() {
        let h = arr1(&[
            0.76006445,
            0.54101887,
            0.42012073,
            0.3157191,
            0.10665804,
            0.04326203,
            0.01315678,
        ]);
        let r = inverse_spectral_fact(&h);
        let h2 = spectral_fact(&r);
        assert_eq!(h.len(), h2.len());
        assert_abs_diff_eq!(h, h2, epsilon = 1e-10);
    }
}
