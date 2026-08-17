//! Seqlock-based atomic pairs for the atomic multi-threading method.
//!
//! `AtomicU128` is unstable on stable Rust, and two separate `AtomicU64`s
//! allow torn (x, y) reads. A seqlock — a version counter plus two bit-pattern
//! `AtomicU64`s — gives torn-free whole-pair reads without 128-bit atomics.
//!
//! Each slot is single-writer, multi-reader: the owning thread calls `store`
//! while every other thread only calls `load`. The writer bumps the sequence
//! to odd while updating, then to even when stable; readers retry until they
//! observe a stable (even) sequence with matching reads of both components.

use crate::vector2::Vector2;
use num_complex::Complex;
use std::sync::atomic::{AtomicU64, Ordering};

/// Seqlock over a pair of `f64` values (single-writer, multi-reader).
///
/// Writer protocol: `seq` becomes odd while `c0`/`c1` are being updated, then
/// even again once the pair is stable. Readers spin until they read a stable
/// even sequence number, guaranteeing the two components are from one write.
pub struct AtomicF64x2 {
    seq: AtomicU64,
    c0: AtomicU64,
    c1: AtomicU64,
}

impl AtomicF64x2 {
    /// Create a seqlock initialized with the given pair.
    pub fn new(x: f64, y: f64) -> Self {
        AtomicF64x2 {
            seq: AtomicU64::new(0),
            c0: AtomicU64::new(x.to_bits()),
            c1: AtomicU64::new(y.to_bits()),
        }
    }

    /// Torn-free read of the pair; retries while the writer is mid-update.
    pub fn load(&self) -> [f64; 2] {
        loop {
            let s1 = self.seq.load(Ordering::Acquire);
            let c0 = f64::from_bits(self.c0.load(Ordering::Relaxed));
            let c1 = f64::from_bits(self.c1.load(Ordering::Relaxed));
            let s2 = self.seq.load(Ordering::Acquire);
            if s1 == s2 && s1 % 2 == 0 {
                return [c0, c1];
            }
        }
    }

    /// Store a new pair; must be called by the single writer of this slot.
    pub fn store(&self, x: f64, y: f64) {
        let s = self.seq.load(Ordering::Relaxed);
        self.seq.store(s + 1, Ordering::Release); // odd: writing
        self.c0.store(x.to_bits(), Ordering::Relaxed);
        self.c1.store(y.to_bits(), Ordering::Relaxed);
        self.seq.store(s + 2, Ordering::Release); // even: stable
    }
}

/// Atomic `Complex<f64>` (seqlock over re/im) — whole-pair atomicity.
pub struct AtomicComplex(AtomicF64x2);

impl AtomicComplex {
    /// Create an atomic complex value.
    pub fn new(z: Complex<f64>) -> Self {
        AtomicComplex(AtomicF64x2::new(z.re, z.im))
    }

    /// Torn-free read of the complex value.
    pub fn load(&self) -> Complex<f64> {
        let [re, im] = self.0.load();
        Complex::new(re, im)
    }

    /// Store a new complex value; single-writer per slot.
    pub fn store(&self, z: Complex<f64>) {
        self.0.store(z.re, z.im);
    }
}

/// Atomic `Vector2<f64>` (seqlock over x/y) — whole-pair atomicity.
pub struct AtomicVec2(AtomicF64x2);

impl AtomicVec2 {
    /// Create an atomic vector value.
    pub fn new(v: Vector2<f64>) -> Self {
        AtomicVec2(AtomicF64x2::new(v.x_, v.y_))
    }

    /// Torn-free read of the vector value.
    pub fn load(&self) -> Vector2<f64> {
        let [x, y] = self.0.load();
        Vector2::new(x, y)
    }

    /// Store a new vector value; single-writer per slot.
    pub fn store(&self, v: Vector2<f64>) {
        self.0.store(v.x_, v.y_);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atomic_f64x2_roundtrip() {
        let slot = AtomicF64x2::new(1.0, 2.0);
        assert_eq!(slot.load(), [1.0, 2.0]);
        slot.store(3.0, 4.0);
        assert_eq!(slot.load(), [3.0, 4.0]);
    }

    #[test]
    fn test_atomic_complex_roundtrip() {
        let slot = AtomicComplex::new(Complex::new(1.0, -2.5));
        assert_eq!(slot.load(), Complex::new(1.0, -2.5));
        slot.store(Complex::new(-3.0, 0.5));
        assert_eq!(slot.load(), Complex::new(-3.0, 0.5));
    }

    #[test]
    fn test_atomic_vec2_roundtrip() {
        let slot = AtomicVec2::new(Vector2::new(1.0, 2.0));
        assert_eq!(slot.load(), Vector2::new(1.0, 2.0));
        slot.store(Vector2::new(3.0, 4.0));
        assert_eq!(slot.load(), Vector2::new(3.0, 4.0));
    }
}
