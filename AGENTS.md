# AGENTS.md - Agent Coding Guide for ginger-rs

## Overview

ginger-rs is a Rust library for parallel polynomial root-finding using Bairstow's method and Aberth's method. It uses the Rust 2021 edition and supports parallel execution via Rayon.

---

## Build, Test, and Lint Commands

### Build
```bash
# Debug build
cargo build

# Release build
cargo build --release

# Run the binary
cargo run
cargo run --release
```

### Test
```bash
# Run all tests
cargo test

# Run all tests with all features
cargo test --all-features --workspace

# Run a specific test
cargo test test_name

# Run tests with output
cargo test -- --nocapture
cargo test --verbose

# Run tests with 4 threads (for parallel tests)
cargo test -- --test-threads 4
```

### Lint and Format
```bash
# Run Clippy lints
cargo clippy

# Run Clippy with all features
cargo clippy --all-targets --all-features --workspace

# Check formatting
cargo fmt --all -- --check

# Format code
cargo fmt --all
```

### Documentation
```bash
# Build docs (fails on warnings)
cargo doc --no-deps --document-private-items --all-features --workspace --examples

# With RUSTDOCFLAGS for warnings as errors
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --document-private-items
```

### Benchmarks
```bash
# Run criterion benchmarks
cargo bench
```

---

## Code Style Guidelines

### General Conventions

- **Edition**: Rust 2021 (`edition = "2021"` in Cargo.toml)
- **Dependencies**: External crates in Cargo.toml (`num-complex`, `num-traits`, `rayon`, `lds-rs`)
- **Dev dependencies**: `approx_eq`, `criterion`

### Imports

- Use `use` statements for traits at the top of modules
- Use external crates directly (e.g., `num_complex::Complex`, `num_traits::Num`, `rayon::prelude::*`)
- For internal modules, use relative paths: `use super::module_name;` or `use crate::module_name;`

```rust
use core::ops::{Add, Div, Mul, Neg, Rem, Sub};
use num_traits::{Num, Signed, Zero};
use lds_rs::lds::Circle;
use num_complex::Complex;
```

### Naming Conventions

- **Structs/Types**: `PascalCase` (e.g., `Vector2`, `Matrix2`, `Options`)
- **Functions/Methods**: `snake_case` (e.g., `initial_guess`, `pbairstow_even`)
- **Fields**: `snake_case` with underscore suffix (e.g., `x_`, `y_`, `max_iters`, `tolerance`)
- **Constants**: `SCREAMING_SNAKE_CASE` or single letter for math constants (e.g., `PI`, `TWO_PI`)

### Type Annotations

- Use explicit type annotations in function signatures for clarity
- Use generics with trait bounds (e.g., `impl<T: Clone + Num> Vector2<T>`)
- Use `type Vec2 = Vector2<f64>;` for common type aliases

```rust
pub struct Vector2<T> {
    pub x_: T,
    pub y_: T,
}

pub struct Options {
    pub max_iters: usize,
    pub tolerance: f64,
    pub tol_ind: f64,
}

type Vec2 = Vector2<f64>;
type Mat2 = Matrix2<f64>;
```

### Formatting

- Follow default `rustfmt` rules (4-space indentation)
- Use `#[inline]` attribute for small, frequently called functions
- Use `#[allow(non_snake_case)]` attribute at module top for math/algorithm code where notation uses non snake_case

```rust
#![allow(non_snake_case)]

#[inline]
pub fn new(x_: T, y_: T) -> Self {
    Vector2 { x_, y_ }
}
```

### Documentation

- Use doc comments (`///`) for public API
- Include Examples section in docs with ````rust` code blocks
- Document argument names, returns, and examples

```rust
/// Creates a new [`Vector2<T>`].
///
/// # Arguments
///
/// * `x_`: The x-coordinate of the vector.
/// * `y_`: The y-coordinate of the vector.
///
/// # Returns
///
/// A new `Vector2<T>` instance.
///
/// # Examples
///
/// ```
/// use ginger::vector2::Vector2;
///
/// let v = Vector2::new(3, 4);
/// assert_eq!(v.x_, 3);
/// ```
#[inline]
pub const fn new(x_: T, y_: T) -> Self {
    Vector2 { x_, y_ }
}
```

### Error Handling

- Use `Result<T, E>` for fallible operations when appropriate
- The library uses tolerance-based convergence checking (not traditional error handling)
- Return `(usize, bool)` tuples: `(number_of_iterations, converged)`

```rust
pub fn pbairstow_even(coeffs: &[f64], vrs: &mut [Vec2], options: &Options) -> (usize, bool) {
    for niter in 1..options.max_iters {
        // ... algorithm logic ...
        if tolerance < options.tolerance {
            return (niter, true);
        }
    }
    (options.max_iters, false)
}
```

### Testing

- Place tests in `#[cfg(test)]` modules at the bottom of source files
- Use `#[test]` attribute for test functions
- Use `approx_eq::assert_approx_eq!` for floating-point comparisons
- Use `super::*` to import tested items

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use approx_eq::assert_approx_eq;

    #[test]
    fn test_name() {
        // test code
    }
}
```

### Performance

- Use `rayon` for parallelization when appropriate
- Use `par_iter_mut()` for parallel iteration
- Use `#[inline]` on small hot functions
- Consider cache locality and memory layout

```rust
use rayon::prelude::*;

pub fn pbairstow_even_mt(coeffs: &[f64], vrs: &mut Vec<Vec2>, options: &Options) -> (usize, bool) {
    // ... setup ...
    let tol_i = vrs
        .par_iter_mut()
        .enumerate()
        .map(|(i, vri)| algorithm_job(coeffs, i, vri))
        .reduce(|| 0.0, |x, y| x.max(y));
    // ...
}
```

### Private Implementation Details

- Prefix with `_` for intentionally unused parameters: `|_idx| { ... }`
- Use internal (private) functions for job/worker functions
- Helper functions often defined within function bodies when localized

---

## Project Structure

```
ginger-rs/
├── src/
│   ├── lib.rs          # Library root, exports
│   ├── main.rs        # Binary entry point
│   ├── vector2.rs     # 2D vector implementation
│   ├── matrix2.rs    # 2x2 matrix implementation
│   ├── horner.rs      # Horner's method for polynomial eval
│   ├── aberth.rs     # Aberth's root-finding method
│   ├── rootfinding.rs # Bairstow's method
│   ├── leja_order.rs # Leja ordering
│   └── ...
├── Cargo.toml
├── CONTRIBUTING.md
└── README.md
```

---

## Common Patterns

### Options/Configuration Struct

```rust
#[derive(Debug)]
pub struct Options {
    pub max_iters: usize,
    pub tolerance: f64,
    pub tol_ind: f64,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            max_iters: 2000,
            tolerance: 1e-12,
            tol_ind: 1e-15,
        }
    }
}
```

### Initial Guess Functions

```rust
pub fn initial_guess(coeffs: &[f64]) -> Vec<Vec2> {
    let degree = coeffs.len() - 1;
    // ... compute initial guesses ...
}
```

### Main Algorithm Pattern

```rust
pub fn algorithm_name(coeffs: &[f64], roots: &mut [Type], options: &Options) -> (usize, bool) {
    for niter in 0..options.max_iters {
        let tolerance = compute_iteration(roots);
        if tolerance < options.tolerance {
            return (niter, true);
        }
    }
    (options.max_iters, false)
}
```

---

## CI/Testing Workflow

The project tests against:
- Rust stable, beta, and nightly
- Runs `cargo test`, `cargo clippy`, `cargo doc` in CI

Before submitting:
1. Run `cargo fmt --all`
2. Run `cargo clippy --all-targets --all-features --workspace`
3. Run `cargo test --all-features --workspace`
4. Build docs: `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --document-private-items`