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
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --document-private-items --all-features --workspace --examples
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
- **Dependencies**: `num-complex`, `num-traits`, `rayon`, `lds-rs`
- **Dev dependencies**: `approx_eq`, `criterion`

### Imports

- Use `use` statements for traits at the top of modules
- Use external crates directly (e.g., `num_complex::Complex`, `num_traits::Num`)
- For internal modules, use relative paths: `use super::module_name;` or `use crate::module_name;`

```rust
use core::ops::{Add, Div, Mul, Neg, Rem, Sub};
use num_traits::{Num, Signed, Zero};
```

### Naming Conventions

- **Structs/Types**: `PascalCase` (e.g., `Vector2`, `Matrix2`, `Options`)
- **Functions/Methods**: `snake_case` (e.g., `initial_guess`, `pbairstow_even`)
- **Fields**: `snake_case` with underscore suffix (e.g., `x_`, `y_`, `max_iters`, `tolerance`)

### Formatting

- Follow default `rustfmt` rules (4-space indentation)
- Use `#[inline]` attribute for small, frequently called functions
- Use `#[allow(non_snake_case)]` for math/algorithm code where notation uses non-snake_case

### Documentation

- Use doc comments (`///`) for public API
- Include Examples section in docs with ````rust` code blocks

### Error Handling

- Return `(usize, bool)` tuples: `(number_of_iterations, converged)`
- Use tolerance-based convergence checking

### Testing

- Place tests in `#[cfg(test)]` modules at the bottom of source files
- Use `#[test]` attribute for test functions
- Use `approx_eq::assert_approx_eq!` for floating-point comparisons

### Performance

- Use `rayon` for parallelization
- Use `par_iter_mut()` for parallel iteration
- Use `#[inline]` on small hot functions

---

## Project Structure

```
ginger-rs/
├── src/
│   ├── lib.rs          # Library root, exports
│   ├── main.rs         # Binary entry point
│   ├── vector2.rs      # 2D vector implementation
│   ├── matrix2.rs    # 2x2 matrix implementation
│   ├── horner.rs      # Horner's method for polynomial eval
│   ├── aberth.rs      # Aberth's root-finding method
│   ├── rootfinding.rs # Bairstow's method
│   └── ...
├── Cargo.toml
└── README.md
```

---

## CI/Testing

The project tests against:
- Rust stable, beta, and nightly
- Runs `cargo test`, `cargo clippy`, `cargo doc`, `cargo fmt` in CI

Before submitting:
1. Run `cargo fmt --all`
2. Run `cargo clippy --all-targets --all-features --workspace`
3. Run `cargo test --all-features --workspace`