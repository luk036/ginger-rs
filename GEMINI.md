# GEMINI.md

## Project Overview

This project, `ginger-rs`, is a Rust library and command-line application for finding the roots of polynomials. It implements the Bairstow and Aberth methods, including parallel versions for improved performance. The library is the primary artifact, with the command-line application serving as a simple entry point.

**Key Technologies:**

*   **Language:** Rust
*   **Core Libraries:**
    *   `num-complex`: For complex number arithmetic.
    *   `rayon`: For data parallelism.
    *   `lds-rs`: For low-discrepancy sequences.
*   **Testing & Benchmarking:**
    *   `criterion`: For benchmarking.
    *   `approx_eq`: For floating-point comparisons in tests.

**Architecture:**

The project is structured as a Rust library with a small binary wrapper. The core logic is in the `src` directory, with different modules for various algorithms and data structures:

*   `aberth.rs`: Implements the Aberth root-finding algorithm.
*   `horner.rs`: Implements Horner's method for polynomial evaluation.
*   `rootfinding.rs`: Implements the Bairstow root-finding algorithm.
*   `matrix2.rs`, `vector2.rs`, `vector2_ref.rs`: Provide simple matrix and vector implementations.

## Building and Running

### Building

To build the project, use the following command:

```bash
cargo build
```

### Running Tests

To run the test suite, use the following command:

```bash
cargo test
```

### Running Benchmarks

To run the benchmarks, use the following command:

```bash
cargo bench
```

## Development Conventions

*   **Coding Style:** The code follows standard Rust conventions.
*   **Testing:** Tests are included in the `src/lib.rs` file and can be run with `cargo test`.
*   **Contributions:** Contribution guidelines are outlined in `CONTRIBUTING.md`.
