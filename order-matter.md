# Order-Independence in Parallel Simultaneous Root-Finding Methods — A Corrected Analysis 🧮

## 1. Introduction & Core Thesis 🎯

In numerical polynomial root-finding, simultaneous methods (which compute all roots at once) offer significant advantages for parallel computing. A key theoretical question is whether the **order** in which individual root approximations are updated affects the final result.

This report investigates two simultaneous algorithms:

- **The Aberth–Ehrlich (AE) method**: Operates directly on complex roots.
- **The Simultaneous Bairstow (PAID) method**: Operates on quadratic factors.

**⚠️ Important correction to a prior draft**: order-independence holds **only for the Jacobi-style (parallel) iteration**, and it holds in a *stronger* sense than previously claimed. It does **not** hold for Gauss–Seidel-style (sequential) iteration. This report:

1. Explains *why* Jacobi is order-independent (bit-for-bit, not merely "up to floating point").
2. Identifies a **flaw in the earlier experiment** (it tested polynomial *multiplication*, not the algorithm's actual operation).
3. Replaces it with experiments on the **real operation** — the sequential 2×2 linear solves (`suppress_old`).
4. Demonstrates that the **single-threaded variants in the library are Gauss–Seidel and are order-dependent**.

---

## 2. Theoretical Background 📚

### 2.1 The Aberth–Ehrlich Method

The Aberth–Ehrlich iteration for finding roots \( z_i \) of a polynomial \( P(z) \):

$$
z_i^{(new)} = z_i - \frac{P(z_i)}{\prod_{j \neq i} (z_i - z_j)} \cdot \frac{1}{1 - \frac{P(z_i)}{\prod_{j \neq i} (z_i - z_j)} \sum_{j \neq i} \frac{1}{z_i - z_j}}
$$

- **Jacobi (Parallel)**: All \( z_j \) on the RHS come from the *previous* iteration. Each update is a pure function of the frozen state → processing order is irrelevant, **bit-for-bit**.
- **Gauss–Seidel (Sequential)**: Updates use freshly-computed values immediately → highly order-dependent.

> 💡 **Key Insight**: For Jacobi, order-independence is *structural*: every job reads only the frozen previous state and writes only its own slot. No commutativity argument is even needed.

---

### 2.2 The Simultaneous Bairstow Method & PAID

Bairstow's method finds quadratic factors:

$$
F_i(x) = x^2 + b_i x + c_i \qquad \text{(in ginger-rs: } x^2 - r_i x - q_i)
$$

The **Parallel Anticipatory Implicit Deflation (PAID)** strategy updates all factors simultaneously. For factor \( i \), the deflated polynomial is:

$$
D_i(x) = \prod_{j \neq i} F_j^{(old)}(x)
$$

and a local \( 2 \times 2 \) system gives the correction:

$$
J_i
\begin{bmatrix} \Delta b_i \\ \Delta c_i \end{bmatrix}
= -D_i(x)
$$

#### The Actual Suppression Operation 🔧

The library does **not** multiply polynomials to obtain \( D_i(x) \). Instead, each job starts from the Horner remainder of \( P \) modulo \( F_i \) and then **composes sequential 2×2 linear solves** (`suppress_old`), one per other factor:

$$
\begin{bmatrix} r p + s & p \\ q p & s \end{bmatrix}
\begin{bmatrix} a \\ b \end{bmatrix}
=
\begin{bmatrix} A \\ B \end{bmatrix}
\qquad
(p, s) = (r_i - r_j,\; q_i - q_j)
$$

Each solve **mutates the remainder in place** (`src/rootfinding.rs:488-490`). So the earlier draft's justification — "D_i(x) is a product of scalars, which is commutative" — **does not describe the code**. The real question is whether *sequential linear solves* commute. The answer turns out to be **yes, exactly** (Section 4).

---

## 3. Jacobi vs Gauss–Seidel in the Library ⚙️

The ginger-rs library contains **both** variants, and they have opposite order-dependence properties:

| Function | Scheme | Data flow |
|----------|--------|-----------|
| `pbairstow_even_mt` / `pbairstow_autocorr_mt` | **Jacobi** 🟢 | `vrsc.copy_from_slice(vrs)` freezes the state (`:435`, `:621`); jobs read only the snapshot, write only their own slot |
| `pbairstow_even` / `pbairstow_autocorr` | **Gauss–Seidel** 🔴 | Live `vrs` passed to each job (`:386`, `:568`); written back in place (`:391`, `:574`) |

The code comment itself acknowledges the second variant: `// Gauss-Seidel fashion` (`:491`, `:681`).

> 🔑 **Consequence**: "Implementing the method with a for loop `0,1,2` or `2,1,0` yields identical results" is **true for `_mt` (Jacobi)** but **false for the single-threaded functions** (Gauss–Seidel).

---

## 4. Experiments 🧪

The prior draft's experiment tested **polynomial multiplication** (`numpy.polynomial.Polynomial` products) — an operation the algorithm never performs. It trivially showed \( 10^{-14} \) differences and proved nothing about the algorithm. Below are corrected experiments on the **actual operation**, all reproducible from `examples/order_experiment.rs` and `examples/trace_rust.rs` in ginger-rs.

### 4.1 E1 — Suppression order in **exact** arithmetic (the decisive test)

Faithful mirror of `suppress_old`/`horner`/`delta` using Python `fractions.Fraction` (no rounding). For a degree-8 palindromic test polynomial and 4 factors, all \( 4! = 24 \) suppression orders were compared:

```
factor i=0: 24 suppression orders -> ALL IDENTICAL (order-independent)
factor i=1: 24 suppression orders -> ALL IDENTICAL (order-independent)
factor i=2: 24 suppression orders -> ALL IDENTICAL (order-independent)
factor i=3: 24 suppression orders -> ALL IDENTICAL (order-independent)
```

> ✅ **The suppression composition is *exactly* commutative in exact arithmetic.** Each 2×2 solve is an exact deflation step; deflating by the product of the other factors is order-independent. The prior draft's *conclusion* survives — for a sounder reason than it gave.

### 4.2 E2 — Suppression order in floating point

Same experiment in f64:

```
worst |diff| in updated vri across all suppression orders: 1.066e-14
```

> ✅ Differences are **machine-epsilon level** (\( \sim 10^{-14} \)). Even with a **near-singular** suppression system (two factors differing by \( 10^{-8} \), determinant \( \sim 10^{-16} \)), the differences stayed **exactly zero**.

### 4.3 E3 — Factor *processing* order in the Jacobi sweep

Process the factor jobs in orders `[0,1,2,3]`, `[3,2,1,0]`, and two shuffles:

```
bit-identical next state for all processing orders: True
bit-identical for 20 consecutive sweeps across orders: True
```

> ✅ **Bit-for-bit identical.** This is the strongest form of order-independence — not "close to equal", **exactly equal**, because each job is a pure function of the frozen snapshot.

### 4.4 E4 — Factor *processing* order in the Gauss–Seidel sweep

Same processing orders, but on the live-array (Gauss–Seidel) sweep:

```
bit-identical next state for all processing orders: False
worst |diff| between natural and other orders: 9.819e+00
```

> 🔴 **Huge divergence — differences of order \( 10 \).** Processing order fundamentally changes the Gauss–Seidel trajectory.

### 4.5 E5 — Full convergence with permuted initial guesses

Run to convergence starting from natural / reversed / shuffled initial guesses:

| Variant | natural | reversed | shuffle-1 | shuffle-2 |
|---------|---------|----------|-----------|-----------|
| **Jacobi** (`pbairstow_even_mt`) | niter = 12 | niter = 12 | niter = 12 | niter = 12 |
| **Gauss–Seidel** (`pbairstow_even`) | niter = 10 | niter = 7 | niter = 9 | niter = 6 |

- **Jacobi**: identical iteration count *and* identical root set for every permutation.
- **Gauss–Seidel**: iteration count varies **6–10** depending on order (roots converge to the same set for this well-conditioned polynomial, but convergence speed differs wildly).

---

## 5. Discussion: What This Means for the Full Algorithm 💡

1. **Suppression order within a job** — *irrelevant*: exactly commutative in exact arithmetic (E1), \( \sim 10^{-14} \) in f64 (E2).
2. **Job processing order in Jacobi** — *irrelevant*: bit-identical next state (E3). Frozen snapshot + disjoint writes ⇒ pure functions.
3. **Job processing order in Gauss–Seidel** — *matters*: next-state differences of order \( 10 \) (E4), iteration counts 6–10 (E5).
4. **Which function you call decides everything**: the `_mt` variants are Jacobi; the single-threaded variants are Gauss–Seidel despite the parallel-sounding name.

---

## 6. Conclusion ✅

Based on theoretical analysis and the corrected experiments:

- **For Aberth–Ehrlich (Jacobi)**: strictly order-independent — each update is a pure function of the frozen state.
- **For Simultaneous Bairstow (Jacobi / `_mt`)**:
  - The **order of suppression** is mathematically irrelevant: the sequential 2×2 solves are **exactly commutative** in exact arithmetic and differ only at the \( 10^{-14} \) level in floating point.
  - The **outer processing order** is irrelevant **bit-for-bit**: jobs read only the frozen snapshot and write disjoint slots.
- **For the single-threaded variants (`pbairstow_even`, `pbairstow_autocorr`)**:
  - These are **Gauss–Seidel**, and are **order-dependent**: next-state differences of order \( 10 \), iteration counts varying 6–10 with order.

> 🚀 **Final Verdict**: The order-independence claim is **true for the parallel (Jacobi) variants** — stronger than the original draft claimed (bit-for-bit, not just \( 10^{-14} \)). But the original draft **overreached**: it claimed the property for "either method". The single-threaded Gauss–Seidel variants are genuinely order-dependent. Use the `_mt` entry points for deterministic, order-free parallelism.

---

## Appendix: Reproducing the Experiments 🔬

- **Rust**: `cargo run --example order_experiment` — permutes initial guesses through the real library functions (`pbairstow_even` vs `pbairstow_even_mt`).
- **Rust trace**: `cargo run --example trace_rust` — per-iteration trajectory comparison.
- **Python mirror**: faithful re-implementation of `horner` / `suppress_old` / `delta` (exact `Fraction` and f64 modes) used for E1–E5.
- **Test polynomial**: degree-8 palindromic `[10, 34, 75, 94, 150, 94, 75, 34, 10]` (used throughout the library's test suite).
