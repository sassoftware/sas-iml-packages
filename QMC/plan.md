# SAS IML QMC Implementation Plan

## User Prompt: Implementation Planning

I want to start a new project to implement a subsystem in SAS IML that will
generate quasi-random sequences in $k$ dimensions in $[0,1]^k$. I want to provide
sufficient functionality so that customers can obtain QMC and randomized QMC
estimates of integrals on $[0,1]^k$.

I have seen implementations in MATLAB, the qmcpy package in Python, the
randtoolbox package in R, and the cooltools package in R. (I do not know if these
are good implementations. They are merely some implementations that I've seen.)
I have read several books, but my favorites are Gentle's *Random Number
Generation and Monte Carlo Methods* and Lemieux's *Monte Carlo and Quasi-Monte
Carlo Sampling*.

The fundamental pieces seem to be lattices, digital nets and sequences
(including Sobol, Faure, Niederreiter, Halton, ...). I want randomized sequences
(either via scrambling or $(u+x) \bmod 1$ tricks, where $u \sim U(0,1)$). I do not
need the functionality to be comprehensive; I want to provide a set of robust
tools that SAS IML customers can use to estimate high-dimensional integrals.
The integrals would have to be transformable to the unit hypercube.

Create an outline of the functionality that I should definitely include for a
minimal viable implementation of QMC and randomized QMC integration. Think about
implementing the features in stages. For example, perhaps the first release uses
only a certain type of QMC sequences, with additional options being implemented
in later stages. Put the "must have" features in the early stages, and the "nice
to have" features in the later stages.

Lastly, I work best when there is existing source code that I can browse and
modify. Currently, I am thinking of two functions:

1. **`QRANDGEN`:** Return $n$ observations from a quasi-random method in $d$
	 dimensions. The syntax is:

	 ```sas
	 x = j(n, d);
	 CALL QRANDGEN(x, "Method", n, param1, param2, ...);
	 ```

	 This will fill `x` with $d$-dimensional quasi-random numbers from the `"Method"`
	 type (e.g., Sobol, Korobov lattice, etc.).

2. **`QMCQUAD`:** Performs quasi-Monte Carlo estimation of an integral. The syntax
	 is:

	 ```sas
	 call QMCQUAD(result, "Funcname", "Method", n||d, <L>, <U>, other_params);
	 ```

	 Here, `L` and `U` are optional limit vectors (defaults $\{0,0,\ldots,0\}$ and
	 $\{1,1,\ldots,1\}$) that specify the bounds for a $d$-dimensional hypercube.
	 The user-defined function `"Funcname"` is sent quasi-random vectors and the
	 function returns an estimate which is the average of the function values
	 evaluated at the quasi-random points. Optionally, return a standard error for
	 the integral.

Help me plan this, Gemini! You're my only hope!

## Gemini Response: Implementation Planning

This is a high-value addition to SAS IML. Providing a robust, native QMC/RQMC
subsystem will significantly upgrade customers' ability to handle
high-dimensional pricing models, Bayesian computations, and risk simulations
where standard Monte Carlo converges too slowly.

Here is a blueprint for your Minimum Viable Product (MVP) and subsequent stages,
along with refinements to your function architectures and the best open-source
repositories to study.

### Staged Implementation Plan

To avoid getting bogged down in the vast landscape of QMC sequences, focus on the
industry standards first.

#### Stage 1: The MVP (Must-Haves)

Your goal here is to deliver a functional, correct subsystem that covers 90% of
practical use cases.

**Sequences:**

- **Sobol:** The undisputed workhorse of QMC. You must include this. Use the Joe
	and Kuo direction numbers (up to 21,201 dimensions).
- **Halton:** Computationally simple to implement using prime bases. Excellent
	for lower dimensions ($d < 20$) and a great baseline for testing your
	architecture.

**Randomization (RQMC):**

- **Random Shift (Cranley-Patterson Rotation):** $y = (x+u) \bmod 1$, where
	$u \sim U(0,1)^d$. It is trivial to implement and works perfectly for lattice
	rules and Halton sequences (though less optimal for Sobol).

**Integration (`QMCQUAD`):**

- **Basic hypercube transformation:** Map $[0,1]^d$ to $[L,U]$.
- Provide the point estimate.

#### Stage 2: Robust RQMC & Error Estimation

Standard QMC cannot easily provide a standard error (SE). Randomized QMC solves
this by creating $m$ independent randomizations of the $n$-length sequence,
yielding $m$ independent integral estimates.

**Sequences:**

- **Rank-1 Lattices (Korobov):** Include a pre-computed table of good generator
	vectors (e.g., from L'Ecuyer or Nuyens). They are exceptionally good for
	periodic functions.

**Randomization:**

- **Digital Scrambling (Linear Matrix Scrambling):** Random shifts break the net
	properties of digital sequences like Sobol. You need a base-2 scrambling method
	(like Owen's or Matoušek's affine scrambling) to properly randomize Sobol nets
	while preserving their low-discrepancy properties.

**Integration (`QMCQUAD`):**

- **Standard Error:** Update the function to accept $m$ replications. Compute the
	integral for each of the $m$ randomized batches, return the mean of the $m$
	estimates as the final integral, and use the sample variance of the $m$
	estimates to compute the SE.

#### Stage 3: Advanced Options (Nice-to-Haves)

**Sequences:** Niederreiter, Faure (often requested by academics, but rarely
outperform optimally scrambled Sobol in practice).

**Integration:**

- **Adaptive Stopping:** Instead of fixed $n$, allow the user to specify an
	absolute or relative tolerance. The routine doubles the number of points (for
	Sobol, $n$ should always be a power of 2) until the SE falls below the tolerance.
- **Non-Rectangular Domains:** Built-in inverse CDF transformations to map
	$[0,1]^d$ to $\mathbb{R}^d$ (e.g., multivariate normal) for financial pricing.

### Refining Your Function Architectures

Your proposed syntax is very close to ideal, but requires a few tweaks to support
RQMC effectively.

#### 1. `QRANDGEN`

Since SAS IML supports dynamic matrix resizing, consider returning the matrix
rather than passing an uninitialized matrix as an argument, or pass $n$ and $d$
explicitly so the function can allocate it.

**Proposed syntax:**

```sas
call QRANDGEN(x, "Method", n, d, <"RandType">, <Seed>);
```

- **`"Method"`:** `"Sobol"`, `"Halton"`, `"Lattice"`.
- **`"RandType"`:** `"None"`, `"Shift"`, `"Scramble"`. (Defaults to `"None"`.)
- **`Seed`:** Optional integer for reproducibility of the randomization.

#### 2. `QMCQUAD`

To calculate a standard error, the function must handle $m$ independent
randomizations.

**Proposed syntax:**

```sas
call QMCQUAD(result, "Funcname", "Method", n, d, <m>, <L>, <U>, <other_params>);
```

- **`result`:** A $1 \times 2$ vector containing `[Estimate, StandardError]`.
- **`n`:** Base number of points per randomization (if using Sobol, enforce or
	warn if $n$ is not a power of 2, as Sobol properties degrade heavily otherwise).
- **`m`:** Number of independent randomizations (default to 10 or 30). If $m = 1$,
	the SE in `result` is set to missing (`.`).

**Implementation detail:** The routine generates $m$ batches of size $n$. It
passes each batch to `"Funcname"`, averages the results for that batch, and then
computes the grand mean and standard deviation across the $m$ batch averages.

### Source Code Repositories to Study

Rather than reading textbooks to decipher the bitwise operations required for
Sobol, rely on these established, battle-tested open-source implementations:

#### SciPy's `scipy.stats.qmc` (Python/Cython)

**Why it's perfect for you:** It is arguably the cleanest, most modern, and
best-documented implementation of QMC available today. It includes Halton, Sobol
(with Joe/Kuo direction numbers), and Owen scrambling.

**Where to look:** Go to the [SciPy GitHub repository](https://github.com/scipy/scipy)
and navigate to `scipy/stats/qmc`. The source code is highly readable and directly
translatable to matrix/array languages.

#### Frances Kuo and Stephen Joe's Sobol Direction Numbers (C++)

**Why it's perfect for you:** This is the academic gold standard for generating
Sobol sequences up to 21,201 dimensions.

**Where to look:** [Sobol direction numbers](https://web.maths.unsw.edu.au/~fkuo/sobol/).
They provide raw C++ code that demonstrates exactly how to apply the XOR bitwise
operations using their pre-calculated direction numbers.

#### Dirk Nuyens' "Magic Point Shop" (MATLAB / Python / C++)

**Why it's perfect for you:** If you implement Lattice rules (Stage 2), you need
good generator vectors. Nuyens provides the algorithms and the pre-computed
vectors.

**Where to look:** Search for "Dirk Nuyens Magic Point Shop". The MATLAB code
provided is incredibly short and will map almost 1:1 with SAS IML matrix
operations.

#### QMCpy (Python)

**Why it's perfect for you:** Developed by the IIT team (Fred Hickernell's group),
it focuses heavily on the adaptive integration side of things (`QMCQUAD`).

**Where to look:** [github.com/QMCSoftware/QMCpy](https://github.com/QMCSoftware/QMCpy).
Study their `stopping_criterion` modules to see how they implement adaptive error
tolerances.

## User Prompt: Testing Correctness

What are some ways that we might test the correctness of a quasi-random sequence?
For example, suppose I tell you that $\{a_1, a_2, \ldots, a_n\}$ is a quasi-random
sequence. What statistical test would you use to determine whether I am telling
the truth? Or, more specifically, suppose I tell you the sequence is a Sobol
sequence. Are there properties that you could check? (Let's not consider
reproducing the numbers from the definition of a Sobol sequence.)

## Gemini Response: Testing Strategies

Standard statistical tests for randomness (like the Kolmogorov-Smirnov, runs, or
Ljung-Box tests) are counterproductive for quasi-random sequences. Because QMC
sequences are designed to be hyper-uniform, they lack the natural clustering and
voids of true randomness. A standard goodness-of-fit test will often reject a QMC
sequence precisely because its empirical CDF matches the theoretical uniform CDF
too perfectly.

To verify a sequence, you must test for hyper-uniformity and structural
stratification rather than stochastic randomness.

### 1. The Structural $(t,m,s)$-Net Test (Sobol-Specific)

Sobol sequences are base-2 $(t,s)$-sequences, meaning their primary defining
characteristic is how perfectly they balance points across binary partitions of
the hypercube.

If you take a block of $N = 2^m$ Sobol points in $s$ dimensions, the sequence forms
a $(t,m,s)$-net. This means you can partition the unit hypercube into elementary
intervals (boxes) of volume $2^{t-m}$, and each box will contain exactly $2^t$
points.

You can programmatically verify this via "Property A" and its higher-dimensional
extensions:

- **1D Balance (Property A):** Take $N = 2^m$ points. Divide the interval $[0,1)$
	into $2^m$ equal sub-intervals. If it is a valid Sobol sequence, every single
	sub-interval will contain exactly one point. You can test this by computing
	$\lfloor N \cdot x_i \rfloor$ and checking that the result is a permutation of
	the integers $0,1,\ldots,N-1$.
- **2D Balance:** Take $N = 2^m$ points. Divide dimension 1 into $2^{k_1}$ segments
	and dimension 2 into $2^{k_2}$ segments, such that $k_1 + k_2 = m$. (For example,
	if $N = 16$, cut the space into $4 \times 4$ grids, or $8 \times 2$, or
	$2 \times 8$.) Every resulting rectangle of area $1/16$ must contain exactly one
	point (assuming the defect parameter $t = 0$ for those dimensions).

If a sequence fails these strict spatial counting tests for $N = 2^m$, it is not
a valid Sobol sequence.

### 2. $L_2$-Star Discrepancy Computation

While calculating the exact $L_\infty$ star discrepancy is an NP-hard geometric
problem, the $L_2$-star discrepancy can be computed analytically in $O(n^2)$ time
using Warnock's formula.

For a set of $n$ points $\{x_1,\ldots,x_n\}$ in $d$ dimensions, Warnock's formula
evaluates the $L_2$ distance between the empirical point distribution and the
uniform distribution. By evaluating this on your sequence and comparing it
against the known discrepancy of a standard pseudo-random generator, the QMC
sequence should yield a significantly lower value. A robust unit test suite
computes the $L_2$ discrepancy for $n \in \{100,500,1000\}$ and verifies it falls
below a threshold derived from the theoretical $O(n^{-1}(\log n)^d)$ bound.

### 3. Empirical Rate of Convergence

The most pragmatic test of a QMC sequence is its behavior in integration. The
error of Monte Carlo integration decays at $O(1/n)$, while QMC decays at roughly
$O(1/n)$ for sufficiently smooth functions.

1. Select a standard suite of integration test functions (the Genz test functions
	 are the industry standard for this, as they allow you to tune continuity,
	 oscillatory behavior, and corner peaks).
2. Compute the integral approximation $\hat{I}_n$ for $N = 2^k$ where
	 $k = 4,5,\ldots,16$.
3. Calculate the absolute error $\epsilon_n = |\hat{I}_n - I|$.
4. Perform a simple linear regression of $\log(\epsilon_n)$ on $\log(n)$.

If the sequence is genuinely quasi-random and the function is smooth, the slope
of that regression line will be distinctly steeper than $-0.5$ (typically between
$-0.8$ and $-1.0$). If the slope hovers around $-0.5$, the sequence has lost its
low-discrepancy properties.

