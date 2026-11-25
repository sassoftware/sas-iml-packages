# MATLAB-to-IML Conversion Progress

## Scope and Conventions

- Translate the ten routines in `HandbookMC_Matlab/quasimc` into
  `HandbookMC_IML/quasimc`, preserving the source algorithms and point ordering.
  Initial translations follow the source loops; subsequent vectorization may
  replace those loops with equivalent matrix operations.
- Each source file defines its IML module and stores it with `store module=(...)`.
- Individual tests live in `HandbookMC_IML/quasimc/tests` and share `check_test`.
  Successful checks are silent; failures print the test name and both results.
- Run tests from `HandbookMC_IML/quasimc/tests` so their relative `%include`
  paths resolve. Include local test files directly and implementation files
  through `../`.
- PROBMVN provides the eventual package model: separate implementation files,
  a definition/include file, an installation script, and a `tests` directory.
  This first pass adds draft modules and tests, not a released package installer.

### IML Indexing Rule

In IML, a linear subscript such as `x[Idx]` **always returns a column vector**,
even when `x` is a row vector. To select columns of a row vector while keeping
the result a row vector, use `x[ , Idx]`.

This resolved the difficult coefficient-shape bug in `sobmat`:

```sas
coeffs = coeffs[, 2:ncol(coeffs)];
```

The former linear selection changed the orientation, making the following
`degree = ncol(coeffs)` incorrect. Preserve explicit two-dimensional subscripts
where later calculations depend on row or column shape.

## Running the Tests

Set the SAS working directory to `QMC/HandbookMC_IML/quasimc/tests`. The tests
include `check_test.sas` locally and implementation files from the parent
directory, for example `../convertbase.sas`. Run the aggregate test file with:

```sas
%include "test_all.sas";
```

No package installation is required. Tests that depend on digit expansion now
include `../convertbase.sas` and then `../cbe.sas` before their dependent module
definitions. Halton and Hammersley still load stored modules; the broader
clean-session setup remains under REVIEW below. The source files and
`check_test.sas` store their modules in the current IML storage library.

Successful `check_test` assertions are silent; failures print
`--- test name --- FAILS ---`, the maximum difference, and the computed and
expected results. Expected missing returns are tested through missing-value
masks because the helper treats a missing computed value as failure.

**REVIEW:** Confirm a reproducible clean-session setup and dependency order for
`test_all.sas` without relying on previously stored modules. The user reports
that the suite passes against the refactored files; whether that run began in
a fresh session remains unconfirmed.

## Progress

Implementation status was reconciled on **2026-11-24** from all eleven SAS
files and the commit history. Validation was updated on **2026-11-25** after
the user ran the test suite against the refactored files and reported that
all tests pass. These are user-confirmed SAS runs, not executions performed
by the assistant in VS Code. The table now includes post-refactor validation.

| Order | MATLAB Source | IML Module | Current State | SAS Validation |
| --- | --- | --- | --- | --- |
| 0 | utilities | `convertbase.sas` | helper functions; bounded integer policy retained | complete |
| 1 | `nbe.m` | `nbe.sas` | NOT USED by the active sequence implementations; retained utility | complete |
| 2 | `vdc.m` | `vdc.sas` | DONE; vectorized over all point indices | complete |
| 3 | `halton.m` | `halton.sas` | DONE; vectorized VDC per base | complete |
| 4 | `hammersley.m` | `hammersley.sas` | DONE; vector assignments and Halton coordinates | complete |
| 5 | `faure.m` | `faure.sas` | DONE; vectorized Pascal construction, digits, and weighted evaluation | complete |
| 6 | `cbe.m` | `cbe.sas` | DONE; vectorized scalar/vector digit conversion | complete |
| 7 | `sobmat.m` | `sobmat.sas` | DONE; packed generating matrices and coefficient-shape fix | complete |
| 8 | `sobol.m` | `sobol.sas` | DONE; vectorized binary digits and weighted evaluation | complete |
| 9 | `korobov.m` | `korobov.sas` | Draft lattice algorithm unchanged; program REVIEW remains | complete; suite passes |
| 10 | `extkorobov.m` | `extkorobov.sas` | Draft lattice algorithm unchanged; program REVIEW remains | complete; suite passes |

## Refactoring Update (2026-11-25)

- On branch `refactor`, `cbe` is now the single active implementation of
  nonnegative-integer base-b expansion. Its most-significant-first output and
  leading-zero padding are unchanged; commented `_orig` functions are untouched.
- `convertFromBase10` delegates to `cbe` after its existing raw-input cutoff
  check and rounding. VDC calls `cbe` directly, checks the scalar missing return,
  and obtains the width from `ncol(digits)` instead of recalculating it.
- Faure and Sobol also call `cbe` and derive their widths from the digit matrix.
  Faure reverses the rows and columns of its small Pascal transform and reverses
  its weight order. Sobol reverses the input columns of each direction block,
  leaving the output weights unchanged. Neither reverses the large digit matrix;
  point ordering and coordinate recurrences are preserved.
- Added `../cbe.sas` to the six test entry points that newly require it.
  Existing test cases and the tests-directory working convention are unchanged.
- Independent exact-arithmetic checks matched VDC fixtures and converter
  rounding cases, Faure old/new outputs and ternary fixtures, and Sobol old/new
  outputs through 20 dimensions and nine digits, including the recurrence
  fixture. Static checks confirmed one active expansion implementation and
  resolved all 58 test include paths. These reference/static checks did not
  execute SAS.
- The user subsequently ran the test suite against the refactored SAS files
  on 2026-11-25 and reported that all tests pass. Post-refactor SAS validation
  is complete for the suite; the remaining REVIEW items are separate from
  that successful run.

## Recent Updates (2026-11-24)

- Reviewed the current sources against the previous documentation commit
  `53ad794` (2026-11-06) and implementation commit `b0f880d` (2026-11-24).
  The earlier VDC, Halton, and Hammersley work is recorded in commits
  `02c4da7` and `4610bbe`. The lattice algorithms remain draft work; argument
  renaming and formatting changes are not treated as algorithm revisions.
- `cbe` now accepts a scalar or vector, normalizes it with `colvec`, and uses
  positional powers, broadcast division, `floor`, and `mod` to extract all
  digits together. Results have one row per value, padded on the left with
  zeros. It obtains digit counts from `numDigBase` and returns missing if that
  result is missing. Tests now include binary and ternary batches.
- `faure` constructs the lower-triangular Pascal transform with row/column
  index matrices and a mask. Base-b digits are extracted in bulk, and a column
  weight vector replaces repeated weights and weighted temporary matrices.
  The genuine recurrence across coordinates remains.
- `sobol` similarly constructs its binary digit matrix with broadcast division
  and evaluates each direction block using a weight column vector. The short
  coordinate loop remains; the per-digit construction loop is gone.
- `sobmat` retains the XOR recurrence and flattened block layout. The crucial
  row-preserving coefficient selection is documented under the IML indexing
  rule above. Tests cover identity and single-digit blocks, recurrences beyond
  eight digits, and the twentieth coordinate.
- Sobol tests cover source-ordered points, prefixes across digit lengths,
  one-dimensional balance, a recurrence point beyond eight digits, and a
  16-point dyadic grid balance check.
- The old `vdc_orig`, `faure_orig`, `sobol_orig`, and `cbe_orig` versions are
  commented reference code, not active definitions or stored comparison
  modules. `nbe` is no longer a dependency of the active sequence pipeline.
- `tests/test_discrepancy.sas` provides an L2-star discrepancy diagnostic for
  Halton and Sobol, including comparisons with a Monte Carlo baseline. It is
  not included by `test_all.sas`; its integration and execution status are
  marked REVIEW below.

## Earlier Updates (2026-10-30)

- Added leading comments to all ten translated SAS files describing each
  function's purpose, input requirements, and output dimensions and ordering.
- Added `HandbookMC_IML/quasimc/convertbase.sas` with five support modules:
  `logbase`, `numDigBase`, `convertFromBase10`, `convertToBase10`, and
  `FirstPrimes`.
- Corrected `logbase`'s positive-input test to `all(x>0)`, made zero require
  one digit in `numDigBase`, and corrected the conversion example's function
  name to `convertFromBase10`.
- Added `HandbookMC_IML/quasimc/tests/test_convertbase.sas`, reusing
  `check_test`, and registered it as the final test in `test_all.sas`.
  Coverage includes the documented examples, default bases, vector and matrix
  shapes, rounding, input preservation, digit validity, and round trips.
- Vectorized `vdc` by converting the full index column to a digit matrix,
  applying reversed-power weights, and summing weighted rows. The former
  point-by-point implementation remains as commented `vdc_orig` reference code.
- Added `FirstPrimes(n)`, which returns an `n` by 1 column vector containing
  the first `n` primes from a fixed lookup table, for integer `1 <= n <= 1000`.
  Inputs outside that range return scalar missing (`.`). The helper is included
  in the utility file's `store module=(...)` statement.
- Extended `test_halton.sas` to use `bases = FirstPrimes(10)` and generate
  1,000 points in ten dimensions. Checks compare the first nonzero point with
  the reciprocals of the bases and the 100th point with explicit expected
  coordinates. This supplies distinct prime bases without hand-written lists;
  availability of 1,000 bases does not set a package-wide dimensionality limit.

### Bounded Integer Policy

The goal is fast support for quasi-Monte Carlo integration, not a
general-purpose number-theory library or arbitrary-precision arithmetic.
Dimensionality limits are a separate policy from the numeric range below;
the existing Sobol tables still support only 20 dimensions.

- `numDigBase`, `convertFromBase10`, and `convertToBase10` each declare
  `safe_cutoff = 1E9;`. Keep these three local settings synchronized if the
  limit changes. The cutoff itself is accepted.
- If any input exceeds `safe_cutoff`, `numDigBase` and `convertFromBase10`
  reject the whole call by returning scalar missing (`.`).
  `convertToBase10` does the same if any reconstructed result exceeds the
  cutoff. Callers must check for missing before continuing.
- `numDigBase` retains the fast, vectorized logarithmic formula for positive
  rounded inputs and returns one digit for zero. `logbase` remains unrestricted
  by the cutoff because digit counting evaluates it at `n+1`.
- SAS/IML uses double precision. `constant('exactint')` is the consecutive
  exact-integer limit, not the largest representable integer; it is `2##53`
  on this platform. Near that limit, logarithmic digit counts can round
  incorrectly even when the input integer is exactly representable.
  The `1E9` cutoff is a conservative practical choice for the logarithmic
  approach, not a claim about the largest theoretically safe cutoff.
- Tests now cover the inclusive cutoff and its predecessor, rejection just
  above the cutoff and in mixed valid/over-limit inputs, and powers of each
  base from 2 through 16 with their neighboring integers. The `EXACTINT` and
  `EXACTINT-1` examples now assert rejection rather than successful conversion.

## Vectorization Status

The shared approach is now implemented in several routines: build digits for
all points together, align positional powers or transform indexing, and use
matrix operations rather than per-point or per-digit loops. Completion here
describes the code, not a measured performance guarantee.

| Functions | Implemented | Remaining Opportunity |
| --- | --- | --- |
| `cbe` | Single source for batch digit extraction with broadcast positional powers | No duplicate active expansion remains |
| `vdc` | CBE digit matrix and reversed-power weights; no point loop or digit-column reversal | Optional matrix-vector weighted sum instead of the weighted temporary matrix |
| `faure`, `sobol` | CBE digit matrices and aligned transforms with weight-column matrix products | Coordinate loops remain intentionally; consider blocks only if profiling justifies them |
| `sobmat` | Uses vectorized `cbe`; coefficient orientation corrected | The final packing loop could batch conversions and align bits per column; keep the dependent XOR recurrence |
| `halton`, `hammersley` | Inherit vectorized VDC; Hammersley already uses vector assignments | Halton retains a small loop over bases; no local Hammersley point loop remains |
| `convertFromBase10` | Rounding and cutoff compatibility wrapper around `cbe` | No separate digit-expansion loop remains |
| `logbase`, `numDigBase`, `convertToBase10` | Array operations, indexed assignments, and row reductions | A matrix-vector product could avoid the weighted temporary in `convertToBase10` |
| `FirstPrimes`, `nbe` | Prime lookup is already vectorized; `nbe` is unused by active sequences | No current optimization priority for the legacy carry loop |
| `korobov`, `extkorobov` | Existing point evaluation is vectorized | Replace replicated matrices with outer products after reviewing the lattice routines; retain modular generator recurrence |

Future changes should preserve single-point behavior, base-power boundaries,
digit orientation, output shape, and source-specific ordering. Benchmark time
and memory: the numeric cutoff is not a memory budget, and complete digit or
multi-coordinate temporary matrices can be large. Block points where needed.

## Validation and Attention Items

- SAS execution is not configured in VS Code. The progress table records the user's
  successful SAS checks of the modified functions as of 2026-11-24.
- `vdc` uses `cbe` for all indices at once and weights
  `base##(-(numDigits:1))` for a row-wise weighted sum. This removes its
  point loop and avoids reversing digit columns. The active implementation
  returns zero immediately for a single point.
- The MATLAB `hammersley` source shifts the Halton coordinates down one row and
  does not use the last supplied base. Preserve this convention in the draft;
  decide later whether to use the usual unshifted Hammersley construction.
- `faure` returns `numPts+1` points, indexed from zero through `numPts`.
  Its digit matrix is most-significant-first from `cbe`, matched to a Pascal
  transform with reversed row/column order and reversed positional weights.
- `sobol` also returns `N+1` points in ordinary binary-index order, not Gray-code
  order. Its CBE digit columns are aligned with reversed input columns of each
  direction block. Request `N=2^m-1` to obtain a complete block of `2^m` points.
- IML has no three-dimensional matrices. `sobmat(d,r)` therefore returns an
  `r` by `r*d` matrix of adjacent blocks. Dimension `k` occupies columns
  `(k-1)*r + (1:r)`; `sobol` reads that layout. This is the draft's one output
  representation change from MATLAB, also called a flattened matrix format.
- The source Sobol table has 53 polynomials but only 20 rows of initial values;
  use dimensions 1 through 20. The tables are retained without replacement.
- `sobmat` uses SAS `bxor` for MATLAB `bitxor` and allocates space for recurrence
  beyond the eight supplied initial values. Its header records unsigned
  32-bit operands (`0..2^32-1`); the supported digit-count limit still needs
  review in conjunction with `cbe`'s cutoff.
- `cbe` converts scalar or vector inputs to a most-significant-first matrix
  without modifying the caller. Zero requires one digit, and shorter values
  are padded on the left. It depends on `numDigBase`, unlike the older scalar
  remainder-based implementation.
- `korobov` retains the source's row order `1:N`, so its origin is last.
  `extkorobov` uses base-2 van der Corput order, so its origin is first. Both
  retain generator construction modulo `N`; dyadic extension is tested.
- The active dependency chains include `halton -> vdc -> cbe`,
  `hammersley -> halton`, `faure -> cbe`, and both `sobol -> cbe` and
  `sobol -> sobmat -> cbe`. The compatibility path is
  `convertFromBase10 -> cbe -> numDigBase -> logbase`.
  `extkorobov` uses `vdc`. Load the utility modules and `cbe` before executing
  dependent routines; `cbe` does not call the compatibility converter.
- Use a prime base at least as large as the dimension for Faure and pairwise
  coprime bases for Halton. Point counts are positive integers; `cbe` accepts
  nonnegative integers. These mathematical prerequisites are not yet validated
  by the draft functions.
- All MATLAB function calls used here have IML translations.

### REVIEW Items

- **REVIEW: Lattice algorithms.** `korobov` and `extkorobov` remain draft
  implementations outside the completed refactoring pass. Their suite tests
  now pass with the revised helper setup; substantive algorithm and generator
  selection review remains future work.
- **REVIEW: Clean-session tests and dependency comments.** Establish the
  compilation/load order for the aggregate and individual tests.
- **REVIEW: Supported numeric limits and missing propagation.** The helper
  cutoff does not by itself validate every generator argument or intermediate
  direction number. Confirm the maximum Sobol/Sobmat digit count with the
  32-bit `bxor` range and `cbe`'s `1E9` bound. Confirm that generators stop
  safely when a bounded helper returns missing, before allocation or subscripting.
- **REVIEW: Discrepancy diagnostic status.** Confirm which cases in
  `tests/test_discrepancy.sas` were run and whether it should remain a manual
  diagnostic or join the suite. It assumes an existing IML session and loaded
  modules and is not registered in `test_all.sas`. Empirical discrepancy
  comparisons are diagnostics, not general integration-error guarantees.

## Checks Completed Here

- The 2026-11-24 documentation review read all eleven implementation files,
  compared the relevant commits, and inspected the affected test setup and
  newer fixtures. It did not rerun SAS or reproduce performance measurements.
- Earlier static checks covered module definitions, `store` statements, and
  `DO`/`END` structure for the ten translations. Earlier independent arithmetic
  checked sequence points, digit expansions, packed Sobol matrix examples and
  recurrences, and lattice ordering and dyadic extension fixtures. Those checks
  do not replace the user's later SAS runs or validate every revised code path.
- The bounded logarithmic digit-count reference check covered 42,080 small,
  cutoff, and base-power cases across bases 2 through 541. It modeled SAS
  `CEIL`'s `1E-12` near-integer tolerance; this was not SAS execution or a proof
  for every supported input. Cutoff digit fixtures and `EXACTINT` rejection
  expectations were also checked independently.
- All ten translated sources have corresponding test files. The aggregate
  retains the original include order and appends the support-utility test;
  the user reports that the suite passes. Clean-session dependency loading
  remains REVIEW rather than assumed complete.
- The copied Sobol tables match all 53 polynomial values and all 160 initial
  values in the MATLAB source.
- The user reports that all suite tests pass against the 2026-11-25 refactored
  SAS files. Post-refactor manual SAS validation is complete; only the separate
  unresolved REVIEW items remain open.
