# MATLAB-to-IML Conversion Progress

## Scope and Conventions

- Translate the ten routines in `HandbookMC_Matlab/quasimc` into
  `HandbookMC_IML/quasimc`, preserving the source algorithms and point ordering.
  Initial translations follow the source loops; subsequent vectorization may
  replace those loops with equivalent matrix operations.
- Each source file defines its IML module and stores it with `store module=(...)`.
- Standalone tests live in `HandbookMC_IML/quasimc/tests` and share `check_test`.
  Successful checks are silent; failures print the test name and both results.
- Run tests from the repository root so their relative `%include` paths resolve.
- PROBMVN provides the eventual package model: separate implementation files,
  a definition/include file, an installation script, and a `tests` directory.
  This first pass adds draft modules and tests, not a released package installer.

## Running the Tests

In SAS, set the working directory to the repository root, then submit:

```sas
%include "QMC/HandbookMC_IML/quasimc/tests/test_all.sas";
```

Each individual `test_FuncName.sas` is also standalone: it opens `PROC IML`,
includes the required module definitions and test helper, performs its checks,
and quits. No package installation is required. Each source file stores its
module in the current IML storage library. Successful assertions are silent;
failures print `--- test name --- FAILS ---` and the actual and expected results.

The current `test_halton.sas` also uses `load module=_all_;`. Store the
modules in `convertbase.sas`, including `FirstPrimes`, in the same IML storage
library before running that test.

## Progress

| Order | MATLAB Source | IML Module | Status |
| --- | --- | --- | --- |
| 1 | `nbe.m` | `nbe.sas` | Translated; tests written, SAS execution complete | NOT USED
| 2 | `vdc.m` | `vdc.sas` | Vectorized; `vdc_orig` retained for comparison, SAS execution complete | DONE
| 3 | `halton.m` | `halton.sas` | Translated; tests written, SAS execution complete | DONE
| 4 | `hammersley.m` | `hammersley.sas` | Translated; tests written, SAS execution complete | DONE
| 5 | `faure.m` | `faure.sas` | Translated; tests written, SAS execution completed | DONE
| 6 | `sobol.m` | `sobol.sas` | Translated; tests written, SAS execution pending |
| 7 | `sobmat.m` | `sobmat.sas` | Translated; tests written, SAS execution pending |
| 8 | `cbe.m` | `cbe.sas` | Translated; tests written, SAS execution pending |
| 9 | `korobov.m` | `korobov.sas` | Translated; tests written, SAS execution pending |
| 10 | `extkorobov.m` | `extkorobov.sas` | Translated; tests written, SAS execution pending |

## Recent Updates (2026-10-30)

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
  point-by-point implementation remains in the same file as `vdc_orig`.
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

## Vectorization Opportunities

The revised `vdc` is the model: operate on all point indices together, reuse a
digit matrix, align positional weights instead of reversing digit columns, and
reduce weighted rows. The following opportunities were identified by scanning
the remaining IML implementations. They are proposals, not implemented or
benchmarked changes.

| Functions | Candidate improvement and constraints |
| --- | --- |
| `convertFromBase10`, `cbe` | Form a matrix or row of quotients using positional powers, then apply `floor` and `mod` to extract all digits together. This can replace the remaining per-digit loops. `convertFromBase10` already processes all input values together; its remaining loop is over digit positions, not points. Preserve leading-zero padding, the zero case, and the numeric-range policy when reusing the newer helper. |
| `faure` | Build digits for `T(0:numPoints)` with the bulk conversion helper instead of repeated division of trailing submatrices. Align the helper's most-significant-first columns with the transform's least-significant-first convention by changing the transform indexing. Replace repeated weight matrices and elementwise weighted sums with a single weight vector and matrix-vector products. The Pascal matrix can also be filled in batches over its valid triangular entries; retain the genuine recurrence across coordinates. |
| `sobol` | Reuse one bulk binary digit matrix and align its ordering with the corresponding columns of each direction block before multiplication. Replace replicated weights and weighted temporary matrices with matrix-vector products. Coordinate blocks are independent and could be combined into a larger blocked calculation, but retaining the short coordinate loop may use less memory. Preserve ordinary binary-index order and the `numPoints+1` rows. |
| `sobmat` | Batch conversion of the completed direction numbers to binary, replacing the final loop that calls `cbe` once per column. Align bits to each column's own digit position before packing the square blocks. The direction-number XOR recurrence depends on earlier values and is not removed by the `vdc` technique; verify intermediate direction numbers fit the helper's cutoff before reusing it. |
| `korobov`, `extkorobov` | Replace the two `repeat` matrices with outer products: `T(1:numPoints) * generator` for Korobov and `sequence * generator` for extended Korobov, retaining the existing normalization and `mod`. This primarily removes temporary arrays because point evaluation is already vectorized. Keep the modular generator recurrence; replacing it with unrestricted multiplier powers risks overflow or precision loss. `extkorobov` also inherits the faster `vdc`. |
| `nbe` | Locate the rightmost digit below `base-1`, increment it, and zero the trailing slice instead of propagating carry one position at a time. Handle all-maximal-digit overflow and empty trailing slices explicitly. Keep the operation on digit vectors rather than converting them to a potentially oversized numeric scalar. |
| `halton` | Immediately benefits from the vectorized `vdc`. Its remaining loop is over bases, not points. A shared digit-position calculation across bases is possible, but different bases need different digit widths; retain the small coordinate loop unless profiling justifies the larger batched representation. |
| `hammersley` | Already uses vector assignments and inherits improvements through `halton`. There is no local point loop to remove. Preserve the shifted Halton rows and the unused final base. |
| `logbase`, `numDigBase`, `convertToBase10` | Already use array operations, indexed assignments, and row reductions. No comparable loop-removal opportunity was found. `convertToBase10` could use a matrix-vector product instead of forming the weighted digit matrix before summing rows. |
| `FirstPrimes` | Already uses a lookup table and vector subscripting to select the first `n` primes. No prime-generation loop or comparable vectorization work is needed. |

Before implementing these candidates:

- Load `convertbase.sas` for the new `vdc` dependencies and `nbe.sas` when
  comparing against `vdc_orig`. The source defines both versions, but its
  current `store` statement stores only `vdc`.
- Compare results at a single point and at point counts just below, equal to,
  and just above base powers. Preserve digit order, output shape, zero rows,
  and the source-specific ordering documented below.
- Check missing returns from bounded helpers before subsequent matrix
  operations. The numeric cutoff is not a memory budget: full digit matrices
  and combined coordinate blocks can be large. Benchmark time and peak memory;
  use blocks of points where needed rather than requiring full materialization.

## Validation and Attention Items

- SAS execution is not configured in VS Code. Checks here are static checks and
  independent reference arithmetic. User-confirmed SAS results are recorded in
  the progress table (`nbe` is confirmed); the remaining SAS runs are pending.
- `vdc_orig` seeds its digit vector with zero, avoiding the earlier empty-seed
  dependency on `nbe`. IML null matrices still need special handling before
  copying them.
- `vdc` uses `convertFromBase10` for all indices at once and weights
  `base##(-(numDigits:1))` for a row-wise weighted sum. This removes its
  point loop and avoids reversing digit columns. `vdc_orig` retains the
  previous point loop and column reversal for comparison. Both return zero
  immediately for a single point.
- The MATLAB `hammersley` source shifts the Halton coordinates down one row and
  does not use the last supplied base. Preserve this convention in the draft;
  decide later whether to use the usual unshifted Hammersley construction.
- The MATLAB `faure` source returns `N+1` points, indexed from zero through `N`.
  Preserve that count rather than silently changing the argument's meaning.
- `faure` uses `repeat` for `repmat`, `comb` for `nchoosek`, and row-sum
  subscripting for `sum(...,2)`. The original matrix and loop operations remain.
- `sobol` also returns `N+1` points in ordinary binary-index order, not Gray-code
  order. Request `N=2^m-1` to obtain a complete block of `2^m` points.
- IML has no three-dimensional matrices. `sobmat(d,r)` therefore returns an
  `r` by `r*d` matrix of adjacent blocks. Dimension `k` occupies columns
  `(k-1)*r + (1:r)`; `sobol` reads that layout. This is the draft's one output
  representation change from MATLAB.
- The source Sobol table has 53 polynomials but only 20 rows of initial values;
  use dimensions 1 through 20. The tables are retained without replacement.
- `sobmat` uses SAS `bxor` for MATLAB `bitxor` and explicitly allocates space for
  direction-number recurrence beyond the eight supplied initial values. Confirm
  SAS bit-width limits before extending the draft to very large digit counts.
- `cbe` returns a zero digit explicitly for input zero, avoiding `log(0)`, and
  uses a local remainder so IML's reference arguments do not modify the caller.
  Tests include exact base-power boundaries and integer reconstruction.
- `korobov` retains the source's row order `1:N`, so its origin is last.
  `extkorobov` uses base-2 van der Corput order, so its origin is first. Both
  retain generator construction modulo `N`; dyadic extension is tested.
- The original ten draft translations still assume valid scalar sizes and
  integer bases as in the MATLAB sources. The support utilities now enforce
  the numeric cutoff above; broader input validation, generator selection,
  and extreme-size robustness remain later package work.
- Use a prime base at least as large as the dimension for Faure and pairwise
  coprime bases for Halton. Point counts are positive integers; `cbe` accepts
  nonnegative integers. These mathematical prerequisites are not yet validated
  by the draft functions.
- All MATLAB function calls used here have IML translations. No unresolved
  calls require `/* RICK NEEDS TO CONVERT THIS FUNCTION */` markers in this pass.

## Checks Completed Here

- Module definitions, `store` statements, and `DO`/`END` structure have been
  checked statically for all ten translations.
- Independent exact arithmetic verifies the expected sequence points, digit
  expansions, packed Sobol matrix examples and recurrence values, and lattice
  ordering and dyadic extension fixtures. This does not execute the IML code.
- The bounded logarithmic digit-count reference check covered 42,080 small,
  cutoff, and base-power cases across bases 2 through 541. It modeled SAS
  `CEIL`'s `1E-12` near-integer tolerance; this was not SAS execution or a proof
  for every supported input. Cutoff digit fixtures and `EXACTINT` rejection
  expectations were also checked independently.
- All ten sources have corresponding standalone tests. The suite preserves the
  requested order, with the support-utility test appended, and all suite and
  dependency include paths resolve.
- The copied Sobol tables match all 53 polynomial values and all 160 initial
  values in the MATLAB source.
- No untranslated MATLAB calls remain. Static structure and editor diagnostic
  checks pass; whitespace is checked separately from runtime correctness.
- The remaining SAS tests, including the new `vdc` comparison, still need to
  run in the user's configured SAS/IML environment. Apart from the user-confirmed
  `nbe` result in the progress table, no successful IML execution is claimed.

### Resolved Reference-Check Error

The earlier statement "an extra hand-written reference matrix had an incorrect
row" referred to a temporary expected matrix in a Python reference check. For
the second Sobol dimension with nine binary digits, the fourth row of that
hand-entered 9-by-9 direction matrix was `000111100`; the correct row is
`000100010`.

The auxiliary expectation was replaced by binomial-coefficient parity:
$V_{ij} = \binom{j-1}{i-1} \bmod 2$ for $i \le j$, and zero otherwise, using
one-based indices. The revised reference check passed. No MATLAB source, IML
module, or SAS test fixture needed changing because of this error. This was not
a SAS test run; execution of the SAS test suite remains pending.
