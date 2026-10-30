# MATLAB-to-IML Conversion Progress

## Scope and Conventions

- Translate the ten routines in `HandbookMC_Matlab/quasimc` into
  `HandbookMC_IML/quasimc`, preserving the source algorithms and loop structure.
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

## Progress

| Order | MATLAB Source | IML Module | Status |
| --- | --- | --- | --- |
| 1 | `nbe.m` | `nbe.sas` | Translated; tests written, SAS execution pending |
| 2 | `vdc.m` | `vdc.sas` | Translated; tests written, SAS execution pending |
| 3 | `halton.m` | `halton.sas` | Translated; tests written, SAS execution pending |
| 4 | `hammersley.m` | `hammersley.sas` | Translated; tests written, SAS execution pending |
| 5 | `faure.m` | `faure.sas` | Translated; tests written, SAS execution pending |
| 6 | `sobol.m` | `sobol.sas` | Translated; tests written, SAS execution pending |
| 7 | `sobmat.m` | `sobmat.sas` | Translated; tests written, SAS execution pending |
| 8 | `cbe.m` | `cbe.sas` | Translated; tests written, SAS execution pending |
| 9 | `korobov.m` | `korobov.sas` | Translated; tests written, SAS execution pending |
| 10 | `extkorobov.m` | `extkorobov.sas` | Translated; tests written, SAS execution pending |

## Validation and Attention Items

- SAS execution is not configured in VS Code. No successful SAS test run has
  been completed; all tests require execution in the user's SAS/IML environment.
  Further checks here will be static checks and independent reference arithmetic.
- `nbe` explicitly handles the empty digit-vector seed before copying its input:
  IML cannot assign an unset/null matrix to another matrix.
- `vdc` reverses its digit vector by column subscripting instead of `fliplr`.
  The single-point case returns before forming a `1:0` subscript range.
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
- These draft translations assume valid scalar sizes and integer bases as in
  the MATLAB sources. Input validation, generator selection, and numerical
  robustness for extreme sizes are later package work, not part of this pass.
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
- All ten sources have corresponding standalone tests. The suite preserves the
  requested order, and all suite and dependency include paths resolve.
- The copied Sobol tables match all 53 polynomial values and all 160 initial
  values in the MATLAB source.
- No untranslated MATLAB calls remain. Static structure and editor diagnostic
  checks pass; whitespace is checked separately from runtime correctness.
- The SAS test suite still needs to be run in the user's configured SAS/IML
  environment. No claim of successful IML execution is made.

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
