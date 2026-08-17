# NMF: Nonnegative Matrix Factorization in SAS/IML

## Description
This project provides a SAS/IML implementation of Nonnegative Matrix Factorization (NMF), including the main `nmf` subroutine and supporting helper functions. The package enables users to factorize nonnegative matrices into low-rank components using several algorithms (ALS, MUP, APG). Example scripts are provided to demonstrate usage and verify performance.

## Installation
To install the NMF package in your SAS/IML session, run the `Install_Pkg.sas` script:

```sas
%include 'Install_Pkg.sas';
```

This will load all necessary modules and definitions.

## Usage
After installation, you can call the main NMF subroutine as follows:

```sas
proc iml;
load module=_all_;

/* Example: Factorize a nonnegative matrix A */
A = {1 2 3, 4 5 6, 7 8 9};
k = 2;
call nmf(w, h, A, k);
print w h;
quit;
```

See the `examples/` directory for more advanced usage and test cases.

## Main Functions
- `nmf` — Main subroutine for NMF
- `proc_nmf` — Wrapper for PROC NMF
- `stdizeWH` — Standardizes factor matrices

## Documentation
Detailed documentation is provided in the source code and example scripts. (Optional: Add `nmf.pdf` for full documentation.)

## License
Copyright (c) 2026 SAS Institute Inc. All rights reserved.
