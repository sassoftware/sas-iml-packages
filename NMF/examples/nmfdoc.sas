/*---------------------------------------------------------------
 |  nmfdoc.sas: Nonnegative Matrix Factorization (NMF) example
 |
 |  Purpose:
 |     Demonstrate the CALL NMF subroutine (SAS/IML implementation)
 |     on real data.
 |
 |  Description:
 |     1. Read the Virginica species of Fisher's Iris data into a
 |        nonnegative data matrix X and print the first five
 |        observations.
 |     2. Approximate X as X ~ W*H of rank 2 using CALL NMF with its
 |        default method.
 |     3. Print the first five rows of the W factor and the H factor.
 |     4. Form the rank-2 approximation Y = W*H and print its first
 |        five rows.
 |     5. Compute and print the overall root-mean-square residual
 |        of the approximation.
 |
 |  Requires:
 |     The NMF modules stored by nmf.sas (run nmf.sas first so that
 |     LOAD MODULE=_ALL_ can find them). Available in SAS 9.4 and later.
 ---------------------------------------------------------------*/
proc iml;
load module=_all_;

/* Apply NMF to the Virginica species of Fisher's Iris data */
use sashelp.iris where(Species='Virginica');
read all var _NUM_ into X[c=varNames];
close;
print (X[1:5,])[c=varNames L="Virginica Species (first five observations)"];

call randseed(1234, 1);           /* Set random number seed for reproducibility */
names = 'NMF 1':'NMF 2';
run nmf(W, H, X, ncol(names));   /* rank-2 non-negative factorization */
print (W[1:5,])[c=names l='W (First five rows)'],
      H[c=varNames r=names l='H'];

Y = W*H;  /* rank-2 approximation of X */
print (Y[1:5,])[c=varNames F=4.1 l="NMF Rank-2 Approximation to First Five Observations"];

RMS = norm(X-W*H)/sqrt(ncol(X)*nrow(X));
print RMS[F=4.2 l="Root Mean Square Residual of the NMF Approximation"];
quit;