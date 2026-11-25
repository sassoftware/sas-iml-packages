proc iml;
%include "../vdc.sas";
%include "../halton.sas";
%include "../hammersley.sas";
%include "../faure.sas";
%include "../sobol.sas";
%include "../sobmat.sas";
%include "../cbe.sas";
%include "../korobov.sas";
%include "../extkorobov.sas";
%include "../convertbase.sas";
%include "check_test.sas";
quit;

proc iml;
load module=_all_;

/*----------------------------------------------------------------------
  1. Test Integrands: Each accepts a 1 x d row vector x and returns scalar
  ----------------------------------------------------------------------*/

/* #1: Sobol' g-function (Sobol, 1990)
   c holds weights a_j >= 0; exact integral = 1 */
start f_sobolG(x, c, u);
   return( prod( (abs(4*x - 2) + c) / (1 + c) ) );
finish;

/* #2: Genz Oscillatory (Genz, 1984) 
   Smooth, periodic sign changes */
start f_oscillatory(x, c, u);
   twopi = 2 * constant('pi');
   return( cos(twopi * u[1] + x * T(c)) );
finish;

/* #3: Genz Product Peak 
   Smooth, sharp interior spike at u */
start f_prodPeak(x, c, u);
   return( prod( 1 / (c##(-2) + (x - u)##2) ) );
finish;

/* #4: Genz Corner Peak 
   Smooth, steep boundary peak at 0 */
start f_cornerPeak(x, c, u);
   d = ncol(x);
   return( (1 + x * T(c)) ## (-(d + 1)) );
finish;

/* #5: Genz Gaussian 
  smooth separable bell peak at u */
start f_gaussian(x, c, u);
   return( exp( -sum( (c # (x - u))##2 ) ) );
finish;

/* #6: Genz Continuous (C^0) 
   Non-differentiable creases along x_j = u_j */
start f_continuous(x, c, u);
   return( exp( -sum( c # abs(x - u) ) ) );
finish;

/*----------------------------------------------------------------------
  2. Exact Analytical Integrals on [0,1]^d
     Depends: cbe() (for Corner Peak vertex enumeration)
  ----------------------------------------------------------------------*/
start exactIntegral(funcName, c, u);
   d = ncol(c);
   if funcName = "f_sobolG" then
      return(1);

   if funcName = "f_oscillatory" then do;
      twopi = 2 * constant('pi');
      halfC = c / 2;
      return( cos(twopi * u[1] + sum(halfC)) * prod(sin(halfC) / halfC) );
   end;

   if funcName = "f_prodPeak" then
      return( prod( c # (atan(c # (1 - u)) + atan(c # u)) ) );

   if funcName = "f_cornerPeak" then do;
      /* Enumerate all 2^d vertices of {0,1}^d via cbe() */
      V = cbe(0:(2##d - 1), 2);
      signs = (-1) ## V[, +];
      denom = 1 + V * T(c);
      return( sum(signs / denom) / (gamma(d + 1) * prod(c)) );
   end;

   if funcName = "f_gaussian" then do;
      rtpi = sqrt(constant('pi'));
      terms = (rtpi / (2 * c)) # (erf(c # (1 - u)) + erf(c # u));
      return( prod(terms) );
   end;

   if funcName = "f_continuous" then do;
      terms = (2 - exp(-c # u) - exp(-c # (1 - u))) / c;
      return( prod(terms) );
   end;

   return(.);
finish;

/*----------------------------------------------------------------------
  3. Quadrature Wrapper: Loops over rows of an M x d point matrix
  ----------------------------------------------------------------------*/
start evalIntegrand(funcName, x, c, u);
   if funcName = "f_sobolG"      then return( f_sobolG(x, c, u) );
   if funcName = "f_oscillatory" then return( f_oscillatory(x, c, u) );
   if funcName = "f_prodPeak"    then return( f_prodPeak(x, c, u) );
   if funcName = "f_cornerPeak"  then return( f_cornerPeak(x, c, u) );
   if funcName = "f_gaussian"    then return( f_gaussian(x, c, u) );
   if funcName = "f_continuous"  then return( f_continuous(x, c, u) );
   return(.);
finish;

/*
   qmcQuadRow: Estimate integral on [0,1]^d from an M x d QMC point set.
   Inputs:
      funcName - Name of integrand module.
      pts      - M x d matrix of quasi-random points in [0,1]^d.
      c, u     - 1 x d parameter row vectors.
   Output: Scalar QMC integral estimate (sample mean of f(pts[i, ])).
*/
start qmcQuadRow(funcName, pts, c, u);
   M = nrow(pts);
   fVals = j(M, 1, .);
   do i = 1 to M;
      fVals[i] = evalIntegrand(funcName, pts[i, ], c, u);
   end;
   return( mean(fVals) );
finish;

/*----------------------------------------------------------------------
  4. Benchmark Driver: Compare Halton vs. Sobol across d in {2,3,5,10,20}
  ----------------------------------------------------------------------*/
primes20 = FirstPrimes(20);
dims = {2, 3, 5, 10, 20};
M = 4096; /* Power of 2 so Sobol forms an exact (t, 12, d)-net */

funcs = {"f_sobolG", "f_oscillatory", "f_prodPeak",
         "f_cornerPeak", "f_gaussian", "f_continuous"};
/* Standard Genz difficulty normalization totals b_k */
difficulty = {1.0, 4.5, 3.625, 0.925, 3.515, 10.2};

call randseed(202611);
summary = j(nrow(dims) * nrow(funcs), 7, .);
rowIdx = 0;

do dIdx = 1 to nrow(dims);
   d = dims[dIdx];
   ptsHalton = halton(primes20[1:d], M);
   ptsSobol  = sobol(d, M - 1); /* M points: indices 0 .. M-1 */

   uBase = randfun(1 || d, "Uniform");
   cBase = randfun(1 || d, "Uniform");

   do fIdx = 1 to nrow(funcs);
      rowIdx = rowIdx + 1;
      fname = funcs[fIdx];

      if fname = "f_sobolG" then
         c = 0:(d - 1); /* Increasing weights a_j = j - 1 */
      else if fname = "f_gaussian" then
         c = cBase * (difficulty[fIdx] / sqrt(ssq(cBase)));
      else
         c = cBase * (difficulty[fIdx] / sum(cBase));

      exactVal  = exactIntegral(fname, c, uBase);
      estHalton = qmcQuadRow(fname, ptsHalton, c, uBase);
      estSobol  = qmcQuadRow(fname, ptsSobol,  c, uBase);

      relErrH = abs(estHalton - exactVal) / abs(exactVal);
      relErrS = abs(estSobol  - exactVal) / abs(exactVal);

      summary[rowIdx, ] = d || fIdx || exactVal || estHalton || estSobol || relErrH || relErrS;
   end;
end;

cnames = {"d" "FuncID" "Exact" "Halton_Est" "Sobol_Est" "Halton_RelErr" "Sobol_RelErr"};
print summary[colname=cnames format=BESTD10.];
