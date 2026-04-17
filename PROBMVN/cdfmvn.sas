/* ------------------------------------------------------------------
   Main Function: cdfmvn_mod
   ------------------------------------------------------------------ */
/* 
   Return the MVN CDF for MVN(Sigma, mu) at each row of b by using quasi-Monte Carlo integration.
   This code implements the lattice rule for numerical integration as described in
   Genz, A. and Bretz, F. (2009), "Computation of Multivariate Normal and t Probabilities", 
   p. 45-48, especially Eqn. (4.16) on p. 48.

   The code uses the Korobov lattice rules with different prime numbers.
   The h value for the Korobov lattice rule is read from a hard-coded matrix, mat.
   I have not located the original reference for this matrix. I assume that Bretz ran experiments to get these h values. 
   The matrix is from Frank Bretz's code.

   Genz, A., & Bretz, F. (2009). "Computation of multivariate normal and t probabilities." 
        Springer Science & Business Media.
*/
   
proc iml;
/* If Sigma is a covariance matrix, then the MVN probability can be computed
   by transforming to the corresponding correlation matrix R.

   If D denotes the diagonal matrix which has the square roots of
   the diagonal entries for Sigma on its diagonal, the correlation matrix R is defined
   by Sigma = D*R*D. Then the transformation x = Dy reduces the general MVN probability to
   Phi_k(b; Sigma, mu) = Phi_k(D^{−1}(b-mu); R)
*/
start cdfmvn_mod(b, Sigma, mu=repeat(0,1,ncol(Sigma)));
   IsValid = mvn_IsValidParmsMVN(b, Sigma, mu);
   if ^IsValid then
      return(.);
   run mvn_StdizeCovToCorr(U, R, b, Sigma, mu);
   if ncol(Sigma) < 10 then
      tol = 1e-4; 
   else
      tol = 1e-3;
   /* use an outer loop over the rows of b matrix, which allows us to compute 
      probabilities for different limits of integration. This is not optimially 
      efficient, since does not re-use the lattice points across different rows of b. */
  prob = j(nrow(b),1,.);
  error = j(nrow(b),1,.);
  do row = 1 to nrow(b);
     U_row = U[row,];
     run cdfmvn_LR(p, e, R, U_row, tol);  
     prob[row] = p;
     error[row] = e;
  end;
  return prob;    /* currently, we don't return the error estimates */
finish;

start cdfmvn_LR(prob, error,       /* output values: probability and error estimate */
             R, _b, eps=0.0001);   /* input values:  correlation matrix R, upper limits b, absolute error eps */
  q=ncol(_b);
  ridge_eps = 1E-12;
  /* apply a ridge factor by adding a small value to the diagonal of R to ensure that it is positive definite. */
  C = t(root(R+ridge_eps*I(q)));        /* C is the lower-triangular Cholesky factor of R */

  /* scale the b vector and the rows of C to prevent needing to divide by C[i,i] inside loops. */
  v = vecdiag(C);
  b = _b / rowvec(v);
  C = C / v;

  /* The p_vector contains a sequence of prime numbers to be used as number 
     of points in the Korobov lattice rule. 
     The code is adaptive. It starts with 157 points. If the error is larger
     than eps, it increases the number of points to the next prime in the list, and so on.
  */
  p_vector={157 313 619 1249 2503 5003 10007 20011};

  /* The mat matrix contains the h values for the Korobov lattice rule
     for dimensions 2 to 32 (rows) and for the different prime numbers (columns).
     That is, mat[,1] contains the multipliers for 157 points, mat[,2] for 313 points, and so on.
     These values for h are precomputed constants. They are the integers that minimize the integration error
     for a given dimension (q) and number of points (p).
  */
  mat={ 1   1   1   1   1    1    1     1,  /* q = 2 (bivariate) */
       46 119 239 512 672 1850 3822  6103,  /* q = 3 (trivariate) */
       46  93 178 136 652 1476 2325  2894,
       17  51  73 197 792  792 1206  8455,
       18  51 104 165 792  380 1927  3629,
       18  80 102 175 253  162 2286  1752,
       11  70 161 303 306  363  343  1920,
       11  70 161 155 153  137  378   652,
       11  93 106  18 288  186   81   146,
       36  62  57  27 288  186  182   156,
       36  15  57  27  29   33   76   136,
       36  15  36  24 128   36   21    44,
        4  19  22  24  64   38   21    31,
       30  15  22  24  16   36   21   161,
       31   9  22  24  16   48   20   161,
       31   9  22  14  16   48   21    11,
        6  20   6  14  16   12   21    11,
        6   9   6  14  64   12   11    13,
        3   9   6  14  16   12   11    13,
        3   9   6   8  16    6   11    13,
        3  16   6   8  16    6    7    13,
        3  16   6   8  16    6    7    22,
        3  16   6   8  16    6    7    13,
        3  16   6   8   8    6    7    13,
        3  16   6   8   8    6    4    13,
        3  16   6   8   8    5    4    16,
        3  16   4   3   8    5    4    16,
        3   4   4   3   8    4    4    13,
        3   4   4   3   8    4    4    13,
        3   4   4   3   8    4    4     8,
        3   4   4   3   8    4    4     8};   /* q = 32 */

  /* The outer DO UNTIL loops control n and p, which are the number of random shifts and the number of points in the lattice rule, respectively. 
     The inner DO (k,j,i) loops compute the lattice rule estimate for a given n and p. */
  error = eps*10;  /* initialize error to a large value to ensure the loop starts */
  n=10;
  do until(n>50 | error<eps);
    index=1;
    do until(index=9 | error<eps);
      p = p_vector[index];
      h = mat[q-1,index];
      /* evaluate the integral by using quasi-Monte Carlo (QMC) method */
      call mvn_qmc_eval(integ, varsum, n, p, h, q, b, C);
      error = 3*sqrt(varsum/(n*(n-1)));  /* 3*SE ==> 99.7% confidence */
      index = index+1;
    end;
    /* Every time we increment n, we perform another p evaluations where p=p_vector[index] (such as 157).
       If we didn't hit error tolerance, increment n by 2, which ensures that we don't do a little more work but 
       we do substantially more work. */
    n=n+2;
  end;
  prob=integ;
  return;
finish;

start mvn_qmc_eval(integ, varsum,   /* output values: integ and varsum for the QMC evaluation of the lattice rule for fixed values of n and p */
               n, p, h, q, b, C);
  /* This module handles only the k, j, and i loops for fixed values of n and p. */
  integ=0;
  varsum=0;
  runif = j(1, q-1, .);
  y = j(1, q, 0);
  e = j(1, q, 0);
  e[1] = cdf("Normal", b[1]);

  /* generator vector z = {1, h, h**2, ..., h**(q-1)}  (MOD p) */
  vec=0:q-2;
  z = mod(j(1,q-1,h)##vec,p);

  do k=1 to n;
    latsum=0;
    call randgen(runif, "Uniform");
    /* (j*z/p) is uniform grid of points in the q-dimensional hypercube.
       Add rr to randomly shift the grid.
       The computation is repeated n times (see the k=1 to n loop) so we can calculate std dev and estimate error.
       The Baker's transformation w(y) = |2*y - 1| is used because lattice rules are most powerful when the 
       function being integrated is periodic and smooth. However, the MVN integrand is not periodic. 
       If you just used the lattice points directly, you would get a convergence rate of roughly O(1/n).
       By applying the Baker’s transformation, we "fold" the integration space. 
       This transformation ensures periodicity b/c the value at 0 is the same as the value at 1.
    */
    do j=1 to p;
      t = mod(runif+j#z/p, 1);
      w=abs( 2*t-1 ); 
      
      /* Robust domain handling and optimized inner product go here */
      do i=2 to q;
        p_val = w[i-1] * e[i-1];
        if p_val <= 0.5 then y[i-1] = quantile("Normal", max(p_val, 1E-15));
        else                 y[i-1] = squantile("Normal", max(1 - p_val, 1E-15));
        e[i] = cdf("Normal", b[i] - C[i, 1:i-1] * y[1:i-1]);
      end;
      f=e[#];
      latsum=latsum+(f-latsum)/j;
    end;
    /*  monitor convergence by using Welford's "online" updates.
        See https://blogs.sas.com/content/iml/2026/04/06/welford-mean-var.html */
    varsum = varsum + (k-1)*(latsum-integ)**2/k;
    integ = integ + (latsum-integ)/k;
  end;
finish;

store module=(cdfmvn_LR mvn_qmc_eval cdfmvn_mod);
QUIT;
