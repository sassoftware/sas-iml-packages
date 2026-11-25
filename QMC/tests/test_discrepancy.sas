/* L2StarDiscrepancy: Computes the L2-star discrepancy for a set of points.

   Purpose:
   Quantifies how well a sequence is uniformly distributed within the unit
   hypercube [0,1)^d. Because QMC sequences are designed to be hyper-uniform,
   standard goodness-of-fit tests are often inappropriate. This metric directly
   evaluates the structural stratification of the point set. Lower values indicate
   better uniformity.

   Formula (Warnock's Formula):
   The analytical L2-star discrepancy is evaluated in O(N^2) time using:
   (D_N^{*(2)})^2 = 3^{-d} - (2^{1-d}/N) * Sum_{i=1 to N} Prod_{j=1 to d} (1 - x_{ij}^2)
                    + (1/N^2) * Sum_{i=1 to N} Sum_{k=1 to N} Prod_{j=1 to d} [1 - max(x_{ij}, x_{kj})]

   References:
   Warnock, T. T. (1972). "Computational investigations of low-discrepancy 
   point sets." In Applications of Number Theory to Numerical Analysis (pp. 319-343).
   Lemieux, C. (2009). Monte Carlo and Quasi-Monte Carlo Sampling, Springer, p. 183.
   
   Inputs:
      x - N x d matrix of points in [0,1)^d.
   Output:
      The L2-star discrepancy (scalar).
*/
start L2StarDiscrepancy(x);
   N = nrow(x);
   d = ncol(x);

   /* Term 1 */
   term1 = 3##(-d);

   /* Term 2: Vectorized product over columns */
   Y = 1 - x##2;
   term2 = (2##(1-d) / N) * sum( Y[,#] );  /* sum of product of each row */

   /* Term 3: N^2 pairwise comparisons */
   term3_sum = 0;
   do i = 1 to N;
      max_x = x <> x[i, ];
      Y = 1 - max_x;
      t3_prod = sum( Y[,#] );  /* sum of product of each row */
      term3_sum = term3_sum + sum(t3_prod);
   end;
   term3 = term3_sum / (N##2);

   L2_sq = term1 - term2 + term3;
   if L2_sq < 0 then return( 0 );
   return( sqrt(L2_sq) );
finish;


/* TEST 1: Empirical Convergence Rate (The "Slope")
   Consider plotting $\log(\text{Discrepancy})$ on the y-axis against $\log(N)$ on the x-axis.
   Because the discrepancy is roughly proportional to $N^k$, taking the log of both sides gives
   $\log(\text{Discrepancy}) \approx k \cdot \log(N) + C$.
   The "slope" we are calculating is that exponent $k$, which represents the rate of convergence.

   What this test does:
   It estimates the rate at which the discrepancy decreases as the number
   of points (N) increases. It calculates the local slope on a log-log
   scale between the current N and the previous N.

   What the output should show:
   Standard pseudo-random numbers converge at O(N^{-0.5}), so their
   slope will stabilize near -0.5. A correctly implemented QMC sequence
   converges faster, so its empirical slope should be noticeably steeper
   (more negative) than -0.5, typically moving toward -0.8 to -1.0.
*/

bases = {2, 3, 5, 7, 11}; /* Test in d=5 dimensions */
N_vals = 100 // T(do(1000, 5000, 1000));
results = j(nrow(N_vals), 3, .);
do i = 1 to nrow(N_vals);
   N = N_vals[i];
   pts = halton(bases, N);        /* generate Halton sequence */
   L2 = L2StarDiscrepancy(pts);   /* compute discrepancy */
   results[i, 1:2] = N || L2;
   
   /* Calculate empirical convergence rate from previous N if i > 1 */
   if i > 1 then do;
      log_N_diff = log(results[i,1]) - log(results[i-1,1]);
      log_L2_diff = log(results[i,2]) - log(results[i-1,2]);
      results[i, 3] = log_L2_diff / log_N_diff;
   end;
end;

cnames = {"N" "L2_Star_Discrepancy" "Empirical_Slope"};
print results[colname=cnames format=BESTD5.];
/*
call series(log10(N_vals), log10(results[,2]))
     label={'log10(N)' 'log10(Discrepancy)'}
     other = "lineparm slope=-0.5 x=2.5 y=-2 /curvelabel='MC Convergence';" grid={x y};
*/


/********************************/
/* compare the L_2 discrepancy for Halton to the expected L_2 discrepancy
   for a random point set  X ~ U(0,1)^d
*/

/*
   TEST 2: Discrepancy Ratio (QMC vs. Expected Monte Carlo)

   What this test does:
   It computes the theoretical expected L2-star discrepancy for a perfectly
   random U(0,1)^d point set of size N. It then divides a QMC sequence's
   discrepancy by this expected baseline.

   What the output should show:
   The ratio should be strictly less than 1.0. Furthermore, as N increases,
   this ratio should steadily decrease, proving that the QMC sequence is
   becoming progressively more uniform than random chance can achieve.
*/
d = nrow(bases) * ncol(bases);
results = j(nrow(N_vals), 5, .);

do i = 1 to nrow(N_vals);
   N = N_vals[i];
   pts = halton(bases, N);        /* generate Halton sequence */
   L2_Halton = L2StarDiscrepancy(pts);   /* compute discrepancy */

   /* Expected L2 Discrepancy for U(0,1)^d */
   E_L2_sq_MC = (2**(-d) - 3**(-d)) / N;
   E_L2_MC = sqrt(E_L2_sq_MC);   
   results[i, 1:4] = N || L2_Halton || E_L2_MC || (L2_Halton / E_L2_MC);
end;
/* The Koksma-Hlawka inequality states that the absolute error of
   numerical integration is bounded by the product of the function's variation
   and the sequence's star discrepancy.
   At $N=1000$, the Halton sequence's discrepancy is 44.4% of standard Monte Carlo.
   Thus, the standard MC method would require N*(1/0.444) ~ 5000 points to match the precision of Halton.
*/
N_MC_needs = N_vals # (1/results[,4]##2);
results[ , 5] = N_MC_Needs;

cnames = {"N" "Halton_L2" "Expected_MC_L2" "Ratio_Halton_to_MC" "Equiv N_MC"};
print results[colname=cnames format=BESTD5.];


/* Test discrepancy of Sobol sequence */
dim = 5;
numPts = 1000;
pts = sobol(dim, numPts-1);
L2_Sobol = L2StarDiscrepancy(pts);
E_L2_MC = sqrt((2##(-dim) - 3##(-dim)) / numPts);
isBetterThanMC = (L2_Sobol < 0.5 * E_L2_MC);
print L2_Sobol  E_L2_MC;
run check_test("sobol: 5D L2-star discrepancy well below MC", isBetterThanMC, 1, 0);