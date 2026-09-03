proc iml;

    
/* Helper module copied from probmvn_tests.sas */
start check_test(test_name, prob, correct, tol=1E-3);
   maxDiff = max(abs(prob-correct));
   if maxDiff = 0.0 then do;
      msg = cat("--- ",test_name, " passes exactly ---");
      print msg[L=""];
   end;
   else if maxDiff > tol then do;
      msg = cat("--- ",test_name, " FAILS ---");
      print msg[L=""], maxDiff prob correct;
   end;
   else do;
      msg = cat("--- ",test_name, " passes ---");
      print msg[L=""];
   end;
finish;

/*==============================================================*/
/*    Clip limit into [-delta, delta], where delta ~ 8.125 is   */
/*    chosen so that                                            */
/*    SDF("Normal", delta) ~ constant("maceps")                 */
/*==============================================================*/
start ClipInterval(x, delta=8.125);
   return( -delta <> (x >< delta) );
finish;

/*==============================================================*/
/* IML Modules from blog post */
/* SECOND ARTICLE: COMPUTE BVN PROBABILITY OVER 
   ANY REGION (a,b) x (c,d), where we allow 
   a = -Infinity, b = +Infinity, c = -Infinity, 
   and d = +Infinity 
   If we use the missing values .M and .I to represent 
   -Infinity and +Infinity, respectively, then there 
   are 16 possible regions to consider.
*/
/***********************************************/

/* DEFINE AND STORE THE CDFBN and PROBBVN FUNCTIONS */
/* Extend the bivariate CDF, which is PROBBNRM(a,b, rho) to support an 
   upper limit of infinity (.I) for either argument:
   Pr(u,.I; rho) = Phi(u) = cdf("Normal", u) is probability over left half plane
   Pr(.I,v; rho) = Phi(v) = cdf("Normal", v) is probability over lower half plane
*/
start CDFBN(u,v,rho);
   if missing(u) & missing(v) then 
      return 1;
   if missing(u) then 
      return cdf("Normal", v);
   if missing(v) then 
      return cdf("Normal", u);
   return probbnrm(u, v, rho);
finish;

start ProbBVN(a,b,c,d,rho);
   ma = missing(a);   mb = missing(b);
   mc = missing(c);   md = missing(d);
   /* 1. complete plane */
   if ma & mb & mc & md then
      return 1;
   /* 2. lower half plane */
   if ma & mb & mc & ^md then
      return CDFBN(.I,d,rho);
   /* 3. upper half plane */
   if ma & mb & ^mc & md then
      return 1 - CDFBN(.I,c,rho);
   /* 4. horiz strip */
   if ma & mb & ^mc & ^md then
      return CDFBN(.I,d,rho) - CDFBN(.I,c,rho);
   /* 5. left half plane */
   if ma & ^mb & mc & md then 
      return CDFBN(b,.I,rho); 
   /* 6. SW quadrant */
   if ma & ^mb & mc & ^md then
      return CDFBN(b,d,rho);
   /* 7. NW quadrant */
   if ma & ^mb & ^mc & md then 
      return CDFBN(b,.I,rho) - CDFBN(b,c,rho); 
   /* 8. left strip (W) */
   if ma & ^mb & ^mc & ^md then
      return CDFBN(b,d,rho) - CDFBN(b,c,rho);
   /* 9. right half plane */
   if ^ma & mb & mc & md then
      return 1 - CDFBN(.I,a,rho);
   /* 10. SE quadrant = lower half - (SW quad)*/
   if ^ma & mb & mc & ^md then
      return CDFBN(.I,d,rho) - CDFBN(a,d,rho);
   /* 11. NE quadrant = right half - (SE quad) */
   if ^ma & mb & ^mc & md then 
      return 1 - CDFBN(.I,a,rho)
               - (CDFBN(.I,c,rho) - CDFBN(a,c,rho));
   /* 12. right strip (E) */
   if ^ma & mb & ^mc & ^md then
      return CDFBN(.I,d,rho) - CDFBN(.I,c,rho) 
           - CDFBN(a,d,rho) + CDFBN(a,c,rho);
   /* 13. vert strip */
   if ^ma & ^mb & mc & md then
      return CDFBN(.I,b,rho) - CDFBN(.I,a,rho);
   /* 14. lower strip (S) */
   if ^ma & ^mb & mc & ^md then
      return CDFBN(b,d,rho) - CDFBN(a,d,rho);
   /* 15. upper strip (N) */
   if ^ma & ^mb & ^mc & md then
      return CDFBN(b,.I,rho) - CDFBN(a,.I,rho) /*upper strip _|_  */
           - CDFBN(b,c,rho) + CDFBN(a,c,rho);
   /* 16. rectangle */
   if ^ma & ^mb & ^mc & ^md then
      return CDFBN(b,d,rho) - CDFBN(a,d,rho) 
           - CDFBN(b,c,rho) + CDFBN(a,c,rho); 
   return( . );  /* should never execute this statement */
finish;

/* Use Monte Carlo simulation to estimate the probability that a 
   MVN random variable is less than a specified value in each coordinate.
   Sigma is a kxk covariance matrix; mu is an option row vector with k elements.
   The row vector b specifies the upper limits of integration. 
   The function returns an estimate of 
   P(X1<b[1] & X2<b[2] & ... & Xk<b[k] | X~MVN(mu, Sigma))
   by simulating N random variates from MVN(mu, Sigma) and returning the proportion
   that are in the specified region.
*/
start MC_CDFMVN(N, b, Sigma, mu=j(1,ncol(Sigma),0));
   X = randnormal(N, mu, Sigma);
   inRegion = j(N,1,1);
   do i = 1 to ncol(b);
      v = (X[,i] < b[i]);
      inRegion = inRegion & v;
   end;
   return mean(inRegion);
finish;

/* Monte Carlo simulation of P( L < X < U | X~MVN(mu,Sigma) )
   where L and U are row vectors and a missing value 
   represents -Infinity in L and represents +Infinity in U.
   NOTE: If N is very large, a straightforward implementation uses a lot of memory.
   Therefore, perform the computations in blocks of size MaxN to
   prevent out-of-memory errors. It actually speeds up the computation, too!
*/
start MC_PROBMVN(N, Lcov, Ucov, Sigma, mu=j(1,ncol(Sigma),0));
   MaxN = 2E5;
   /* standardize the parameters to the correlation scale */
   L = Xform_Limits_Cov2Corr(Lcov, Sigma, mu);
   U = Xform_Limits_Cov2Corr(Ucov, Sigma, mu);
   R = cov2corr(Sigma);
   N_remain = N;
   nIter = ceil(N / MaxN);
   sum = 0;
   do j = 1 to nIter;
      K = min(MaxN, N_remain);
      inRegion = j(K,1,1);
      X = randnormal(K, j(1,ncol(R),0), R);
      /* Note: The '<' operator works correctly with a missing value
         on the left. The expression (a<y & y<b) is correct if a=.;
         However, if b=., you need to use (a<y). */
      do i = 1 to ncol(L);
         if U[i]=. then 
            v = (L[i] < X[,i]);
         else
            v = ((L[i] < X[,i]) & (X[,i] < U[i]));
         inRegion = inRegion & v;
      end;
      N_remain = N_remain - K;
      sum = sum + sum(inRegion);
   end;
   return sum / N;
finish;

/* return a list with the MC est and a 95% CL. The list looks like
   [prob, lower95, upper95] */
start MC_PROBMVN_CL(N, L, U, Sigma, mu=j(1,ncol(Sigma),0));
   prob_MC = MC_PROBMVN(N, L, U, Sigma, mu);
   SE_MC = sqrt( prob_MC * (1-prob_MC)/N );
   Lower95 = prob_MC - 1.96*SE_MC;
   Upper95 = prob_MC + 1.96*SE_MC;
   return( [prob_MC, Lower95, Upper95] );
finish;

start CreateAR1(N, rho);
   Sigma = rho##distance(T(1:N), T(1:N), "L1"); /* AR(1) */
   return( Sigma );
finish;

store module=(
check_test
ClipInterval
CDFBN
ProbBVN
MC_CDFMVN
MC_PROBMVN
MC_PROBMVN_CL
CreateAR1
);
QUIT;