/* BEFORE running this example, store the modules in the package, as shown in Install_Pkg.sas */

/* Analytical tests for special cases of multivariate normal CDF, dim > 3.
   For dim > 3, the CDFMVN function uses a quasi-Monte Carlo method, so 
   the function returns approximate probabilities that might vary from run to run.
   In dimension d:
   Specify b = {b1 b2 ... bd} and Sigma, which is a d-D covariance matrix.
   Then the call 
   prob = CDFMVN(b, Sigma, mu) returns the probability that a 
   d-dimensional normal vector with mean mu and covariance Sigma is less than or equal to b:
   prob = P(X1 < b1 & X2 < b2 & ... & Xd < b_d | X ~ MVN(mu, Sigma)), where MVN stands for multivariate normal.
   If not specified, the default value of mu is {0 0 ... 0}.
*/
proc iml;
load module=_all_;

/* --- TEST SUITE --- */
/* Helper module to format test results */
start check_test(test_name, prob, correct, tol=0.001);
   maxDiff = max(abs(prob-correct));
   if maxDiff > tol then do;
      msg = cat("--- ",test_name, " FAILS ---");
      print msg[L=""], maxDiff prob correct;
   end;
   else do; 
      msg = cat("--- ",test_name, " passes ---");
      print msg[L=""];
   end;
finish;

call randseed(12345);

print "--- Starting Test Suite for CDFMVN ---";

/* 1. 4-D Identity Matrix */
test_name = "Test 1: 4-D Identity Matrix";
R = I(4);
b = {0 -1 -2 3};
correct = prod( cdf("Normal", b) );
prob = cdfmvn(b, R);
run check_test(test_name, prob, correct);

/* 2. 5-D Identity Matrix */
test_name = "Test 2: 5-D Identity Matrix";
R = I(5);
b = {1 1 1 1 1};
correct = prod( cdf("Normal", b) );
prob = cdfmvn(b, R);
run check_test(test_name, prob, correct);

/* 3. 5-D Rank-1 Update (Equicorrelated) 
   Let R = 0.5*I + 0.5*1*1'. This is equivalent to rho=0.5 
   For rho=0.5, b=0, the prob is 1/(dim+1) = 1/6 = 0.16666...
*/
test_name = "Test 3: 5-D Equicorrelated (rho=0.5)";
v = j(5,1, sqrt(0.5));
R = 0.5*I(5) + v*v`;
b = {0 0 0 0 0};
correct = 1/6;
prob = cdfmvn(b, R);
run check_test(test_name, prob, correct);

/* 4. 5-D Min Matrix */
test_name = "Test 4: 5-D Min Matrix";
Sigma = {1 1 1 1 1, 
         1 2 2 2 2, 
         1 2 3 3 3, 
         1 2 3 4 4, 
         1 2 3 4 5};
b = 0:4;
correct = 0.4597946;
prob = cdfmvn(b, Sigma);
run check_test(test_name, prob, correct);

/* 5. 8-D Min Matrix */
test_name = "Test 5: 8-D Min Matrix";
Sigma = {1 1 1 1 1 1 1 1, 
         1 2 2 2 2 2 2 2, 
         1 2 3 3 3 3 3 3, 
         1 2 3 4 4 4 4 4, 
         1 2 3 4 5 5 5 5, 
         1 2 3 4 5 6 6 6, 
         1 2 3 4 5 6 7 7, 
         1 2 3 4 5 6 7 8};
b = 0:7;
correct = 0.4590496;
prob = cdfmvn(b, Sigma);
run check_test(test_name, prob, correct);

/* 6. A kxk equicorrelated matrix with rho=0.5 and b=0.
      Theoretical prob at (0,0,...,0) is 1/(k+1) 
*/
kk = {9, 12, 15};
do i=1 to nrow(kk); 
   k = kk[i];
   v = j(k,1,sqrt(0.5));
   R = 0.5*I(k) + v*v`;
   b = j(1,k,0);
   correct = 1/(k+1);
   prob = cdfmvn(b, R);
   /* for large matrices, reduce the desired precision */
   test_name = cat("Test 6: Equicorrelated (rho=0.5), dim=",strip(char(k,2)));
   run check_test(test_name, prob, correct, 0.01);
end;


/* 7. Rank-1 singular correlation matrix.
   If R is a matrix of all 1s, then X1=X2=...=Xn.
   P(X1 < 1, X2 < 2, ..., X8 < 8) = P(X1 < min(b)) = P(X1 < 1).
*/
test_name = "Test 7: 8-D Singular (All 1s)";
k = 8;
R = j(k,k,1);
b = 1:k;
correct = cdf("Normal", min(b));
prob = cdfmvn(b, R);
run check_test(test_name, prob, correct);


/* 8. Block-diagonal correlation matrix. If the blocks are 1x1 and 2x2,
      then the MVN probability is the product of univariate and bivariate probs.
*/
test_name = "Test 8: Block-Diagonal Matrix";
/* Unit Test: 5-D Block Diagonal Correlation Matrix */
R =  { 1.0  0.5  0.0  0.0  0.0,
       0.5  1.0  0.0  0.0  0.0,
       0.0  0.0  1.0  0.0  0.0,
       0.0  0.0  0.0  1.0 -0.3,
       0.0  0.0  0.0 -0.3  1.0 };
b = {0.5  0.8  0.4  1.0 -0.2};
/* Calculate the correct value using the product of components */
p1 = probbnrm(b[,1], b[,2], 0.5);   /* Bivariate Block 1 */
p2 = cdf("Normal", b[,3]);          /* Univariate Block 2 */
p3 = probbnrm(b[,4], b[,5], -0.3); /* Bivariate Block 3 */
correct = p1 # p2 # p3;
prob = cdfmvn(b, R);
run check_test(test_name, prob, correct);

/* 9. 4-D example from Bretz */
test_name = "Test 9: 4-D Example from Bretz";
R={1          0.7071068    0          0,
   0.7071068  1            0.5        0,
   0          0.5          1          0.3333333,
   0          0            0.3333333  1};
b={1 1 1 1};
correct = 0.583122;
prob = cdfmvn(b, R);
run check_test(test_name, prob, correct);

/* 10. Same example, but with a covariance matrix instead of correlation. 
       Also, add a mean vector. The result should be the same as the previous test. 
*/
test_name = "Test 10: 4-D Example from Bretz on Covariance Scale with Mean";
v = {0.5, 1, 2, 4};
D = diag(v);
mu = {10 20 30 40};
Sigma = D*R*D;
b_new = mu + v`#b;  /* b_new are limits on the covariance scale that corresponds to b on the correlation scale */
prob = cdfmvn(b_new, Sigma, mu);
/* correct is the same as in the previous problem */
run check_test(test_name, prob, correct);

/* 11. Diagonal covariance matrix with non-zero mean. 
       The result should be the product of univariate probabilities. 
*/
test_name = "Test 11: 4-D Diagonal Covariance+Mean";
Sigma = diag({1, 4, 9, 16});
mu    = {10 20 30 40};
b     = mu + rowvec(sqrt(vecdiag(Sigma)));  /* Evaluates at z=1 for all */
correct = cdf("Normal", 1)**4;
prob    = cdfmvn(b, Sigma, mu);
run check_test(test_name, prob, correct);


test_name = "Test 12: 4-D Block-Diagonal Covariance+Mean";
Sigma = { 2  1   0   0,
          1  2   0   0,
          0  0   3 -1.5,
          0  0 -1.5  3 };
mu    = {-5  5  0 10};
b     = mu + rowvec(sqrt(vecdiag(Sigma)));
/* Block 1 rho: 1/2 = 0.5. Block 2 rho: -1.5/3 = -0.5 */
correct = probbnrm(1, 1, 0.5) * probbnrm(1, 1, -0.5);
prob    = cdfmvn(b, Sigma, mu);
run check_test(test_name, prob, correct);


test_name = "Test 13: 5-D Block-Diagonal Covariance+Mean";
Sigma = { 4  2  0  0  0,
          2  4  0  0  0,
          0  0  9  0  0,
          0  0  0  5  4,
          0  0  0  4  5 };
mu    = {1  2  3  4  5};
b     = mu + rowvec(sqrt(vecdiag(Sigma)));
/* Block 1 rho: 0.5. Block 2: Univariate. Block 3 rho: 4/5 = 0.8 */
correct = probbnrm(1, 1, 0.5) * cdf("Normal", 1) * probbnrm(1, 1, 0.8);
prob    = cdfmvn(b, Sigma, mu);
run check_test(test_name, prob, correct);


test_name = "Test 14: 6-D Block-Diagonal Covariance+Mean";
Sigma = { 2  1.6  0  0  0   0,
          1.6 2   0  0  0   0,
          0  0    3  0  0   0,
          0  0    0  4  0   0,
          0  0    0  0  5 -2.5,
          0  0    0  0 -2.5  5 };
mu    = {10  20  30  40  50  60};
b     = mu + rowvec(sqrt(vecdiag(Sigma)));
/* Block 1 rho: 0.8. Blocks 2/3: Univariate. Block 4 rho: -0.5 */
correct = probbnrm(1, 1, 0.8) * cdf("Normal", 1)**2 * probbnrm(1, 1, -0.5);
prob    = cdfmvn(b, Sigma, mu);
call check_test(test_name, prob, correct);
 

test_name = "Test 15: Multiple Rows of b (Limits of Integration)";
/* 5-D Block Diagonal Correlation Matrix */
R =  { 1.0  0.5  0.0  0.0  0.0,
       0.5  1.0  0.0  0.0  0.0,
       0.0  0.0  1.0  0.0  0.0,
       0.0  0.0  0.0  1.0 -0.3,
       0.0  0.0  0.0 -0.3  1.0 };

b = {0.5  0.8  0.4  1    -0.2 ,
     0.5  0.8  0.4  1    -0.2 ,
     0.6  0.9  0.5  1.1  -0.1 ,
     1    1.3  0.9  1.5   0.3 ,
     0.8  0.4  1   -0.2   0.5 ,
     0.4  1   -0.2  0.5   0.8 ,
     1   -0.2  0.5  0.8   0.4 ,
    -0.2  0.5  0.8  0.4   1    };

prob = cdfmvn(b, R);
/* compare with calling the rows one at a time */
correct = j(nrow(b),1,.);
do i = 1 to nrow(b);
   correct[i] = cdfmvn(b[i,], R);
end;
/* abs(prob - correct) should be small. They use different 
   quasi-random numbers, so they aren't exactly equal */
run check_test(test_name, prob, correct);


QUIT;
