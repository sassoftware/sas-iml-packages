/* BEFORE running this example, store the modules in the package, as shown in Install_Pkg.sas */

/* Analytical tests for special cases of bivariate normal CDF.
   Specify b = {b1 b2} and Sigma, which is a 2-D covariance matrix.
   Then the call 
   prob = CDFMVN(b, Sigma, mu) returns the probability that a 
   bivariate normal vector with mean mu and covariance Sigma is less than or equal to b:
   prob = P(X1 < b1 & X2 < b2 | X ~ BVN(mu, Sigma)), where BVN stands for bivariate normal.
   If not specified, the default value of mu is {0 0}.
*/
proc iml;
load module=_all_;

tol = 1E-14; /* Tolerance for checking results */
/* Analytical Test 1: Orthant probabilities */
pi = constant('pi');
R = {1 0.5, 0.5 1};
b = {0 0};
prob = cdfmvn(b, R);
correct = 1/4 + arsin(R[1,2]) / (2*pi);
if abs(prob - correct) > tol then 
   print "ERROR: Orthant probability", (R[1,2])[L='rho'] prob correct;

/* Analytical Test 2: Diagonal matrices ==> Identity correlation */
b = {-1.2 0.7};
R = I(2);
prob = cdfmvn(b, R);
correct = cdf('normal',b[1]) # cdf('normal',b[2]);
if abs(prob - correct) > tol then 
   print "ERROR: Diagonal matrix", prob correct;

/* Test 3: Comparison Test with PROBBNRM */
b = { 1.2 -0.7};
R = {1 -0.34, -0.34 1};
prob = cdfmvn(b, R);
correct = probbnrm(b[1], b[2], R[1,2]);
if abs(prob - correct) > tol then 
   print "ERROR: Compare with PROBBNRM", prob correct;
print "--- DONE ---";

/* 2-D ANALYTICAL TESTS */
/* Analytical Test 1: Orthant probabilities */ 
rho_vec = do(-0.9, 0.9, 0.1);
R = I(2);
b = {0 0};
do i = 1 to ncol(rho_vec);
   rho = rho_vec[i];
   R[1,2] = rho;
   R[2,1] = R[1,2];
   prob = cdfmvn(b, R);
   correct = 1/4 + arsin(rho) / (2*pi);
   if abs(prob - correct) > tol then 
      print "ERROR: Orthant probability", rho prob correct;
end;

/* Analytical Test 2: Diagonal matrices ==> Identity correlation */
xy = expandgrid( -2:2, -2:2 );
R = I(2);
do j = 1 to nrow(xy);
   b = xy[j,];
   prob = cdfmvn(b, R);
   correct = cdf('normal',b[1]) # cdf('normal',b[2]);
   if abs(prob - correct) > tol then 
      print "ERROR: Diagonal matrix", prob correct;
end;

/* Comparison Test with PROBBNRM: Correlation {1 rho, rho 1} */
rho_vec = do(-0.9, 0.9, 0.1);
xy = expandgrid( -2:2, -2:2 );
R = I(2);
do j = 1 to nrow(xy);
   b = xy[j,];
   do i = 1 to ncol(rho_vec);
      rho = rho_vec[i];
      R[1,2] = rho;
      R[2,1] = R[1,2];
      prob = cdfmvn(b, R);
      correct = probbnrm(b[1],b[2], R[1,2]);
      if abs(prob - correct) > tol then 
         print "ERROR: Compare with PROBBNRM", prob correct; 
   end;
end;

print "--- DONE TESTS FOR BIVARIATE CDF ---";
QUIT; 
