/*
IML tests that call the PROBMVN for the following situations:
1. 2-D problems. Include finite and semi-infinite domains of integration.
2. Degenerate problem. Create 3-D and 4-D examples that have a coordinate where the domain of integration is (-Infty, Infty), and thus it reduces to a 0-D, 1-D, or 2-D problems. Again, include both finite and semi-infinite limits for the variables that remain.
3. Degenerate problems that reduce to a diagonal matrix when the (-Infty, Infity) coordinates are excluded.
*/

/* Bivariate normal probabilities on rectangular domains. */
options ps=32000 nodate nonumber;

/* Tabulate how the built-in functions CDFMVN and PROBMVN work:
   Invalid Opt   CDFMVN/PROBMVN return a missing value.

   TESTNAME      Function  dim=2        dim=3         R=diag(5)    General dim>3
   opt={.,1}     CDFMVN    NO err est;  NO err est;   err est=0;   YES err est
                 PROBMVN   NO err est;  YES err est;  err est=0;   YES err est
*/

proc iml;
load module = _all_;

/* define the problems for various sizes */
* 2-D;
R2 = {1 0.5,
     0.5 1};
L2 = {-1.0 -0.5};
U2 = { 0.8  1.2};
correctCDF2 =  0.7327305;
correctPROB2 = 0.3849111;

* 3-D;
R3 = {1 0.5 0.3,
      0.5 1 0.4,
      0.3 0.4 1};
L3 = {-1.0 -0.5 -0.2};
U3 = { 0.8  1.2  0.5};
correctCDF3 =  0.5511043;
correctPROB3 = 0.1130838;

* 5-D with diagonal Sigma;
D5 = diag({5,4,3,2,1});
L5 = {-1.0 -0.5 -0.2 -0.1 -0.3};
U5 = { 0.8  1.2  0.5  0.7  0.9};
correctCDFD5 = 0.1603164;
correctPROBD5 = 0.0015286;

* General 5-D;
R5 = {1 0.8 0.6 0.4 0.3,
      0.8 1 0.7 0.5 0.4,
      0.6 0.7 1 0.6 0.5,
      0.4 0.5 0.6 1 0.7,
      0.3 0.4 0.5 0.7 1};
correctCDF5 =  0.5036553;
correctPROB5 = 0.0334044;

/*==============================*/
/* 1. 2-D test cases            */
/*==============================*/
print "--- Starting Tests for CDFMVN_MOD and PROBNRM_MOD option ---";

/* INVALID OPT */
opt = {1,.};
testName = "Test 1: Invalid Option: opt={1,.}";
cdf = cdfmvn_mod (U2, R2, , opt);
if cdf ^= . then
   print "CDFMVN_MOD did not catch an invalid option (dim=2).";
prob = probmvn_mod(L2, U2, R2, , opt);
if prob ^= . then
   print "PROBMVN_MOD did not catch an invalid option (dim=2).";
cdf = cdfmvn_mod (U5, R5, , opt);
if cdf ^= . then
   print "CDFMVN_MOD did not catch an invalid option (dim=5).";
prob = probmvn_mod(L5, U5, R5, , opt);
if prob ^= . then
   print "PROBMVN_MOD did not catch an invalid option (dim=5).";

/* RETURN ERROR ESTIMATE FOR LARGE-ENOUGH PROBLEM SIZES */
opt = {.,1};

testName = "Test 2: 2-D CDFMVN/PROBMVN does not return an error estimate: opt={.,1}";
cdf = cdfmvn_mod(U2, R2, , opt);
if ncol(cdf) = 2 then
   print "CDFMVN_MOD 2-D did not return one column.";
run check_test(testName, cdf[,1], correctCDF2);
prob = probmvn_mod(L2, U2, R2, , opt);
if ncol(prob) = 2 then
   print "PROBMVN_MOD 2-D did not return one column.";
run check_test(testName, prob[,1], correctPROB2);

testName = "Test 3: 3-D CDFMVN does not return an error estimate: opt={.,1}";
cdf = cdfmvn_mod(U3, R3, , opt);
if ncol(cdf) ^= 1 then
   print "CDFMVN_MOD 2-D did not return one column.";
run check_test(testName, cdf[,1], correctCDF3);
prob = probmvn_mod(L3, U3, R3, , opt);
if ncol(prob) ^= 2 then
   print "PROBMVN_MOD 3-D did not return two columns.";
run check_test(testName, prob[,1], correctPROB3);

testName = "Test 4: Diagonal Covariance returns exact zero for error estimate";
cdf = cdfmvn_mod(U5, D5, , opt);
run check_test(testName, cdf[,1], correctCDFD5, 1E-6);
if ncol(cdf) ^= 2 then
   print "CDFMVN_MOD Diagonal 5-D did not return two columns.";
if cdf[,2] ^= 0 then
   print "CDFMVN_MOD Diagonal 5-D did not return an exact zero for the error estimate.";
prob = probmvn_mod(L5, U5, D5, , opt);
run check_test(testName, prob[,1], correctPROBD5, 1E-6);
if ncol(prob) ^= 2 then
   print "PROBMVN_MOD Diagonal 5-D did not return two columns.";
if prob[,2] ^= 0 then
   print "PROBMVN_MOD Diagonal 5-D did not return an exact zero for the error estimate.";

testName = "Test 5: General Covariance returns a QMC error estimate";
cdf = cdfmvn_mod(U5, R5, , opt);
run check_test(testName, cdf[,1], correctCDF5);
if ncol(cdf) ^= 2 then
   print "CDFMVN_MOD 5-D did not return two columns.";
if cdf[2] <= 0 then
   print "CDFMVN_MOD 5-D did not return a positive error estimate.", cdf;
prob = probmvn_mod(L5, U5, R5, , opt);
run check_test(testName, prob[,1], correctPROB5);
if ncol(prob) ^= 2 then
   print "PROBMVN_MOD 5-D did not return two columns.";
if prob[2] <= 0 then
   print "PROBMVN_MOD 5-D did not return a positive error estimate.", prob;

print "--- Completed Tests for CDFMVN_MOD and PROBNRM_MOD option ---";
quit;

