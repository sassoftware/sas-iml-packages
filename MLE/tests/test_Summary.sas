proc iml;
%include "MLE_define.sas";
QUIT;

proc iml;
load module=_ALL_;
use sashelp.heart;
   read all var "Systolic";
close;

L_gamma = MLE_Fit("Gamma", Systolic);

run MLE_Summary(L_gamma);     /* print basic results */
run MLE_Summary(L_gamma, 2);  /* print all results */
QUIT;

/*******************************************/
/* Test Suite for MLE_Summary Function    */
/* Tests all printOpt levels (0-2)        */
/*******************************************/
proc iml;
load module=_ALL_;
call randseed(54321);

/* ========================================
   TEST 1: Normal Distribution - printOpt=0 (Minimal)
   Expected: Parameter estimates only
   ======================================== */
/* Generate clean normal data */
N1 = 500;
mu_true1 = 10;
sigma_true1 = 3;
Y1 = j(N1, 1, .);
call randgen(Y1, "Normal", mu_true1, sigma_true1);

L1 = MLE_Fit("Normal", Y1);
print "TEST 1: Normal Distribution, True parameters: mu=10, sigma=3";
run MLE_Summary(L1, 0);


/* ========================================
   TEST 2: Normal Distribution - printOpt=1 (Basic)
   Expected: Estimates + Fit Criteria
   ======================================== */
/* Use data with missing values to test CompleteCases */
Y2 = {3, 4, 8, 3, 7, ., 5, 3, ., 3, 6, 6, 1, ., 2, 4, 5};
print "TEST 2: Data with missing values (will be removed)";
print Y2[format=4.1];

L2 = MLE_Fit("Normal", Y2);
run MLE_Summary(L2, 1);


/* ========================================
   TEST 3: Normal Distribution - printOpt=2 (Extended)
   Expected: Basic + Optimization details
   ======================================== */
/* Large sample for good convergence */
N3 = 1000;
mu_true3 = 5;
sigma_true3 = 2;
Y3 = j(N3, 1, .);
call randgen(Y3, "Normal", mu_true3, sigma_true3);

L3 = MLE_Fit("Normal", Y3);
print "TEST 3: Normal Distribution, True parameters: mu=5, sigma=2";
run MLE_Summary(L3, 2);


/* ========================================
   TEST 4: Exponential Distribution - printOpt=1
   ======================================== */
/* Single parameter distribution */
N4 = 800;
scale_true = 4;
Y4 = j(N4, 1, .);
call randgen(Y4, "Expo", scale_true);

L4 = MLE_Fit("Expo", Y4);
print "TEST 4: Exponential Distribution, True parameter: scale=4";
run MLE_Summary(L4, 1);


/* ========================================
   TEST 5: Gamma Distribution - printOpt=2
   ======================================== */
/* Two-parameter distribution with different characteristics */
N5 = 600;
alpha_true = 3;
lambda_true = 2;
Y5 = j(N5, 1, .);
call randgen(Y5, "Gamma", alpha_true, lambda_true);

L5 = MLE_Fit("Gamma", Y5);
print "TEST 5: Gamma Distribution, True parameters: alpha=3, lambda=2";
run MLE_Summary(L5, 2);


/* ========================================
   TEST 6: Weibull Distribution - printOpt=1
   ======================================== */
/* Reliability/survival analysis distribution */
N6 = 500;
c_true = 2.5;
lambda_true_weib = 5;
Y6 = j(N6, 1, .);
call randgen(Y6, "Weibull", c_true, lambda_true_weib);

L6 = MLE_Fit("Weibull", Y6);
print "TEST 6: Weibull Distribution, True parameters: c=2.5, lambda=5";
run MLE_Summary(L6, 1);


/* ========================================
   TEST 7: Beta Distribution - printOpt=2
   ======================================== */
/* Bounded distribution (0,1) */
N7 = 400;
alpha_beta = 2;
beta_beta = 5;
Y7 = j(N7, 1, .);
call randgen(Y7, "Beta", alpha_beta, beta_beta);

L7 = MLE_Fit("Beta", Y7);
print "TEST 7: Beta Distribution, True parameters: alpha=2, beta=5";
run MLE_Summary(L7, 2);


/* ========================================
   TEST 8: Small Sample Normal - printOpt=1
   ======================================== */
/* Test with very small sample to see behavior */
Y8 = {12.5, 14.3, 11.8, 13.9, 12.1, 15.2, 13.5, 14.8};
print "TEST 8: Small sample (n=8)";
print Y8[format=5.1];

L8 = MLE_Fit("Normal", Y8);
run MLE_Summary(L8, 1);


/* ========================================
   TEST 9: Lognormal Distribution - printOpt=1
   ======================================== */
/* Right-skewed distribution */
N9 = 500;
mu_ln = 1;
sigma_ln = 0.5;
Y9 = j(N9, 1, .);
call randgen(Y9, "Lognormal", mu_ln, sigma_ln);

L9 = MLE_Fit("Lognormal", Y9);
print "TEST 9: Lognormal Distribution, True parameters: mu=1, sigma=0.5";
run MLE_Summary(L9, 1);


/* ========================================
   TEST 10: Normal Distribution - printOpt=2 (Extended)
   Expected: Estimates + fit criteria + optimization details
   ======================================== */
/* Normal data with good sample size */
N10 = 500;
Y10 = j(N10, 1, .);
call randgen(Y10, "Normal", 20, 5);

L10 = MLE_Fit("Normal", Y10);
print "TEST 10: Normal Distribution, True parameters: mu=20, sigma=5";
run MLE_Summary(L10, 2);


/* ========================================
   TEST 11: Skewed Data - printOpt=2
   Expected: Extended output for right-skewed distribution
   ======================================== */
/* Generate right-skewed data using Gamma */
N11 = 400;
Y11 = j(N11, 1, .);
call randgen(Y11, "Gamma", 0.5, 2);  /* Low alpha creates right skew */

L11 = MLE_Fit("Gamma", Y11);
print "TEST 11: Gamma Distribution, True parameters: alpha=0.5, lambda=2 (Low alpha creates right skew)";
run MLE_Summary(L11, 2);


/* ========================================
   TEST 12: Data with Outliers - printOpt=2
   Expected: Extended output with optimization details
   ======================================== */
/* Generate data with some outliers */
N12 = 100;
Y12 = j(N12, 1, .);
call randgen(Y12, "Normal", 10, 2);
/* Add some outliers */
Y12[95] = 25;  /* High outlier */
Y12[96] = 26;  /* High outlier */
Y12[97] = -5;  /* Low outlier */

L12 = MLE_Fit("Normal", Y12);
print "TEST 12: Normal Distribution with 3 artificial outliers";
run MLE_Summary(L12, 2);


/* ========================================
   TEST 13: Comparison Across All printOpt Levels
   Expected: Progressive detail from level 0 to 2
   ======================================== */
/* Same data, different detail levels */
N13 = 300;
Y13 = j(N13, 1, .);
call randgen(Y13, "Normal", 0, 1);

L13 = MLE_Fit("Normal", Y13);

print " ";
print "--- Level 0: Minimal ---";
run MLE_Summary(L13, 0);

print " ";
print "--- Level 1: Basic ---";
run MLE_Summary(L13, 1);

print " ";
print "--- Level 2: Extended ---";
run MLE_Summary(L13, 2);


/* ========================================
   TEST 14: Beta Distribution - printOpt=2
   Expected: Extended output for bounded distribution (0,1)
   ======================================== */
/* Beta distribution has bounded support */
N14 = 300;
Y14 = j(N14, 1, .);
call randgen(Y14, "Beta", 2, 5);

L14 = MLE_Fit("Beta", Y14);
print "TEST 14: Beta Distribution, True parameters: alpha=2, beta=5 (Support on (0,1))";
run MLE_Summary(L14, 2);


/* ========================================
   TEST 15: Lognormal - Heavy-tailed Data with printOpt=2
   Expected: Extended output with optimization details
   ======================================== */
/* Lognormal with large sigma creates heavy tails */
N15 = 400;
Y15 = j(N15, 1, .);
call randgen(Y15, "Lognormal", 1, 1.5);  /* Large sigma creates heavy tails */

L15 = MLE_Fit("Lognormal", Y15);
print "TEST 15: Lognormal Distribution, True parameters: mu=1, sigma=1.5 (Large sigma creates heavy tails)";
run MLE_Summary(L15, 2);


/* ========================================
   TEST 16: Small Sample with printOpt=2
   Expected: Extended output with small sample
   ======================================== */
Y16 = {12.5, 14.3, 11.8, 13.9, 12.1, 15.2, 13.5, 14.8, 12.9, 13.2};

L16 = MLE_Fit("Normal", Y16);
print "TEST 16: Small sample (n=10)";
run MLE_Summary(L16, 2);


/* ========================================
   TEST 17: Confidence Intervals - 95% (default)
   Expected: Parameter table with CI columns
   ======================================== */
Y17 = {9.8, 10.2, 10.1, 9.9, 10.0, 10.3};
L17 = MLE_Fit("Normal", Y17);
run MLE_Summary(L17, 1, 1);   /* printOpt=1, showCI=1 */


/* ========================================
   TEST 18: Confidence Intervals - 90% (alpha=0.10)
   Expected: Wider CI than 95%
   ======================================== */
Y18 = {4.8, 5.1, 4.9, 5.0, 5.2};
L18 = MLE_Fit("Normal", Y18);
run MLE_Summary(L18, 1, 1, 0.10);  /* alpha=0.10 */


/* ========================================
   TEST 19: Custom Title Label
   Expected: Custom title in output
   ======================================== */
Y19 = {2, 3, 4, 5, 6, 7};
L19 = MLE_Fit("Expo", Y19);
run MLE_Summary(L19, 1, 1, 0.05, "Exponential Fit (Custom Title)");


/* ========================================
   TEST 20: Single Parameter Distribution with CI
   Expected: CI for single parameter (Exponential)
   ======================================== */
Y20 = {1.2, 0.8, 0.9, 1.5, 1.1};
L20 = MLE_Fit("Expo", Y20);
run MLE_Summary(L20, 1, 1);


/* ========================================
   TEST 21: showZ=0, showP=0 Options
   Expected: No Z-statistics or p-values in table
   ======================================== */
Y21 = j(100, 1, .);
call randgen(Y21, "Normal", 5, 2);
L21 = MLE_Fit("Normal", Y21);
run MLE_Summary(L21, 1, 0, 0.05, "", 0, 0);  /* showZ=0, showP=0 */


/* ========================================
   TEST 22: All Options Combined - printOpt=2 + CI
   Expected: Extended output with all features
   ======================================== */
Y22 = j(200, 1, .);
call randgen(Y22, "Gamma", 3, 2);
L22 = MLE_Fit("Gamma", Y22);
run MLE_Summary(L22, 2, 1, 0.05, "Comprehensive Gamma Analysis");


QUIT;
