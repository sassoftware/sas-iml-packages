proc iml;
%include "&MLE_path./MLE_define.sas";
QUIT;

/* test calls to the top-level subroutine */
PROC IML;
load module=_ALL_;

use sashelp.heart;
   read all var "Systolic";
close;

L_gamma = MLE_Fit("Gamma", Systolic);

run MLE_Summary(L_gamma);     /* print basic results */
run MLE_Summary(L_gamma, 2);  /* print all results */
QUIT;

proc iml;
%include "MLE_define.sas";
load module=_ALL_;

call randseed(12345);

print "========================================";
print "MLE_FIT FUNCTION TESTS";
print "========================================";
print " ";


/*******************************************/
/* TEST 1: Normal Distribution - Basic Functionality */
/*******************************************/
print "========================================";
print "TEST 1: Normal Distribution - Basic Functionality";
print "========================================";
N1 = 1000;
mu_true1 = 10;
sigma_true1 = 5;
Y1 = j(N1, 1, .);
call randgen(Y1, "Normal", mu_true1, sigma_true1);

print "True parameters: mu=10, sigma=5";
print "Sample size: 1000";
print " ";

L1 = MLE_Fit("Normal", Y1);

/* Verify FitObj structure */
print "--- FitObj Contents ---";
ListNames = ListGetName(L1);
print ListNames[label="List Components"];

/* Extract and display key results */
estimate1 = L1$"Estimate"; 
parmNames1 = L1$"ParmNames";
stdErr1 = L1$"StdErr";
LL1 = L1$"LL";
grad1 = L1$"Grad";
hess1 = L1$"Hessian";
crit1 = L1$"Crit";
critNames1 = L1$"CritNames";

print " ";
print "--- Parameter Estimates ---";
results1 = (mu_true1 || sigma_true1) // estimate1;
print results1[colname=parmNames1 rowname={"True" "Estimated"} format=10.6];

print " ";
print "--- Standard Errors ---";
print stdErr1[colname=parmNames1 format=10.6];

print " ";
print "--- Log-Likelihood ---";
print LL1[format=12.4];

print " ";
print "--- Gradient (should be near zero) ---";
print grad1[colname=parmNames1 format=12.8];

print " ";
print "--- Information Criteria ---";
print crit1[colname=critNames1 format=12.4];


/*******************************************/
/* TEST 2: Different Distributions */
/*******************************************/
print " ";
print "========================================";
print "TEST 2: Multiple Distributions";
print "========================================";

/* 2a. Exponential Distribution */
print " ";
print "--- 2a. Exponential Distribution ---";
N2a = 800;
scale_true = 3.5;
Y2a = j(N2a, 1, .);
call randgen(Y2a, "Expo", scale_true);

print "True parameter: scale=3.5";
L2a = MLE_Fit("Expo", Y2a);
est2a = L2a$"Estimate";
se2a = L2a$"StdErr";
print est2a[label="Estimate" format=10.6];
print se2a[label="Std Error" format=10.6];

/* 2b. Gamma Distribution */
print " ";
print "--- 2b. Gamma Distribution ---";
N2b = 600;
alpha_true = 4;
lambda_true = 2;
Y2b = j(N2b, 1, .);
call randgen(Y2b, "Gamma", alpha_true, lambda_true);

print "True parameters: alpha=4, lambda=2";
L2b = MLE_Fit("Gamma", Y2b);
est2b = L2b$"Estimate";
se2b = L2b$"StdErr";
parmNames2b = L2b$"ParmNames";
results2b = (alpha_true || lambda_true) // est2b;
print results2b[colname=parmNames2b rowname={"True" "Estimated"} format=10.6];

/* 2c. Beta Distribution */
print " ";
print "--- 2c. Beta Distribution ---";
N2c = 500;
alpha_beta = 2.5;
beta_beta = 4.5;
Y2c = j(N2c, 1, .);
call randgen(Y2c, "Beta", alpha_beta, beta_beta);

print "True parameters: alpha=2.5, beta=4.5";
L2c = MLE_Fit("Beta", Y2c);
est2c = L2c$"Estimate";
parmNames2c = L2c$"ParmNames";
results2c = (alpha_beta || beta_beta) // est2c;
print results2c[colname=parmNames2c rowname={"True" "Estimated"} format=10.6];

/* 2d. Weibull Distribution */
print " ";
print "--- 2d. Weibull Distribution ---";
N2d = 700;
c_true = 2.2;
lambda_weib = 6;
Y2d = j(N2d, 1, .);
call randgen(Y2d, "Weibull", c_true, lambda_weib);

print "True parameters: c=2.2, lambda=6";
L2d = MLE_Fit("Weibull", Y2d);
est2d = L2d$"Estimate";
parmNames2d = L2d$"ParmNames";
results2d = (c_true || lambda_weib) // est2d;
print results2d[colname=parmNames2d rowname={"True" "Estimated"} format=10.6];


/*******************************************/
/* TEST 3: Custom Initial Values */
/*******************************************/
print " ";
print "========================================";
print "TEST 3: Custom Initial Values";
print "========================================";

N3 = 500;
Y3 = j(N3, 1, .);
call randgen(Y3, "Normal", 5, 2);

/* 3a. Default initialization (MoM) */
print " ";
print "--- 3a. Default MoM Initialization ---";
L3a = MLE_Fit("Normal", Y3);
est3a = L3a$"Estimate";
print "True: mu=5, sigma=2";
print est3a[colname={"mu" "sigma"} format=10.6];

/* 3b. Good custom initial guess */
print " ";
print "--- 3b. Good Custom Initial Guess ---";
init_good = {4.9, 2.1};
L3b = MLE_Fit("Normal", Y3, init_good);
est3b = L3b$"Estimate";
print "Initial guess:" init_good[format=10.4];
print "Estimate:" est3b[format=10.6];

/* 3c. Poor custom initial guess */
print " ";
print "--- 3c. Poor Custom Initial Guess ---";
init_poor = {50, 20};
L3c = MLE_Fit("Normal", Y3, init_poor);
est3c = L3c$"Estimate";
print "Initial guess:" init_poor[format=10.4];
print "Estimate (should still converge):" est3c[format=10.6];


/*******************************************/
/* TEST 4: Different Optimization Methods */
/*******************************************/
print " ";
print "========================================";
print "TEST 4: Different Optimization Methods";
print "========================================";

N4 = 800;
Y4 = j(N4, 1, .);
call randgen(Y4, "Normal", 10, 3);

print "True parameters: mu=10, sigma=3";
print " ";

/* 4a. NLPQN (default) */
print "--- 4a. NLPQN Method ---";
L4a = MLE_Fit("Normal", Y4, , "NLPQN");
est4a = L4a$"Estimate";
print est4a[colname={"mu" "sigma"} format=10.6];

/* 4b. NLPNRA */
print " ";
print "--- 4b. NLPNRA Method ---";
L4b = MLE_Fit("Normal", Y4, , "NLPNRA");
est4b = L4b$"Estimate";
print est4b[colname={"mu" "sigma"} format=10.6];

/* 4c. ACTIVE (Viya only) */
print " ";
print "--- 4c. ACTIVE Method ---";
L4c = MLE_Fit("Normal", Y4, , "ACTIVE");
est4c = L4c$"Estimate";
print est4c[colname={"mu" "sigma"} format=10.6];

/* 4d. Comparison table */
print " ";
print "--- Comparison Across Methods ---";
comparison = est4a // est4b // est4c;
print comparison[colname={"mu" "sigma"} rowname={"NLPQN" "NLPNRA" "ACTIVE"} format=10.6];


/*******************************************/
/* TEST 5: Missing Data Handling */
/*******************************************/
print " ";
print "========================================";
print "TEST 5: Missing Data Handling";
print "========================================";

/* Data with missing values */
Y5 = {5.2, 4.8, ., 6.1, 5.5, 4.9, ., 5.8, 6.2, 5.1, ., 5.4, 4.7, 6.0, 5.3};
print "Original data (n=15, with 3 missing values):";
print Y5[format=5.2];

L5 = MLE_Fit("Normal", Y5);
Y5_clean = L5$"y";
est5 = L5$"Estimate";

print " ";
print "Clean data used in estimation (n=12):";
print Y5_clean[format=5.2];

print " ";
print "Estimates:";
print est5[colname={"mu" "sigma"} format=10.6];


/*******************************************/
/* TEST 6: Small vs Large Sample Sizes */
/*******************************************/
print " ";
print "========================================";
print "TEST 6: Small vs Large Sample Sizes";
print "========================================";

/* 6a. Very small sample */
print " ";
print "--- 6a. Small Sample (n=10) ---";
Y6a = {4.5, 5.2, 4.8, 5.5, 4.9, 5.1, 4.7, 5.3, 5.0, 4.6};
L6a = MLE_Fit("Normal", Y6a);
est6a = L6a$"Estimate";
se6a = L6a$"StdErr";
print est6a[colname={"mu" "sigma"} format=10.6 label="Estimates"];
print se6a[colname={"mu" "sigma"} format=10.6 label="Std Errors"];

/* 6b. Medium sample */
print " ";
print "--- 6b. Medium Sample (n=200) ---";
Y6b = j(200, 1, .);
call randgen(Y6b, "Normal", 5, 1);
L6b = MLE_Fit("Normal", Y6b);
est6b = L6b$"Estimate";
se6b = L6b$"StdErr";
print est6b[colname={"mu" "sigma"} format=10.6 label="Estimates"];
print se6b[colname={"mu" "sigma"} format=10.6 label="Std Errors"];

/* 6c. Large sample */
print " ";
print "--- 6c. Large Sample (n=5000) ---";
Y6c = j(5000, 1, .);
call randgen(Y6c, "Normal", 5, 1);
L6c = MLE_Fit("Normal", Y6c);
est6c = L6c$"Estimate";
se6c = L6c$"StdErr";
print est6c[colname={"mu" "sigma"} format=10.6 label="Estimates"];
print se6c[colname={"mu" "sigma"} format=10.6 label="Std Errors (note smaller SEs)"];

/* 6d. Standard error comparison */
print " ";
print "--- Standard Error Comparison ---";
se_comparison = se6a // se6b // se6c;
print se_comparison[colname={"StdErr(mu)" "StdErr(sigma)"} 
                     rowname={"n=10" "n=200" "n=5000"} format=10.6
                     label="Standard Errors Decrease with Sample Size"];


/*******************************************/
/* TEST 7: Hessian and Covariance Matrix */
/*******************************************/
print " ";
print "========================================";
print "TEST 7: Hessian and Covariance Matrix";
print "========================================";

N7 = 1000;
Y7 = j(N7, 1, .);
call randgen(Y7, "Normal", 5, 2);

L7 = MLE_Fit("Normal", Y7);
hess7 = L7$"Hessian";
se7 = L7$"StdErr";

print "--- Hessian Matrix ---";
print hess7[colname={"mu" "sigma"} rowname={"mu" "sigma"} format=12.4];

print " ";
print "--- Covariance Matrix (inv(-Hess)) ---";
cov7 = inv(-hess7);
print cov7[colname={"mu" "sigma"} rowname={"mu" "sigma"} format=12.6];

print " ";
print "--- Standard Errors ---";
print "From FitObj:";
print se7[colname={"mu" "sigma"} format=10.6];

se7_check = rowvec(sqrt(vecdiag(cov7)));
print " ";
print "From sqrt(diag(Cov)) (should match):";
print se7_check[colname={"mu" "sigma"} format=10.6];


/*******************************************/
/* TEST 8: Integration with MLE_Summary */
/*******************************************/
print " ";
print "========================================";
print "TEST 8: Integration with MLE_Summary";
print "========================================";

N8 = 800;
Y8 = j(N8, 1, .);
call randgen(Y8, "Lognormal", 1, 0.5);

print "Fitting Lognormal distribution";
print "True parameters: mu=1, sigma=0.5";
print " ";

L8 = MLE_Fit("Lognormal", Y8);

print "--- Summary Output (printOpt=2) ---";
run MLE_Summary(L8, 2);


QUIT;