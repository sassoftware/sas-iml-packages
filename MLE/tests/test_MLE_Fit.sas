/* BEFORE running this example, store the modules in the MLE package, as shown in test_Install.sas */

/* test calling the top-level subroutines */
PROC IML;
load module=_ALL_;

use sashelp.heart;
   read all var "Systolic";
close;

L_gamma = MLE_Fit("Gamma", Systolic);

run MLE_Summary(L_gamma);     /* print basic results */
run MLE_Summary(L_gamma, 2);  /* print all results */
QUIT;

PROC IML;
load module=_ALL_;
call randseed(98765);

/*******************************************/
/* TEST 1: Normal Distribution - Basic Functionality */
/*******************************************/
title "TEST 1: Normal Distribution - Basic Functionality";
N1 = 1000;
mu = 10;
sigma = 5;
Y1 = j(N1, 1, .);
call randgen(Y1, "Normal", mu, sigma);
L1 = MLE_Fit("Normal", Y1);
run MLE_Summary(L1) printOpt=2 showCI=1;

/* Verify FitObj structure */
ListNames = ListGetName(L1);
reset wide;
print ListNames[label="Components of Fit Object"];

/*******************************************/
/* TEST 2: Different Distributions */
/*******************************************/
/* 2a. Exponential Distribution */
title "--- 2a. Exponential Distribution ---";
print /;
N2a = 800;
sigma = 3.5;
Y2a = j(N2a, 1, .);
call randgen(Y2a, "Expo", sigma);

L2a = MLE_Fit("Expo", Y2a);
run MLE_Summary(L2a);

/* 2b. Gamma Distribution */
title "--- 2b. Gamma Distribution ---";
print /;
N2b = 600;
alpha = 4;
lambda = 2;
Y2b = j(N2b, 1, .);
call randgen(Y2b, "Gamma", alpha, lambda);
L2b = MLE_Fit("Gamma", Y2b);
run MLE_Summary(L2b);

/* 2c. Beta Distribution */
title "--- 2c. Beta Distribution ---";
print /;
N2c = 500;
alpha = 2.5;
beta = 4.5;
Y2c = j(N2c, 1, .);
call randgen(Y2c, "Beta", alpha, beta);
L2c = MLE_Fit("Beta", Y2c);
run MLE_Summary(L2c);

/* 2d. Weibull Distribution */
title "--- 2d. Weibull Distribution ---";
print /;
N2d = 700;
c = 2.2;
lambda = 6;
Y2d = j(N2d, 1, .);
call randgen(Y2d, "Weibull", c, lambda);
L2d = MLE_Fit("Weibull", Y2d);
run MLE_Summary(L2d);


/*******************************************/
/* TEST 3: Custom Initial Values */
/*******************************************/
/* 3a. Default initialization (MoM) */
title "--- 3a. Default MoM Initialization ---";
print /;
N3 = 500;
Y3 = j(N3, 1, .);
call randgen(Y3, "Normal", 5, 2);
L3a = MLE_Fit("Normal", Y3);
run MLE_Summary(L3a);

/* 3b. Good custom initial guess */
title "--- 3b. Good Custom Initial Guess ---";
print /;
init_good = {4.9, 2.1};
L3b = MLE_Fit("Normal", Y3, init_good);
run MLE_Summary(L3b);

/* 3c. Poor custom initial guess */
title "--- 3c. Poor Custom Initial Guess ---";
print /;
init_poor = {50, 20};
L3c = MLE_Fit("Normal", Y3, init_poor);
run MLE_Summary(L3c);
/* TEMP TEMP TEMP */
STOP;

/*******************************************/
/* TEST 4: Different Optimization Methods */
/*******************************************/
title "--- 4a. NLPQN Method (default) ---";
print /;
N4 = 800;
Y4 = j(N4, 1, .);
call randgen(Y4, "Normal", 10, 3);
L4a = MLE_Fit("Normal", Y4, , "NLPQN");
run MLE_Summary(L4a);

title "--- 4b. NLPNRA Method ---";
print /;
L4b = MLE_Fit("Normal", Y4, , "NLPNRA");
run MLE_Summary(L4b);

if ^IsSAS9() then do;
   /* 4c. ACTIVE (Viya only) */
   title "--- 4c. ACTIVE Method ---";
   print /; 
   L4c = MLE_Fit("Normal", Y4, , "ACTIVE");
   run MLE_Summary(L4c);
end;

/*******************************************/
/* TEST 5: Missing Data Handling */
/*******************************************/
/* Data with missing values */
title "--- 5a. Data (n=15) with 3 missing values ---";
print /; 
Y5 = {5.2, 4.8, ., 6.1, 5.5, 4.9, ., 5.8, 6.2, 5.1, ., 5.4, 4.7, 6.0, 5.3};
L5 = MLE_Fit("Normal", Y5);
Y5_clean = L5$"y";
est5 = L5$"Estimate";
print (rowvec(Y5_clean))[format=4.1 L="Clean data (n=12)" c=(1:12)] ;
print est5[colname={"mu" "sigma"} format=10.6 L="Estimates"];

/*******************************************/
/* TEST 6: Hessian and Covariance Matrix */
/*******************************************/
title "--- TEST 6: Hessian and Covariance Matrix";
print /; 
N7 = 1000;
Y7 = j(N7, 1, .);
call randgen(Y7, "Normal", 5, 2);
L7 = MLE_Fit("Normal", Y7);
ods select  grad grad_norm hessian eigval;
run MLE_Summary(L7, 2);
ods select all;

/*******************************************/
/* TEST 7: Degenerate data */
/*******************************************/
title "--- TEST 7: Constant data";
print /; 
Y8 = repeat({9}, 20);
L8 = MLE_Fit("Normal", Y8);
run MLE_Summary(L8);

QUIT;
