proc iml;
%include "&MLE_path./MLE_define.sas";
QUIT;

/* test calls to the top-level subroutine */
PROC IML;
load module=_ALL_;

use sashelp.heart;
   read all var "Systolic";
close;

/* Primary use Case: call top-level MLE routine to get estimates */
/* parameter estimates */
gamma_est = MLE("Gamma", Systolic);  /* default guess is MoM */
print gamma_est;

/* you can get the MoM directly and use it (or another guess)
   est_MoM = MLE_MoM("Dist", y);
*/
gamma_Mom = MLE_MoM("Gamma", Systolic);  
print gamma_MoM;
gamma_est = MLE("Gamma", Systolic, gamma_Mom);  /* specify a guess */
print gamma_est;
QUIT;


proc iml;
load module=_ALL_;

/*******************************************/
/* MLE function tests */
/*******************************************/

call randseed(12345);
print "=== Testing MLE Function ===";

/*******************************************/
/* Test 1: Normal distribution */
/*******************************************/
print "Test 1a: Normal Distribution - Basic MLE";
N = 1000;
mu_true = 5;
sigma_true = 2;
Y_normal = j(N, 1, .);
call randgen(Y_normal, "Normal", mu_true, sigma_true);
/* Test MLE with default parameters (MoM initial guess, NLPQN method) */
est_normal = MLE("Normal", Y_normal);
res_normal = (mu_true || sigma_true) // est_normal;
print res_normal[rowname={"True" "MLE_Est"} colname={"mu" "sigma"}];

/*******************************************/
/* Customized initial guess - Normal Distribution */
/*******************************************/
print "Test 1b: Normal Distribution - Customized initial guess";

good_guess_normal = {4.8 2.1};  /* Close to true values */
est_normal_good = MLE("Normal", Y_normal, good_guess_normal);
poor_guess_normal = {100 50};   /* Far from true values */
est_normal_poor = MLE("Normal", Y_normal, poor_guess_normal);
res_normal_custom =  (good_guess_normal || est_normal_good) // (poor_guess_normal || est_normal_poor);
print res_normal_custom[rowname={"Good Initial Guess" "Poor Initial Guess" "MLE_Est"} colname={"mu_guess" "sigma_guess" "mu_est" "sigma_est"}];

   /*******************************************/
/* Different Optimization Methods */
   /*******************************************/
print "Test 1c: Normal Distribution - Different Optimization Methods";
   
/* Test NLP method */
   est_normal_NLPQN = MLE("Normal", Y_normal, , "NLPQN");
   est_normal_NLPNRA = MLE("Normal", Y_normal, , "NLPNRA");
   
/* Test NLPSOLVE method (Viya only) */
   est_normal_ACTIVE = MLE("Normal", Y_normal, , "ACTIVE");
est_normal_IP = MLE("Normal", Y_normal, , "IP");
est_normal_IPDIRECT = MLE("Normal", Y_normal, , "IPDIRECT");

res_normal_methods = (est_normal_NLPQN // est_normal_NLPNRA // est_normal_ACTIVE // est_normal_IP // est_normal_IPDIRECT);
print res_normal_methods[rowname={"NLPQN" "NLPNRA" "ACTIVE" "IP" "IPDIRECT"} colname={"mu" "sigma"}];


/*******************************************/
/* Test 2: Beta Distribution */
/*******************************************/
print " ";
print "========================================";
print "Test 2: Beta Distribution (lik_LL_Beta)";
print "========================================";

N2 = 800;
alpha_true_beta = 3;
beta_true_beta = 2;
Y_beta = j(N2, 1, .);
call randgen(Y_beta, "Beta", alpha_true_beta, beta_true_beta);

print "True parameters: alpha=3, beta=2";
print "Data range: [0, 1] (bounded)";
print " ";

/* Test 2a: Default MLE */
print "--- Test 2a: Default MLE ---";
est_beta = MLE("Beta", Y_beta);
res_beta = (alpha_true_beta || beta_true_beta) // est_beta;
print res_beta[rowname={"True" "MLE_Est"} colname={"alpha" "beta"} format=8.4];

/* Test 2b: Custom initial values */
print " ";
print "--- Test 2b: Custom Initial Values ---";
init_beta_good = {2.8, 1.9};
init_beta_poor = {10, 10};
est_beta_good = MLE("Beta", Y_beta, init_beta_good);
est_beta_poor = MLE("Beta", Y_beta, init_beta_poor);
res_beta_init = (init_beta_good || est_beta_good) // (init_beta_poor || est_beta_poor);
print res_beta_init[rowname={"Good Init" "Poor Init"} 
                     colname={"alpha_init" "beta_init" "alpha_est" "beta_est"} format=8.4];

/* Test 2c: Different optimization methods */
print " ";
print "--- Test 2c: Different Optimization Methods ---";
est_beta_NLPQN = MLE("Beta", Y_beta, , "NLPQN");
est_beta_NLPNRA = MLE("Beta", Y_beta, , "NLPNRA");
est_beta_ACTIVE = MLE("Beta", Y_beta, , "ACTIVE");
est_beta_IP = MLE("Beta", Y_beta, , "IP");
est_beta_IPDIRECT = MLE("Beta", Y_beta, , "IPDIRECT");
res_beta_methods = (est_beta_NLPQN // est_beta_NLPNRA // est_beta_ACTIVE // est_beta_IP // est_beta_IPDIRECT);
print res_beta_methods[rowname={"NLPQN" "NLPNRA" "ACTIVE" "IP" "IPDIRECT"} colname={"alpha" "beta"} format=8.4];


/*******************************************/
/* Test 3: Exponential Distribution */
/*******************************************/
print " ";
print "========================================";
print "Test 3: Exponential Distribution (lik_LL_Expo)";
print "========================================";

N3 = 1000;
scale_true = 4.5;
Y_expo = j(N3, 1, .);
call randgen(Y_expo, "Expo", scale_true);

print "True parameter: scale=4.5";
print "Data range: (0, infinity)";
print " ";

/* Test 3a: Default MLE */
print "--- Test 3a: Default MLE ---";
est_expo = MLE("Expo", Y_expo);
res_expo = scale_true // est_expo;
print res_expo[rowname={"True" "MLE_Est"} label="Scale Parameter" format=8.4];

/* Test 3b: Different scales */
print " ";
print "--- Test 3b: Small and Large Scales ---";
scale_small = 0.5;
scale_large = 50;
Y_expo_small = j(600, 1, .);
Y_expo_large = j(600, 1, .);
call randgen(Y_expo_small, "Expo", scale_small);
call randgen(Y_expo_large, "Expo", scale_large);

est_expo_small = MLE("Expo", Y_expo_small);
est_expo_large = MLE("Expo", Y_expo_large);
res_expo_scale = (scale_small || est_expo_small) // (scale_large || est_expo_large);
print res_expo_scale[rowname={"Small Scale" "Large Scale"} 
                      colname={"True" "Estimated"} format=8.4];

/* Test 3c: Different optimization methods */
print " ";
print "--- Test 3c: Different Optimization Methods ---";
est_expo_NLPQN = MLE("Expo", Y_expo, , "NLPQN");
est_expo_NLPNRA = MLE("Expo", Y_expo, , "NLPNRA");
est_expo_ACTIVE = MLE("Expo", Y_expo, , "ACTIVE");
est_expo_IP = MLE("Expo", Y_expo, , "IP");
est_expo_IPDIRECT = MLE("Expo", Y_expo, , "IPDIRECT");
res_expo_methods = (est_expo_NLPQN // est_expo_NLPNRA // est_expo_ACTIVE // est_expo_IP // est_expo_IPDIRECT);
print res_expo_methods[rowname={"NLPQN" "NLPNRA" "ACTIVE" "IP" "IPDIRECT"} label="Scale Estimates" format=8.4];


/*******************************************/
/* Test 4: Gamma Distribution */
/*******************************************/
print " ";
print "========================================";
print "Test 4: Gamma Distribution (lik_LL_Gamma)";
print "========================================";

N4 = 900;
alpha_true_gamma = 3.5;
lambda_true_gamma = 2.5;
Y_gamma = j(N4, 1, .);
call randgen(Y_gamma, "Gamma", alpha_true_gamma, lambda_true_gamma);

print "True parameters: alpha=3.5, lambda=2.5";
print "Mean = alpha*lambda = 8.75";
print " ";

/* Test 4a: Default MLE */
print "--- Test 4a: Default MLE ---";
est_gamma = MLE("Gamma", Y_gamma);
res_gamma = (alpha_true_gamma || lambda_true_gamma) // est_gamma;
print res_gamma[rowname={"True" "MLE_Est"} colname={"alpha" "lambda"} format=8.4];

/* Test 4b: Different shape parameters */
print " ";
print "--- Test 4b: Different Shape Parameters ---";
/* Small alpha (more skewed) */
alpha_small = 1.5;
lambda_gamma2 = 2;
Y_gamma_small = j(600, 1, .);
call randgen(Y_gamma_small, "Gamma", alpha_small, lambda_gamma2);
est_gamma_small = MLE("Gamma", Y_gamma_small);

/* Large alpha (more symmetric) */
alpha_large = 10;
Y_gamma_large = j(600, 1, .);
call randgen(Y_gamma_large, "Gamma", alpha_large, lambda_gamma2);
est_gamma_large = MLE("Gamma", Y_gamma_large);

res_gamma_shape = (alpha_small || lambda_gamma2 || est_gamma_small) // 
                   (alpha_large || lambda_gamma2 || est_gamma_large);
print res_gamma_shape[rowname={"Small Alpha (skewed)" "Large Alpha (symmetric)"} 
                       colname={"alpha_true" "lambda_true" "alpha_est" "lambda_est"} format=8.4];

/* Test 4c: Different optimization methods */
print " ";
print "--- Test 4c: Different Optimization Methods ---";
est_gamma_NLPQN = MLE("Gamma", Y_gamma, , "NLPQN");
est_gamma_NLPNRA = MLE("Gamma", Y_gamma, , "NLPNRA");
est_gamma_ACTIVE = MLE("Gamma", Y_gamma, , "ACTIVE");
est_gamma_IP = MLE("Gamma", Y_gamma, , "IP");
est_gamma_IPDIRECT = MLE("Gamma", Y_gamma, , "IPDIRECT");
res_gamma_methods = (est_gamma_NLPQN // est_gamma_NLPNRA // est_gamma_ACTIVE // est_gamma_IP // est_gamma_IPDIRECT);
print res_gamma_methods[rowname={"NLPQN" "NLPNRA" "ACTIVE" "IP" "IPDIRECT"} colname={"alpha" "lambda"} format=8.4];


/*******************************************/
/* Test 5: Gumbel Distribution */
/*******************************************/
print " ";
print "========================================";
print "Test 5: Gumbel Distribution (lik_LL_Gumbel)";
print "========================================";

N5 = 800;
mu_true_gumbel = 2;
sigma_true_gumbel = 1.5;
Y_gumbel = j(N5, 1, .);
call randgen(Y_gumbel, "Gumbel", mu_true_gumbel, sigma_true_gumbel);

print "True parameters: mu=2, sigma=1.5";
print "Distribution: Extreme value Type I";
print " ";

/* Test 5a: Default MLE */
print "--- Test 5a: Default MLE ---";
est_gumbel = MLE("Gumbel", Y_gumbel);
res_gumbel = (mu_true_gumbel || sigma_true_gumbel) // est_gumbel;
print res_gumbel[rowname={"True" "MLE_Est"} colname={"mu" "sigma"} format=8.4];

/* Test 5b: Different location parameters */
print " ";
print "--- Test 5b: Different Locations ---";
mu_negative = -5;
mu_positive = 10;
sigma_gumbel2 = 2;

Y_gumbel_neg = j(600, 1, .);
Y_gumbel_pos = j(600, 1, .);
call randgen(Y_gumbel_neg, "Gumbel", mu_negative, sigma_gumbel2);
call randgen(Y_gumbel_pos, "Gumbel", mu_positive, sigma_gumbel2);

est_gumbel_neg = MLE("Gumbel", Y_gumbel_neg);
est_gumbel_pos = MLE("Gumbel", Y_gumbel_pos);

res_gumbel_loc = (mu_negative || sigma_gumbel2 || est_gumbel_neg) // 
                  (mu_positive || sigma_gumbel2 || est_gumbel_pos);
print res_gumbel_loc[rowname={"Negative mu" "Positive mu"} 
                      colname={"mu_true" "sigma_true" "mu_est" "sigma_est"} format=8.4];

/* Test 5c: Different optimization methods */
print " ";
print "--- Test 5c: Different Optimization Methods ---";
est_gumbel_NLPQN = MLE("Gumbel", Y_gumbel, , "NLPQN");
est_gumbel_NLPNRA = MLE("Gumbel", Y_gumbel, , "NLPNRA");
est_gumbel_ACTIVE = MLE("Gumbel", Y_gumbel, , "ACTIVE");
est_gumbel_IP = MLE("Gumbel", Y_gumbel, , "IP");
est_gumbel_IPDIRECT = MLE("Gumbel", Y_gumbel, , "IPDIRECT");
res_gumbel_methods = (est_gumbel_NLPQN // est_gumbel_NLPNRA // est_gumbel_ACTIVE // est_gumbel_IP // est_gumbel_IPDIRECT);
print res_gumbel_methods[rowname={"NLPQN" "NLPNRA" "ACTIVE" "IP" "IPDIRECT"} colname={"mu" "sigma"} format=8.4];


/*******************************************/
/* Test 6: Inverse Gaussian Distribution */
/*******************************************/
print " ";
print "========================================";
print "Test 6: Inverse Gaussian Distribution (lik_LL_IGauss)";
print "========================================";

N6 = 700;
lambda_true_ig = 3;
mu_true_ig = 2;
Y_igauss = j(N6, 1, .);
call randgen(Y_igauss, "IGauss", lambda_true_ig, mu_true_ig);

print "True parameters: lambda=3, mu=2";
print "Also known as: Wald distribution";
print " ";

/* Test 6a: Default MLE */
print "--- Test 6a: Default MLE ---";
est_igauss = MLE("IGauss", Y_igauss);
res_igauss = (lambda_true_ig || mu_true_ig) // est_igauss;
print res_igauss[rowname={"True" "MLE_Est"} colname={"lambda" "mu"} format=8.4];

/* Test 6b: Alternative parameterizations */
print " ";
print "--- Test 6b: Different Parameter Values ---";
lambda_ig2 = 5;
mu_ig2 = 4;
Y_igauss2 = j(600, 1, .);
call randgen(Y_igauss2, "IGauss", lambda_ig2, mu_ig2);

est_igauss2 = MLE("IGauss", Y_igauss2);
res_igauss2 = (lambda_ig2 || mu_ig2) // est_igauss2;
print res_igauss2[rowname={"True" "MLE_Est"} colname={"lambda" "mu"} format=8.4];

/* Test 6c: Different optimization methods */
print " ";
print "--- Test 6c: Different Optimization Methods ---";
est_igauss_NLPQN = MLE("IGauss", Y_igauss, , "NLPQN");
est_igauss_NLPNRA = MLE("IGauss", Y_igauss, , "NLPNRA");
est_igauss_ACTIVE = MLE("IGauss", Y_igauss, , "ACTIVE");
est_igauss_IP = MLE("IGauss", Y_igauss, , "IP");
est_igauss_IPDIRECT = MLE("IGauss", Y_igauss, , "IPDIRECT");
res_igauss_methods = (est_igauss_NLPQN // est_igauss_NLPNRA // est_igauss_ACTIVE // est_igauss_IP // est_igauss_IPDIRECT);
print res_igauss_methods[rowname={"NLPQN" "NLPNRA" "ACTIVE" "IP" "IPDIRECT"} colname={"lambda" "mu"} format=8.4];


/*******************************************/
/* Test 7: Lognormal Distribution */
/*******************************************/
print " ";
print "========================================";
print "Test 7: Lognormal Distribution (lik_LL_LN2)";
print "========================================";

N7 = 1000;
mu_true_ln = 1.5;
sigma_true_ln = 0.8;
Y_lognormal = j(N7, 1, .);
call randgen(Y_lognormal, "Lognormal", mu_true_ln, sigma_true_ln);

print "True parameters: mu=1.5, sigma=0.8";
print "Median = exp(mu) = 4.48";
print " ";

/* Test 7a: Default MLE */
print "--- Test 7a: Default MLE ---";
est_ln = MLE("Lognormal", Y_lognormal);
res_ln = (mu_true_ln || sigma_true_ln) // est_ln;
print res_ln[rowname={"True" "MLE_Est"} colname={"mu" "sigma"} format=8.4];

/* Test 7b: Different dispersion levels */
print " ";
print "--- Test 7b: Low vs High Dispersion ---";
mu_ln2 = 2;
sigma_ln_low = 0.3;   /* Low dispersion */
sigma_ln_high = 1.5;  /* High dispersion */

Y_ln_low = j(600, 1, .);
Y_ln_high = j(600, 1, .);
call randgen(Y_ln_low, "Lognormal", mu_ln2, sigma_ln_low);
call randgen(Y_ln_high, "Lognormal", mu_ln2, sigma_ln_high);

est_ln_low = MLE("Lognormal", Y_ln_low);
est_ln_high = MLE("Lognormal", Y_ln_high);

res_ln_disp = (mu_ln2 || sigma_ln_low || est_ln_low) // 
               (mu_ln2 || sigma_ln_high || est_ln_high);
print res_ln_disp[rowname={"Low Dispersion" "High Dispersion"} 
                   colname={"mu_true" "sigma_true" "mu_est" "sigma_est"} format=8.4];

/* Test 7c: Different optimization methods */
print " ";
print "--- Test 7c: Different Optimization Methods ---";
est_ln_NLPQN = MLE("Lognormal", Y_lognormal, , "NLPQN");
est_ln_NLPNRA = MLE("Lognormal", Y_lognormal, , "NLPNRA");
est_ln_ACTIVE = MLE("Lognormal", Y_lognormal, , "ACTIVE");
est_ln_IP = MLE("Lognormal", Y_lognormal, , "IP");
est_ln_IPDIRECT = MLE("Lognormal", Y_lognormal, , "IPDIRECT");
res_ln_methods = (est_ln_NLPQN // est_ln_NLPNRA // est_ln_ACTIVE // est_ln_IP // est_ln_IPDIRECT);
print res_ln_methods[rowname={"NLPQN" "NLPNRA" "ACTIVE" "IP" "IPDIRECT"} colname={"mu" "sigma"} format=8.4];


/*******************************************/
/* Test 8: Weibull Distribution */
/*******************************************/
print " ";
print "========================================";
print "Test 8: Weibull Distribution (lik_LL_Weib2)";
print "========================================";

N8 = 850;
c_true = 2.5;
lambda_true_weib = 5;
Y_weibull = j(N8, 1, .);
call randgen(Y_weibull, "Weibull", c_true, lambda_true_weib);

print "True parameters: c=2.5, lambda=5";
print "Used in: Reliability analysis, survival data";
print " ";

/* Test 8a: Default MLE */
print "--- Test 8a: Default MLE ---";
est_weibull = MLE("Weibull", Y_weibull);
res_weibull = (c_true || lambda_true_weib) // est_weibull;
print res_weibull[rowname={"True" "MLE_Est"} colname={"c" "lambda"} format=8.4];

/* Test 8b: Different shape parameters */
print " ";
print "--- Test 8b: Different Shape Parameters ---";
/* c < 1: Decreasing hazard */
c_decreasing = 0.8;
lambda_weib2 = 4;
Y_weib_dec = j(600, 1, .);
call randgen(Y_weib_dec, "Weibull", c_decreasing, lambda_weib2);
est_weib_dec = MLE("Weibull", Y_weib_dec);

/* c = 1: Constant hazard (Exponential) */
c_constant = 1.0;
Y_weib_const = j(600, 1, .);
call randgen(Y_weib_const, "Weibull", c_constant, lambda_weib2);
est_weib_const = MLE("Weibull", Y_weib_const);

/* c > 1: Increasing hazard */
c_increasing = 3.5;
Y_weib_inc = j(600, 1, .);
call randgen(Y_weib_inc, "Weibull", c_increasing, lambda_weib2);
est_weib_inc = MLE("Weibull", Y_weib_inc);

res_weib_shape = (c_decreasing || lambda_weib2 || est_weib_dec) // 
                  (c_constant || lambda_weib2 || est_weib_const) //
                  (c_increasing || lambda_weib2 || est_weib_inc);
print res_weib_shape[rowname={"Decreasing Hazard (c<1)" "Constant Hazard (c=1)" "Increasing Hazard (c>1)"} 
                      colname={"c_true" "lambda_true" "c_est" "lambda_est"} format=8.4];

/* Test 8c: Different optimization methods */
print " ";
print "--- Test 8c: Different Optimization Methods ---";
est_weibull_NLPQN = MLE("Weibull", Y_weibull, , "NLPQN");
est_weibull_NLPNRA = MLE("Weibull", Y_weibull, , "NLPNRA");
est_weibull_ACTIVE = MLE("Weibull", Y_weibull, , "ACTIVE");
est_weibull_IP = MLE("Weibull", Y_weibull, , "IP");
est_weibull_IPDIRECT = MLE("Weibull", Y_weibull, , "IPDIRECT");
res_weibull_methods = (est_weibull_NLPQN // est_weibull_NLPNRA // est_weibull_ACTIVE // est_weibull_IP // est_weibull_IPDIRECT);
print res_weibull_methods[rowname={"NLPQN" "NLPNRA" "ACTIVE" "IP" "IPDIRECT"} colname={"c" "lambda"} format=8.4];

QUIT;

