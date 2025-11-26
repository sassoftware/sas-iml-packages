proc iml;
%include "&MLE_path./MLE_define.sas";
QUIT;

/* test calls to the top-level subroutine */
PROC IML;
load module=_ALL_;
print "--- MLE: A successfult test prints only 'TEST DONE' ---";

use sashelp.heart;
   read all var "Systolic";
close;

/* Primary use Case: call top-level MLE routine to get estimates */
/* parameter estimates */
gamma_est = MLE("Gamma", Systolic);  /* default guess is MoM */
correct = {
37.08868 ,
3.6914116 
};
if max(abs(gamma_est-correct)) > 1e-4  then 
   print "ERROR: Gamma results are incorrect";

/* you can get the MoM directly and use it (or another guess)
   est_MoM = MLE_MoM("Dist", y);
*/
gamma_MoM = MLE_MoM("Gamma", Systolic);  
correct = {
33.259906 ,
4.116355 
};
if max(abs(gamma_MoM-correct)) > 1e-4  then 
   print "ERROR: Gamma MoM results are incorrect";

gamma_est2 = MLE("Gamma", Systolic, gamma_Mom);  /* specify a guess */
correct = gamma_est;
if max(abs(gamma_est2-correct)) > 1e-4  then 
   print "ERROR: Gamma estimates are incorrect";
   
print "TEST DONE";
QUIT;


/*******************************************/
/* MLE function tests */
/*******************************************/

proc iml;
load module=_ALL_;
call randseed(12345, 1);

printFlag = 0;
if printFlag then 
   print "--- A successfult test prints multiple tables and the text 'TEST DONE' ---";
else
   print "--- A successfult test prints only 'TEST DONE' ---";

/*******************************************/
/* Normal distribution */
/*******************************************/
title "Test: Normal Distribution - Basic MLE";
N = 1000;
mu = 5;
sigma = 2;
Y = j(N, 1, .);
call randgen(Y, "Normal", mu, sigma);
/* Test MLE with default parameters (MoM initial guess, NLPQN method) */
est = MLE("Normal", Y);
result = (mu // sigma) || est;
if printFlag then 
   print result[c={"Parameter" "MLE_Est"} r={"mu" "sigma"} L="Basic MLE: Normal"];
correct = {
5 4.9788067 ,
2 2.0177194 
};
if max(abs(result-correct)) > 1e-4 then 
   print "ERROR: Basic MLE: Normal";

/*******************************************/
/* Customized initial guess - Normal Distribution */
/*******************************************/
title "Test: Normal Distribution - Customized initial guess";

guess = {4.8 2.1};  /* Close to Parameter values */
est = MLE("Normal", Y, guess);
result =  colvec(guess) || est;
if printFlag then 
   print result[c={"Initial Guess" "MLE_Est"} r={"mu" "sigma"} L="Close Custom Guess: Normal"];
correct = {
4.8 4.9788049 ,
2.1 2.0177162 
};
if max(abs(result-correct)) > 1e-4 then 
   print "ERROR: Close Custom Guess: Normal";

guess = {40 50};   /* Far from Parameter values */
est = MLE("Normal", Y, guess);
result =  colvec(guess) || est;
if printFlag then 
   print result[c={"Initial Guess" "MLE_Est"} r={"mu" "sigma"} L="Far-away Custom Guess: Normal"];
correct = {
40 4.9788067 ,
50 2.0177163 
};
if max(abs(result-correct)) > 1e-4 then 
   print "ERROR: Far-away Custom Guess: Normal";

/*******************************************/
/* Different Optimization Methods */
/*******************************************/
title "Test: Normal Distribution - Different Optimization Methods";

/* Test NLPQN and NLPNRA methods */
est_NLPQN = MLE("Normal", Y, , "NLPQN");
est_NLPNRA = MLE("Normal", Y, , "NLPNRA");

/* Test NLPSOLVE methods: ACTIVE, IP, and IPDIRECT (Viya only) */
%if (&sysver. ne 9.4) %then %do;
   est_ACTIVE = MLE("Normal", Y, , "ACTIVE");
   est_IP = MLE("Normal", Y, , "IP");
   est_IPDIRECT = MLE("Normal", Y, , "IPDIRECT");
   result = (est_NLPQN || est_NLPNRA || est_ACTIVE || est_IP || est_IPDIRECT);
   if printFlag then 
      print result[c={"NLPQN" "NLPNRA" "ACTIVE" "IP" "IPDIRECT"} r={"mu" "sigma"} L="Optim Methods: Normal"];
   correct = {
        4.9788072 4.9788067 4.9788067 4.9788067 4.9788067 ,
        2.0177164 2.0177164 2.0177164 2.0177164 2.0177164
   };
   if max(abs(result-correct)) > 1e-4 then 
      print "ERROR: Optim Methods: Normal";
%end;

%else %do;
   result = (est_NLPQN || est_NLPNRA);
   if printFlag then 
      print result[c={"NLPQN" "NLPNRA"} r={"mu" "sigma"} L="Optim Methods: Normal"];
   correct = {
   4.9788072 4.9788067 ,
   2.0177164 2.0177164 
   };
   if max(abs(result-correct)) > 1e-4 then 
      print "ERROR: Optim Methods: Normal";
%end;

/*******************************************/
/* Beta Distribution */
/*******************************************/
title "Test: Beta Distribution - Default MLE";

N2 = 800;
alpha = 3;
beta = 2;
Y_beta = j(N2, 1, .);
call randgen(Y_beta, "Beta", alpha, beta);

est_beta = MLE("Beta", Y_beta);
res_beta = (alpha // beta) || est_beta;
if printFlag then 
   print res_beta[c={"Parameter" "MLE_Est"} r={"alpha" "beta"} L="Default MLE: Beta"];
correct = {
3 2.9234179 ,
2 1.9278731 
};
if max(abs(res_beta-correct)) > 1e-4 then 
   print "ERROR: Default MLE: Beta";

/* Custom initial values */
title "Test: Beta Distribution - Custom Initial Values";
guess_beta = {2.8, 1.9};
est_beta = MLE("Beta", Y_beta, guess_beta);
result =  colvec(guess_beta) || est_beta;
if printFlag then 
   print result[c={"Guess" "MLE_Est"} r={"alpha" "beta"} L="Close Custom Guess: Beta"];
correct = {
2.8 2.9234184 ,
1.9 1.9278734 
};
if max(abs(result-correct)) > 1e-4 then 
   print "ERROR: Close Custom Guess: Beta";

guess_beta = {10, 10};
est_beta = MLE("Beta", Y_beta, guess_beta);
result =  colvec(guess_beta) || est_beta;
if printFlag then 
   print result[c={"Guess" "MLE_Est"} r={"alpha" "beta"} L="Far-away Custom Guess: Beta"];
correct = {
10 2.9234156 ,
10 1.9278713 
};
if max(abs(result-correct)) > 1e-4 then 
   print "ERROR: Far-away Custom Guess: Beta";

/* Different optimization methods */
title "Test: Beta Distribution - Different Optimization Methods";
est_beta_NLPQN = MLE("Beta", Y_beta, , "NLPQN");
est_beta_NLPNRA = MLE("Beta", Y_beta, , "NLPNRA");

/* NLPSOLVE methods: ACTIVE, IP, and IPDIRECT (Viya only) */
%if (&sysver. ne 9.4) %then %do;
   est_beta_ACTIVE = MLE("Beta", Y_beta, , "ACTIVE");
   est_beta_IP = MLE("Beta", Y_beta, , "IP");
   est_beta_IPDIRECT = MLE("Beta", Y_beta, , "IPDIRECT");

   res_beta_methods = est_beta_NLPQN || est_beta_NLPNRA || est_beta_ACTIVE || est_beta_IP || est_beta_IPDIRECT;
   if printFlag then 
      print res_beta_methods[c={"NLPQN" "NLPNRA" "ACTIVE" "IP" "IPDIRECT"} r={"alpha" "beta"} L="Optim Methods: Beta"];
   correct = {
               2.9234184 2.9234184 2.9234184 2.9234184 2.9234184 ,
               1.9278735 1.9278735 1.9278735 1.9278735 1.9278735 
   };
   if max(abs(res_beta_methods-correct)) > 1e-4 then 
      print "ERROR: Optim Methods: Beta";
%end;

%else %do;
   res_beta_methods = (est_beta_NLPQN || est_beta_NLPNRA);
   if printFlag then 
      print res_beta_methods[rowname={"alpha" "beta"} colname={"NLPQN" "NLPNRA"}];
   correct = {
   2.9234184 2.9234184 ,
   1.9278735 1.9278735 
   };
   if max(abs(res_beta_methods-correct)) > 1e-4 then 
      print "ERROR: Optim Methods: Beta";
%end;

/*******************************************/
/* Exponential Distribution */
/*******************************************/

title "Test: Exponential Distribution - Default MLE";

N3 = 1000;
sigma = 4.5;
Y_expo = j(N3, 1, .);
call randgen(Y_expo, "Expo", sigma);

est_expo = MLE("Expo", Y_expo);
res_expo = sigma || est_expo;
if printFlag then 
   print res_expo[c={"Parameter" "MLE_Est"} r={"sigma"} L="Default MLE: Exponential"];
correct = {
4.5 4.3838043
};
if max(abs(res_expo-correct)) > 1e-4 then 
   print "ERROR: Default MLE: Exponential";

/* Different scales */
title "Test: Exponential Distribution - Small and Large Scales";
scale_small = 0.5;
Y_expo_small = j(600, 1, .);
call randgen(Y_expo_small, "Expo", scale_small);
est_expo_small = MLE("Expo", Y_expo_small);
res_expo_small = scale_small || est_expo_small;
if printFlag then 
   print res_expo_small[c={"Parameter" "MLE_Est"} r={"sigma"} L="Small Scale: Exponential"];
correct = {
0.5 0.4988719
};
if max(abs(res_expo_small-correct)) > 1e-4 then 
   print "ERROR: Small Scale: Exponential";

scale_large = 50;
Y_expo_large = j(600, 1, .);
call randgen(Y_expo_large, "Expo", scale_large);
est_expo_large = MLE("Expo", Y_expo_large);
res_expo_large = scale_large || est_expo_large;
if printFlag then 
   print res_expo_large[c={"Parameter" "MLE_Est"} r={"sigma"} L="Large Scale: Exponential"];
correct = {
50 50.017493
};
if max(abs(res_expo_large-correct)) > 1e-4 then 
   print "ERROR: Large Scale: Exponential";  


/* Different optimization methods */
title "Test: Exponential Distribution - Different Optimization Methods";
est_expo_NLPQN = MLE("Expo", Y_expo, , "NLPQN");
est_expo_NLPNRA = MLE("Expo", Y_expo, , "NLPNRA");
%if (&sysver. ne 9.4) %then %do;
   est_expo_ACTIVE = MLE("Expo", Y_expo, , "ACTIVE");
   est_expo_IP = MLE("Expo", Y_expo, , "IP");
   est_expo_IPDIRECT = MLE("Expo", Y_expo, , "IPDIRECT");
   res_expo_methods = est_expo_NLPQN || est_expo_NLPNRA || est_expo_ACTIVE || est_expo_IP || est_expo_IPDIRECT;
   if printFlag then 
      print res_expo_methods[c={"NLPQN" "NLPNRA" "ACTIVE" "IP" "IPDIRECT"} r={"sigma"} L="Optim Methods: Exponential"];
   correct = {
   4.3838043 4.3838043 4.3838043 4.3838043 4.3838043
   };
   if max(abs(res_expo_methods-correct)) > 1e-4 then 
      print "ERROR: Optim Methods: Exponential";
%end; 
%else %do;
   res_expo_methods = (est_expo_NLPQN || est_expo_NLPNRA);
   if printFlag then 
      print res_expo_methods[c={"NLPQN" "NLPNRA"} r={"sigma"} L="Optim Methods: Exponential"];
   correct = {
   4.3838043 4.3838043
   };
   if max(abs(res_expo_methods-correct)) > 1e-4 then 
      print "ERROR: Optim Methods: Exponential";
%end;

/*******************************************/
/* Gamma Distribution */
/*******************************************/
title "Test: Gamma Distribution - Default MLE";

N4 = 900;
alpha = 3.5;
lambda = 2.5;
Y_gamma = j(N4, 1, .);
call randgen(Y_gamma, "Gamma", alpha, lambda);

est_gamma = MLE("Gamma", Y_gamma);
res_gamma = (alpha // lambda) || est_gamma;
if printFlag then 
   print res_gamma[c={"Parameter" "MLE_Est"} r={"alpha" "lambda"} L="Default MLE: Gamma"];
correct = {
3.5 3.4019092 ,
2.5 2.5537189
};
if max(abs(res_gamma-correct)) > 1e-4 then 
   print "ERROR: Default MLE: Gamma";

/* Shape parameters */
title "Test: Gamma Distribution - Different Shape Parameters";
alpha_small = 1.5;
lambda_gamma2 = 2;
Y_gamma_small = j(600, 1, .);
call randgen(Y_gamma_small, "Gamma", alpha_small, lambda_gamma2);
est_gamma_small = MLE("Gamma", Y_gamma_small);
res_gamma_small = (alpha_small // lambda_gamma2) || est_gamma_small;
if printFlag then 
   print res_gamma_small[c={"Parameter" "MLE_Est"} r={"alpha" "lambda"} L="Small Alpha: Gamma"];
correct = {
1.5 1.5911678 ,
2 1.9176423
};
if max(abs(res_gamma_small-correct)) > 1e-4 then 
   print "ERROR: Small Alpha: Gamma";

alpha_large = 10;
lambda_gamma2 = 2;
Y_gamma_large = j(600, 1, .);
call randgen(Y_gamma_large, "Gamma", alpha_large, lambda_gamma2);
est_gamma_large = MLE("Gamma", Y_gamma_large);
res_gamma_large = (alpha_large // lambda_gamma2) || est_gamma_large;
if printFlag then 
   print res_gamma_large[c={"Parameter" "MLE_Est"} r={"alpha" "lambda"} L="Large Alpha: Gamma"];
correct = {
10 11.399901 ,
2 1.7405274
};
if max(abs(res_gamma_large-correct)) > 1e-4 then 
   print "ERROR: Large Alpha: Gamma";


/* Different optimization methods */
title "Test: Gamma Distribution - Different Optimization Methods";
est_gamma_NLPQN = MLE("Gamma", Y_gamma, , "NLPQN");
est_gamma_NLPNRA = MLE("Gamma", Y_gamma, , "NLPNRA");
%if (&sysver. ne 9.4) %then %do;
   est_gamma_ACTIVE = MLE("Gamma", Y_gamma, , "ACTIVE");
   est_gamma_IP = MLE("Gamma", Y_gamma, , "IP");
   est_gamma_IPDIRECT = MLE("Gamma", Y_gamma, , "IPDIRECT");
   res_gamma_methods = est_gamma_NLPQN || est_gamma_NLPNRA || est_gamma_ACTIVE || est_gamma_IP || est_gamma_IPDIRECT;
   if printFlag then 
      print res_gamma_methods[c={"NLPQN" "NLPNRA" "ACTIVE" "IP" "IPDIRECT"} r={"alpha" "lambda"} L="Optim Methods: Gamma"];
   correct = {
      3.4019092 3.4019098 3.4019097 3.4019097 3.4019097 ,
      2.5537189 2.5537184 2.5537186 2.5537186 2.5537186
   };
   if max(abs(res_gamma_methods-correct)) > 1e-4 then 
      print "ERROR: Optim Methods: Gamma";
%end; 
%else %do;
   res_gamma_methods = (est_gamma_NLPQN || est_gamma_NLPNRA);
   if printFlag then 
      print res_gamma_methods[c={"NLPQN" "NLPNRA"} r={"alpha" "lambda"} L="Optim Methods: Gamma"];
   correct = {
      3.4019092 3.4019098 ,
      2.5537189 2.5537184
   };
   if max(abs(res_gamma_methods-correct)) > 1e-4 then 
      print "ERROR: Optim Methods: Gamma";
%end;



/*******************************************/
/* Inverse Gaussian Distribution */
/*******************************************/
title "Test: Inverse Gaussian Distribution - Default MLE";

N6 = 700;
lambda = 3;
mu = 2;
Y_igauss = j(N6, 1, .);
call randgen(Y_igauss, "IGauss", lambda, mu);

est_igauss = MLE("IGauss", Y_igauss);
res_igauss = (lambda // mu) || est_igauss;
if printFlag then 
   print res_igauss[c={"Parameter" "MLE_Est"} r={"lambda" "mu"} L="Default MLE: Inverse Gaussian"];
correct = {
   3 3.0986385 ,
   2 2.0489423
};
if max(abs(res_igauss-correct)) > 1e-4 then 
   print "ERROR: Default MLE: Inverse Gaussian";

/* Alternative parameterizations */
title "Test: Inverse Gaussian Distribution - Different Parameter Values";
lambda_ig2 = 5;
mu_ig2 = 4;
Y_igauss2 = j(600, 1, .);
call randgen(Y_igauss2, "IGauss", lambda_ig2, mu_ig2);

est_igauss2 = MLE("IGauss", Y_igauss2);
res_igauss2 = (lambda_ig2 // mu_ig2) || est_igauss2;
if printFlag then 
   print res_igauss2[c={"Parameter" "MLE_Est"} r={"lambda" "mu"} L="Different Parameter Values: Inverse Gaussian"];
correct = {
   5 4.9975411 ,
   4 3.8397393
};
if max(abs(res_igauss2-correct)) > 1e-4 then 
   print "ERROR: Different Parameter Values: Inverse Gaussian";

/* Different optimization methods */
title "Test: Inverse Gaussian Distribution - Different Optimization Methods";
est_igauss_NLPQN = MLE("IGauss", Y_igauss, , "NLPQN");
est_igauss_NLPNRA = MLE("IGauss", Y_igauss, , "NLPNRA");
%if (&sysver. ne 9.4) %then %do;
   est_igauss_ACTIVE = MLE("IGauss", Y_igauss, , "ACTIVE");
   est_igauss_IP = MLE("IGauss", Y_igauss, , "IP");
   est_igauss_IPDIRECT = MLE("IGauss", Y_igauss, , "IPDIRECT");
   res_igauss_methods = est_igauss_NLPQN || est_igauss_NLPNRA || est_igauss_ACTIVE || est_igauss_IP || est_igauss_IPDIRECT;
   if printFlag then 
      print res_igauss_methods[c={"NLPQN" "NLPNRA" "ACTIVE" "IP" "IPDIRECT"} r={"lambda" "mu"} L="Optim Methods: Inverse Gaussian"];
   correct = {
   3.0986385 3.0986393 3.0986393 3.0986393 3.0986393 ,
   2.0489423 2.0489423 2.0489423 2.0489423 2.0489423
   };
   if max(abs(res_igauss_methods-correct)) > 1e-4 then 
      print "ERROR: Optim Methods: Inverse Gaussian";
%end; 
%else %do;
   res_igauss_methods = (est_igauss_NLPQN || est_igauss_NLPNRA);
   if printFlag then 
      print res_igauss_methods[c={"NLPQN" "NLPNRA"} r={"lambda" "mu"} L="Optim Methods: Inverse Gaussian"];
   correct = {
   3.0986385    3.0986393 ,
   2.0489423    2.0489423
   };
   if max(abs(res_igauss_methods-correct)) > 1e-4 then 
      print "ERROR: Optim Methods: Inverse Gaussian";
%end;

/*******************************************/
* Lognormal Distribution */
/*******************************************/
title "Test: Lognormal Distribution - Default MLE";

N7 = 1000;
mu = 1.5;
sigma = 0.8;
Y_lognormal = j(N7, 1, .);
call randgen(Y_lognormal, "Lognormal", mu, sigma);

est_ln = MLE("Lognormal", Y_lognormal);
res_ln = (mu // sigma) || est_ln;
if printFlag then 
   print res_ln[c={"Parameter" "MLE_Est"} r={"mu" "sigma"} L="Default MLE: Lognormal"];
correct = {
   1.5 1.5004915 ,
   0.8 0.8154118
};
if max(abs(res_ln-correct)) > 1e-4 then 
   print "ERROR: Default MLE: Lognormal";

/* Different dispersion levels */
title "Test: Lognormal Distribution - Low vs High Dispersion";
mu_ln2 = 2;
sigma_ln_low = 0.3;   /* Low dispersion */
Y_ln_low = j(600, 1, .);
call randgen(Y_ln_low, "Lognormal", mu_ln2, sigma_ln_low);
est_ln_low = MLE("Lognormal", Y_ln_low);
res_ln_low = (mu_ln2 // sigma_ln_low) || est_ln_low;
if printFlag then 
   print res_ln_low[c={"Parameter" "MLE_Est"} r={"mu" "sigma"} L="Low Dispersion: Lognormal"];
correct = {
   2 2.0126728 ,
   0.3 0.3043297
};
if max(abs(res_ln_low-correct)) > 1e-4 then 
   print "ERROR: Low Dispersion: Lognormal";

sigma_ln_high = 1.5;  /* High dispersion */
Y_ln_high = j(600, 1, .);
call randgen(Y_ln_high, "Lognormal", mu_ln2, sigma_ln_high);
est_ln_high = MLE("Lognormal", Y_ln_high);
res_ln_high = (mu_ln2 // sigma_ln_high) || est_ln_high;
if printFlag then 
   print res_ln_high[c={"Parameter" "MLE_Est"} r={"mu" "sigma"} L="High Dispersion: Lognormal"];
correct = {
   2  2.0003243 ,
   1.5 1.447963
};
if max(abs(res_ln_high-correct)) > 1e-4 then 
   print "ERROR: High Dispersion: Lognormal";

/* Different optimization methods */
title "Test: Lognormal Distribution - Different Optimization Methods";
est_ln_NLPQN = MLE("Lognormal", Y_lognormal, , "NLPQN");
est_ln_NLPNRA = MLE("Lognormal", Y_lognormal, , "NLPNRA");
%if (&sysver. ne 9.4) %then %do;
   est_ln_ACTIVE = MLE("Lognormal", Y_lognormal, , "ACTIVE");
   est_ln_IP = MLE("Lognormal", Y_lognormal, , "IP");
   est_ln_IPDIRECT = MLE("Lognormal", Y_lognormal, , "IPDIRECT");
   res_ln_methods = est_ln_NLPQN || est_ln_NLPNRA || est_ln_ACTIVE || est_ln_IP || est_ln_IPDIRECT;
   if printFlag then 
      print res_ln_methods[c={"NLPQN" "NLPNRA" "ACTIVE" "IP" "IPDIRECT"} r={"mu" "sigma"} L="Optim Methods: Lognormal"];
   correct = {
   1.5004915 1.5004915 1.5004916 1.5004916 1.5004916 ,
   0.8154118 0.8154118 0.8154118 0.8154118 0.8154118
   };
   if max(abs(res_ln_methods-correct)) > 1e-4 then 
      print "ERROR: Optim Methods: Lognormal";
%end; 
%else %do;
   res_ln_methods = (est_ln_NLPQN || est_ln_NLPNRA);
   if printFlag then 
      print res_ln_methods[c={"NLPQN" "NLPNRA"} r={"mu" "sigma"} L="Optim Methods: Lognormal"];
   correct = {
   1.5004916 1.5004915 ,
   0.8154118 0.8154117
   };
   if max(abs(res_ln_methods-correct)) > 1e-4 then 
      print "ERROR: Optim Methods: Lognormal";
%end;

/*******************************************/
/* Weibull Distribution */
/*******************************************/
title "Test: Weibull Distribution - Default MLE";

N8 = 850;
c = 2.5;
lambda = 5;
Y_weibull = j(N8, 1, .);
call randgen(Y_weibull, "Weibull", c, lambda);

est_weibull = MLE("Weibull", Y_weibull);
res_weibull = (c // lambda) || est_weibull;
if printFlag then 
   print res_weibull[c={"Parameter" "MLE_Est"} r={"c" "lambda"} L="Default MLE: Weibull"];
correct = {
   2.5 2.4787196 ,
   5 4.8807104
};
if max(abs(res_weibull-correct)) > 1e-4 then 
   print "ERROR: Default MLE: Weibull";

/* Different shape parameters */
title "Test: Weibull Distribution - Different Shape Parameters";
/* c < 1: Decreasing hazard */
c_decreasing = 0.8;
lambda_weib2 = 4;
Y_weib_dec = j(600, 1, .);
call randgen(Y_weib_dec, "Weibull", c_decreasing, lambda_weib2);
est_weib_dec = MLE("Weibull", Y_weib_dec);
res_weib_dec = (c_decreasing // lambda_weib2) || est_weib_dec;
if printFlag then 
   print res_weib_dec[c={"Parameter" "MLE_Est"} r={"c" "lambda"} L="Decreasing Hazard: Weibull"];
correct = {
   0.8 0.7739147 ,
   4 4.1694258
};
if max(abs(res_weib_dec-correct)) > 1e-4 then 
   print "ERROR: Decreasing Hazard: Weibull";

/* c = 1: Constant hazard (Exponential) */
c_constant = 1.0;
Y_weib_const = j(600, 1, .);
call randgen(Y_weib_const, "Weibull", c_constant, lambda_weib2);
est_weib_const = MLE("Weibull", Y_weib_const);
res_weib_const = (c_constant // lambda_weib2) || est_weib_const;
if printFlag then 
   print res_weib_const[c={"Parameter" "MLE_Est"} r={"c" "lambda"} L="Constant Hazard: Weibull"];
correct = {
   1 0.9644618 ,
   4 3.8627185
};
if max(abs(res_weib_const-correct)) > 1e-4 then 
   print "ERROR: Constant Hazard: Weibull";

/* c > 1: Increasing hazard */
c_increasing = 3.5;
Y_weib_inc = j(600, 1, .);
call randgen(Y_weib_inc, "Weibull", c_increasing, lambda_weib2);
est_weib_inc = MLE("Weibull", Y_weib_inc);
res_weib_inc = (c_increasing // lambda_weib2) || est_weib_inc;
if printFlag then 
   print res_weib_inc[c={"Parameter" "MLE_Est"} r={"c" "lambda"} L="Increasing Hazard: Weibull"];
correct = {
   3.5 3.4951947 ,
   4 3.9559502
};
if max(abs(res_weib_inc-correct)) > 1e-4 then 
   print "ERROR: Increasing Hazard: Weibull";

/* Different optimization methods */
title "Test: Weibull Distribution - Different Optimization Methods";
est_weibull_NLPQN = MLE("Weibull", Y_weibull, , "NLPQN");
est_weibull_NLPNRA = MLE("Weibull", Y_weibull, , "NLPNRA");
%if (&sysver. ne 9.4) %then %do;
   est_weibull_ACTIVE = MLE("Weibull", Y_weibull, , "ACTIVE");
   est_weibull_IP = MLE("Weibull", Y_weibull, , "IP");
   est_weibull_IPDIRECT = MLE("Weibull", Y_weibull, , "IPDIRECT");
   res_weibull_methods = est_weibull_NLPQN || est_weibull_NLPNRA || est_weibull_ACTIVE || est_weibull_IP || est_weibull_IPDIRECT;
   if printFlag then 
      print res_weibull_methods[c={"NLPQN" "NLPNRA" "ACTIVE" "IP" "IPDIRECT"} r={"c" "lambda"} L="Optim Methods: Weibull"];
   correct = {
   2.4787196 2.4787196 2.4787196 2.4787196 2.4787196 ,
   4.8807104 4.8807103 4.8807104 4.8807104 4.8807104
   };
   if max(abs(res_weibull_methods-correct)) > 1e-4 then 
      print "ERROR: Optim Methods: Weibull";
%end; 
%else %do;
   res_weibull_methods = (est_weibull_NLPQN || est_weibull_NLPNRA);
   if printFlag then 
      print res_weibull_methods[c={"NLPQN" "NLPNRA"} r={"c" "lambda"} L="Optim Methods: Weibull"];
   correct = {
   2.4787197   2.4787195 ,
   4.8807105   4.8807104
   };
   if max(abs(res_weibull_methods-correct)) > 1e-4 then 
      print "ERROR: Optim Methods: Weibull";
%end;


/*******************************************/
/* Gumbel Distribution */
/*******************************************/
title "Test: Gumbel Distribution - Default MLE";

N5 = 800;
mu = 2;
sigma = 1.5;
Y_gumbel = j(N5, 1, .);
call randgen(Y_gumbel, "Gumbel", mu, sigma);

est_gumbel = MLE("Gumbel", Y_gumbel);
res_gumbel = (mu // sigma) || est_gumbel;
if printFlag then 
   print res_gumbel[c={"Parameter" "MLE_Est"} r={"mu" "sigma"} L="Default MLE: Gumbel"];
correct = {
   2 2.0732519 ,
   1.5 1.5377421
};
if max(abs(res_gumbel-correct)) > 1e-4 then 
   print "ERROR: Default MLE: Gumbel";

/* Different optimization methods */
title "Test: Gumbel Distribution - Different Optimization Methods";
est_gumbel_NLPQN = MLE("Gumbel", Y_gumbel, , "NLPQN");
est_gumbel_NLPNRA = MLE("Gumbel", Y_gumbel, , "NLPNRA");

%if (&sysver. ne 9.4) %then %do;
   est_gumbel_ACTIVE = MLE("Gumbel", Y_gumbel, , "ACTIVE");
   est_gumbel_IP = MLE("Gumbel", Y_gumbel, , "IP");
   est_gumbel_IPDIRECT = MLE("Gumbel", Y_gumbel, , "IPDIRECT");
   res_gumbel_methods = est_gumbel_NLPQN || est_gumbel_NLPNRA || est_gumbel_ACTIVE || est_gumbel_IP || est_gumbel_IPDIRECT;
   if printFlag then 
      print res_gumbel_methods[c={"NLPQN" "NLPNRA" "ACTIVE" "IP" "IPDIRECT"} r={"mu" "sigma"} L="Optim Methods: Gumbel"];
   correct = {
   2.0732519  2.073252  2.073252  2.073252  2.073252 ,
   1.5377421 1.5377421 1.5377421 1.5377421 1.5377421
   };
   if max(abs(res_gumbel_methods-correct)) > 1e-4 then 
      print "ERROR: Optim Methods: Gumbel";
%end; 
%else %do;
   res_gumbel_methods = (est_gumbel_NLPQN || est_gumbel_NLPNRA);
   if printFlag then 
      print res_gumbel_methods[c={"NLPQN" "NLPNRA"} r={"mu" "sigma"} L="Optim Methods: Gumbel"];
   correct = {
   2.0732522   2.0732518 ,
   1.5377422 1.537742
   };
   if max(abs(res_gumbel_methods-correct)) > 1e-4 then 
      print "ERROR: Optim Methods: Gumbel";
%end;

print "TEST DONE";

QUIT;
