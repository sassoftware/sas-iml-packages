/**************************************************/
/* Test: Compare MLE package with PROC NLMIXED    */
/**************************************************/

/* Generate sample data for multiple distributions */
%let nObs = 50;
data testdata;
   call streaminit(12345);
   do i = 1 to &nObs;
      beta_y = rand("Beta", 2, 5);           /* alpha=2, beta=5 */
      gamma_y = rand("Gamma", 3, 2);         /* shape=3, scale=2 */
      logn_y = rand("Lognormal", 1, 0.5);    /* mu=1, sigma=0.5 */
      expo_y = rand("Exponential", 4);       /* scale=4 */
      igauss_y = rand("IGAUSS", 2, 3);       /* lambda=2, mu=3 */
      weib_y = rand("Weibull", 2, 5);        /* shape=2, scale=5 */
      output;
   end;
   drop i;
run;

/* Load MLE package */
*%include "&MLE_path/MLE_Define.sas";

proc iml;
load module=_ALL_;

/* Read all test data */
use testdata;
read all var {beta_y gamma_y logn_y expo_y igauss_y weib_y};
close testdata;

/*=================================================================*/
/* TEST 1: BETA (alpha=2, beta=5)                                  */
/*=================================================================*/
print "=============== TEST 1: BETA (alpha=2, beta=5), nObs=&nObs ===============";

Fit = MLE_Fit("Beta", beta_y);
call MLE_Summary(Fit, 0, 1, 0.05);

submit;
   proc nlmixed data=testdata DF=%eval(&nObs-2);
      parms alpha =2 beta= 5;
      bounds alpha > 0, beta > 0;
      ll = logpdf("Beta", beta_y, alpha, beta);
      model beta_y ~ general(ll);
      ods select ParameterEstimates;
   run;
endsubmit;

/*=================================================================*/
/* TEST 2: GAMMA (alpha=3, lambda=2)                               */
/*=================================================================*/
print "=============== TEST 2: GAMMA (alpha=3, lambda=2), nObs=&nObs ===============";
Fit = MLE_Fit("Gamma", gamma_y);
call MLE_Summary(Fit, 0, 1, 0.05);

submit;
   proc nlmixed data=testdata DF=%eval(&nObs-2);
      parms alpha=3 lambda=2;
      bounds alpha > 0, lambda > 0;
      ll = logpdf("gamma", gamma_y, alpha, lambda);
      model gamma_y ~ general(ll);
      ods select ParameterEstimates;
   run;
endsubmit;

/*=================================================================*/
/* TEST 3: LOGNORMAL (mu=1, sigma=0.5)                             */
/*=================================================================*/
print "=============== TEST 3: LOGNORMAL (mu=1, sigma=0.5), nObs=&nObs ===============";
Fit = MLE_Fit("Lognormal", logn_y);
call MLE_Summary(Fit, 0, 1, 0.05);

submit;
   proc nlmixed data=testdata DF=%eval(&nObs-2);
      parms mu=1 sigma=0.5;
      bounds sigma > 0;
      ll = logpdf("lognormal", logn_y, mu, sigma);
      model logn_y ~ general(ll);
      ods select ParameterEstimates;
   run;
endsubmit;

/*=================================================================*/
/* TEST 4: EXPONENTIAL (sigma=4)                                   */
/*=================================================================*/
print "=============== TEST 4: EXPONENTIAL (sigma=4), nObs=&nObs ===============";
Fit = MLE_Fit("Expo", expo_y);
call MLE_Summary(Fit, 0, 1, 0.05);

submit;
   proc nlmixed data=testdata DF=%eval(&nObs-1);
      parms sigma=4;
      bounds sigma > 0;
      ll = logpdf("expo", expo_y, sigma);
      model expo_y ~ general(ll);
      ods select ParameterEstimates;
   run;
endsubmit;

/*=================================================================*/
/* TEST 5: IGAUSS (lambda=2, mu=3)                                 */
/*=================================================================*/
print "=============== TEST 5: IGAUSS (lambda=2, mu=3), nObs=&nObs ===============";
Fit = MLE_Fit("IGauss", igauss_y);
call MLE_Summary(Fit, 0, 1, 0.05);

submit;
   proc nlmixed data=testdata DF=%eval(&nObs-2) FD=CENTRAL;
      parms lambda=2 mu=3;
      bounds lambda > 0, mu > 0;
      ll = logpdf("IGauss", igauss_y, lambda, mu);
      model igauss_y ~ general(ll);
      ods select ParameterEstimates;
   run;
endsubmit;

/*=================================================================*/
/* TEST 6: WEIBULL (c=2, lambda=5)                                 */
/*=================================================================*/
print "=============== TEST 6: WEIBULL (c=2, lambda=5), nObs=&nObs ===============";
Fit = MLE_Fit("Weibull", weib_y);
call MLE_Summary(Fit, 0, 1, 0.05);

submit;
   proc nlmixed data=testdata DF=%eval(&nObs-2);
      parms c=2 lambda=5;
      bounds c > 0, lambda > 0;
      ll = logpdf("Weibull", weib_y, c, lambda);
      model weib_y ~ general(ll);
      ods select ParameterEstimates;
   run;
endsubmit;

ods select all;
title;
quit;
