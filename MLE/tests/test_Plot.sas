proc iml;
%include "&MLE_path./MLE_define.sas";
QUIT;

proc iml;
load module=_all_;     /* load the MLE library */
use sashelp.heart;
   read all var "Systolic";
close;

/* Primary use Case: call top-level MLE routine to get estimates */
/* parameter estimates */
gamma_est = MLE("Gamma", Systolic);  /* default guess is MoM */
print gamma_est;

title "Plot from low-level routine";
run mle_Plot_Overlay(Systolic, "Gamma", rowvec(gamma_est));

L_gamma = MLE_Fit("Gamma", Systolic);
run MLE_Summary(L_gamma);

title "Plot from high-level routine";
run MLE_Plot(L_gamma);/* overlay the curve on a histogram */


L_weib = MLE_Fit("Weibull", Systolic);
run MLE_Summary(L_weib);

title "Histogram with Multiple Density Estimates";
run MLE_Plot(L_gamma, L_weib); /* overlay curves on the histogram */


*QUIT;

