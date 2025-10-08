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
