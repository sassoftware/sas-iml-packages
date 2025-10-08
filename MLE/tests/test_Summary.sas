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
