proc iml;
%include "MLE_define.sas";
%include "MLE_Fit.sas";
QUIT;

/* test calls to the top-level subroutine */
PROC IML;
load module=(MLE);

QUIT;
