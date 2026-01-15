/* Example of using the MLE library */

/* Step 0. Download library from GitHub. 
   An example is given in the file test_Install.sas in this directory.
   For details, see
   https://blogs.sas.com/content/iml/2023/02/06/git-share-sas-programs.html

   -OR-

   Define the MLE_Path variable to point to the directory that contains 
   the MLE source code and manually define the path as follows:

%let MLE_Path = u:\gitpp\DEV\sas-iml-packages\MLE;    * a Windows-style path uses '\' separators;
proc iml;
   %include "&MLE_path\MLE_Define.sas"; 
quit;
*/


proc iml;
load module=_all_;     /* load the MLE library */
use sashelp.heart;
   read all var "Systolic";
close;

/* Primary use case: call top-level MLE routine to get estimates */
/* parameter estimates */
gamma_est = MLE("Gamma", Systolic);  /* default guess is MoM */
parmNames = lik_dist_parmnames("Gamma");
print gamma_est[r=parmNames];

/* you can get the MoM directly and use it (or another guess)
   est_MoM = MLE_MoM("Dist", y);
*/
gamma_Mom = MLE_MoM("Gamma", Systolic);  
gamma_est2 = MLE("Gamma", Systolic, gamma_Mom);  /* specify a guess */
print gamma_MoM[r=parmNames], gamma_est2[r=parmNames];

/* The MLE_Fit function returns a list with many items */
L_gamma = MLE_Fit("Gamma", Systolic);
run MLE_Summary(L_gamma);     /* print basic results */
run MLE_Summary(L_gamma, 2);  /* print details about the MLE solution */

title "Overlay the Histogram and Fitted Model";
run MLE_Plot(L_gamma);        /* overlay the curve on a histogram */

/* You can fit other distributional models and overlay their estimates */
L_LogN = MLE_Fit("Lognormal", Systolic);
L_Weibull = MLE_Fit("Weibull", Systolic);
title "Overlay the Three Fitted Models";
run MLE_Plot(L_gamma, L_LogN, L_Weibull);        /* overlay multiple curves */

/* you can retreive and combine information from summary objects */
LL_Table = L_gamma$"Crit" || L_LogN$"Crit" || L_Weibull$"Crit";
print LL_Table[r=(L_gamma$"CritNames") c={"Gamma" "Lognormal" "Weibull"}];

QUIT;