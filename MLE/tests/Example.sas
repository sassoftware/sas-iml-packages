/* Top-down design of MLE library. How will customers use?
*/

/* Step 0. download library from GitHub. See
   https://blogs.sas.com/content/iml/2023/03/13/metalog-sas.html
*/
/* download  package from GitHub */
/* clone repo to WORK, or use permanent libref */
/* clone repository; if repository exists, skip download */
/*
options dlcreatedir;
%let repoPath = %sysfunc(getoption(WORK))/sas-iml-packages;  
%let MLE_path = &repoPath/MLE; 

data _null_;
if fileexist("&repoPath.") then 
   put 'Repository already exists; skipping the clone operation'; 
else do;
   put "Cloning repository 'sas-iml-packages'";
   rc = gitfn_clone("https://github.com/sassoftware/sas-iml-packages/", "&repoPath." ); 
end;
run;
*/

/* NOTE: You must define the MLE_Path variable to point to the directory that contains 
   the MLE source code BEFORE you %include the MLE_Define.sas file */
%let MLE_Path = u:\gitpp\DEV\sas-iml-packages\MLE;

/* Use %INCLUDE to read source code and STORE functions to current storage library */
proc iml;
%include "&MLE_path/MLE_Define.sas";  /* one file with all modules */
quit;


proc iml;
load module=_all_;     /* load the MLE library */
use sashelp.heart;
   read all var "Systolic";
close;

/* Primary use Case: call top-level MLE routine to get estimates */
/* parameter estimates */
gamma_est = MLE("Gamma", Systolic);  /* default guess is MoM */
parmNames = lik_dist_parmnames("Gamma");
print gamma_est[r=parmNames];

/* you can get the MoM directly and use it (or another guess)
   est_MoM = MLE_MoM("Dist", y);
*/
gamma_Mom = MLE_MoM("Gamma", Systolic);  
print gamma_MoM[r=parmNames];
gamma_est = MLE("Gamma", Systolic, gamma_Mom);  /* specify a guess */
print gamma_est[r=parmNames];

/* The MLE_Fit function returns a list with many items */
L_gamma = MLE_Fit("Gamma", Systolic);
run MLE_Summary(L_gamma);     /* print basic results */
run MLE_Summary(L_gamma, 2);  /* print details about the MLE solution */

title "Overlay the Histogram and Fitted Model";
run MLE_Plot(L_gamma);        /* overlay the curve on a histogram */

/* For low-level routines, we can call low-level routines:
   est_MoM = lik_MoM_DistName(Y)
   and 
   ll = lik_LL_DistName(param) 
   For the the LL we must first execute
   isValid = MLE_Init(Y);
   which defines the global variable gMLE_y.
   For completeness, we also include
   run MLE_End(Y);
   which frees the global variable gMLE_y.
   */
isValid = MLE_Init(Systolic);              /* create global variable */
gamma_MOM = lik_MoM_Gamma(Systolic); /* copies Systolic to gMLE_y */
ll = lik_LL_Gamma(gamma_MOM);        /* accesses gMLE_y */
print gamma_MOM[r=parmNames], ll;
run MLE_End();                       /* frees global variable */

QUIT;