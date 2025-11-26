/* Examples in the documentation for the MLE package */

/* NOTE: You must define the MLE_Path variable to point to the directory that contains 
   the MLE source code BEFORE you %include the MLE_Define.sas file */
%let MLE_Path = u:\gitpp\DEV\sas-iml-packages\MLE;

/* Use %INCLUDE to read source code and STORE functions to current storage library */
proc iml;
%include "&MLE_path/MLE_Define.sas";  /* one file with all modules */
quit;

proc iml;
load module=_all_;     /* load the MLE library */

/* MLE example */
print "-----------------";
y = {3,4,8,3,7,5,3,3,6,6,1,2};
Normal_est = MLE("Normal", y);           /* default call */
Gamma_est = MLE("Gamma", y, {3.5 1.5});  /* initial guess */
LN_est = MLE("LN2", y, "MoM", "NLPNRA"); /* specify optimization */
print Normal_est[r={'mu' 'sigma'}],
      Gamma_est[r={'alpha' 'lambda'}],
      LN_est[r={'mu' 'sigma'}];

/* MLE_Fit and MLE_Summary examples */
print "-----------------";
y = {3,4,8,3,7,.,5,3,.,3,6,6,1,.,2};
L = MLE_Fit("Lognormal", y);
run MLE_Summary(L);

/* MLE_MoM example */
print "-----------------";
y = {3,4,8,3,7,.,5,3,.,3,6,6,1,.,2};
LogNormal_MoM = MLE_MoM("Lognormal", y);
print LogNormal_MoM[c={'mu' 'sigma'}];

/* MLE_LL example */
print "-----------------";
y = {3,4,8,3,7,.,5,3,.,3,6,6,1,.,2};
est = MLE("LN2", y);
LL = MLE_LL("LN2", est);
print LL;

/* MLE_Init and MLE_End examples */
print "-----------------";
y = {3,4,8,3,7,.,5,3,.,3,6,6,1,.,2};
isValid = MLE_Init(y, "Gamma");  /* remove missing values and check y > 0 */
LL = MLE_LL("Gamma", {4,1});  /* LL for these parameters */
print LL[L="Gamma LL for param={4,1}"];
run MLE_End();

y = {3,4,8,3,7,.,5,3,.,3,6,6,1,.,2,10,3,4,5,4,3};
L = MLE_Fit("Gamma", y);
Title "Example of Histogram and PDF for Gamma MLE";
run MLE_Plot(L);


*QUIT;

