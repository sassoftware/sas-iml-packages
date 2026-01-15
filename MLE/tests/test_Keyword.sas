/* BEFORE running this example, store the modules in the MLE package, as shown in test_Install.sas */

/* tests for keywords */
proc iml;
load module= _ALL_;
print "--- MLE_Keywords: A successfult test prints only 'TEST DONE' ---";

/* MLE_Init */
y = {1, ., 2, ., ., 3, .};
isValid = MLE_Init(y);
if ^isValid | ^all(gMLE_y = {1,2,3}) then print "gMLE_y is wrong";

isValid = MLE_Init({. . . });  /* ERROR */
if isValid then print "ERROR: should be invalid data";
run MLE_End();
if ^(ncol(gMLE_y)=0) then print "gMLE_y is not empty";

/* lik_dist_keywords */
inputs = {
      'Beta'        
      'Exponential' 
      'Exp'         
      'Expo'        
      'Gamma'       
      'Gumbel'      
      'IGauss'      
      'InvGauss'    
      'Wald'        
      'Lognormal'   
      'LN'          
      'LN2'         
      'LN3'         
      'Normal'      
      'Weibull'     
      'Weib2'       
      'Weib3'       
      'Custom1'
      'T'
}`;

/* Helper function to check for equality between two character matrices */
start CharMatEqual(A, B);
   if type(A)^='C' | type(B)^='C' |
      (nrow(A)^=nrow(B)) | (ncol(A)^=ncol(B)) then 
      return(0);
   return all(upcase(A)=upcase(B));
finish;


/* get keywords and suffixes */
results = J(nrow(inputs), 4, BlankStr(15));
mattrib results[c={'Dist' 'Keyword' 'Suffix' 'Name'}];
do i = 1 to nrow(inputs);
   results[i,1] = inputs[i];
   results[i,2] = lik_dist_keyword(inputs[i]);
   results[i,3] = lik_dist_suffix(inputs[i]);
   results[i,4] = lik_dist_name(inputs[i]);
end;
CORRECT = {
'Beta'        'BETA' 'Beta' 'Beta', 
'Exponential' 'EXPO' 'Expo' 'Expo', 
'Exp'         'EXPO' 'Expo' 'Expo', 
'Expo'        'EXPO' 'Expo' 'Expo', 
'Gamma'       'GAMM' 'Gamma' 'Gamma', 
'Gumbel'      'GUMB' 'Gumbel' 'Gumbel', 
'IGauss'      'IGAU' 'IGauss' 'IGauss', 
'InvGauss'    'IGAU' 'IGauss' 'IGauss', 
'Wald'        'IGAU' 'IGauss' 'IGauss', 
'Lognormal'   'LN2'  'LN2' 'Lognormal', 
'LN'          'LN2'  'LN2' 'Lognormal', 
'LN2'         'LN2'  'LN2' 'Lognormal', 
'LN3'         'LN3'  'LN3' 'Lognormal', 
'Normal'      'NORM' 'Normal' 'Normal', 
'Weibull'     'WEI2' 'Weib2' 'Weibull', 
'Weib2'       'WEI2' 'Weib2' 'Weibull', 
'Weib3'       'WEI3' 'Weib3' 'Weibull', 
'Custom1'     ' '  ' '  ' ',
'T'           ' '  ' '  ' '
};
if ^CharMatEqual(results, CORRECT) then 
   print "ERROR: lik_dist functions not correct", results, CORRECT;
else;

/* construct name of function */
results = J(nrow(inputs), 3, BlankStr(25));
mattrib results[c={'Dist' 'LL Func' 'MoM Func'}];
do i = 1 to nrow(inputs);
   results[i,1] = inputs[i];
   results[i,2] = lik_func_name(inputs[i], "LL");
   results[i,3] = lik_func_name(inputs[i], "MoM");
end;
CORRECT = {
'Beta'        'lik_LL_Beta' 'lik_MOM_Beta',
'Exponential' 'lik_LL_Expo' 'lik_MOM_Expo',
'Exp'         'lik_LL_Expo' 'lik_MOM_Expo',
'Expo'        'lik_LL_Expo' 'lik_MOM_Expo',
'Gamma'       'lik_LL_Gamma' 'lik_MOM_Gamma',
'Gumbel'      'lik_LL_Gumbel' 'lik_MOM_Gumbel',
'IGauss'      'lik_LL_IGauss' 'lik_MOM_IGauss',
'InvGauss'    'lik_LL_IGauss' 'lik_MOM_IGauss',
'Wald'        'lik_LL_IGauss' 'lik_MOM_IGauss',
'Lognormal'   'lik_LL_LN2' 'lik_MOM_LN2',
'LN'          'lik_LL_LN2' 'lik_MOM_LN2',
'LN2'         'lik_LL_LN2' 'lik_MOM_LN2',
'LN3'         'lik_LL_LN3' 'lik_MOM_LN3',
'Normal'      'lik_LL_Normal' 'lik_MOM_Normal',
'Weibull'     'lik_LL_Weib2' 'lik_MOM_Weib2',
'Weib2'       'lik_LL_Weib2' 'lik_MOM_Weib2',
'Weib3'       'lik_LL_Weib3' 'lik_MOM_Weib3',
'Custom1'     'LL_Custom1' 'MOM_Custom1',
'T'           'LL_T' 'MOM_T'
};
if ^CharMatEqual(results, CORRECT) then 
   print "ERROR: lik_func_name function not correct", results, CORRECT;
else;

/* Test getting ParmNames */
results = J(nrow(inputs), 5, BlankStr(25));
mattrib results[c={'Dist' 'Num Parm' 'Parm 1' 'Parm 2' 'Parm 3'}];
do i = 1 to nrow(inputs);
   results[i,1] = inputs[i];
   names = lik_dist_parmnames(inputs[i]);
   numParms = ncol(names);
   results[i,2] = kstrip(char(numParms));
   results[i,3:2+numParms] = names;
end;
CORRECT = {
'Beta'        '2' 'alpha' 'beta'  ' ',   
'Exponential' '1' 'sigma'  ' '  ' ',
'Exp'         '1' 'sigma'  ' '  ' ',
'Expo'        '1' 'sigma'  ' '  ' ',    
'Gamma'       '2' 'alpha' 'lambda'  ' ',   
'Gumbel'      '2' 'mu'    'sigma'  ' ',   
'IGauss'      '2' 'lambda' 'mu'  ' ',   
'InvGauss'    '2' 'lambda' 'mu'  ' ',   
'Wald'        '2' 'lambda' 'mu'  ' ',   
'Lognormal'   '2' 'mu' 'sigma'  ' ',   
'LN'          '2' 'mu' 'sigma'  ' ',   
'LN2'         '2' 'mu' 'sigma'  ' ',   
'LN3'         '3' 'mu' 'sigma' 'theta',
'Normal'      '2' 'mu' 'sigma'  ' ',   
'Weibull'     '2' 'c'  'sigma'  ' ',   
'Weib2'       '2' 'c'  'sigma'  ' ',   
'Weib3'       '3' 'c'  'sigma' 'theta',
'Custom1'     '1' ' '  ' '  ' ',
'T'           '1' ' '  ' '  ' '
};
if ^CharMatEqual(results, CORRECT) then 
   print "ERROR: lik_dist_parmnames function not correct", results, CORRECT;
else;

print "--- TEST DONE ---";
*QUIT;

