proc iml;
%include "&MLE_path.\MLE_util.sas";
%include "&MLE_path.\MLE_keywords.sas";
QUIT;

/* tests */
proc iml;
load module= _ALL_;

/* MLE_Init */
y = {1, ., 2, ., ., 3, .};
run MLE_Init(y);
if ^all(gMLE_y = {1,2,3}) then print "gMLE_y is wrong";

run MLE_Init({. . . });  /* ERROR */
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

/* get keywords and suffixes */
results = J(nrow(inputs), 4, BlankStr(15));
mattrib results[c={'Dist' 'Keyword' 'Suffix' 'Name'}];
do i = 1 to nrow(inputs);
   results[i,1] = inputs[i];
   results[i,2] = lik_dist_keyword(inputs[i]);
   results[i,3] = lik_dist_suffix(inputs[i]);
   results[i,4] = lik_dist_name(inputs[i]);
end;
print results;


/* construct name of function */
results = J(nrow(inputs), 3, BlankStr(25));
mattrib results[c={'Dist' 'LL Func' 'MoM Func'}];
do i = 1 to nrow(inputs);
   results[i,1] = inputs[i];
   results[i,2] = lik_func_name(inputs[i], "LL");
   results[i,3] = lik_func_name(inputs[i], "MoM");
end;
print results;


/* Test getting ParmNames */
results = J(nrow(inputs), 5, BlankStr(25));
mattrib results[c={'Dist' 'Num Parm' 'Parm 1' 'Parm 2' 'Parm 3'}];
do i = 1 to nrow(inputs);
   results[i,1] = inputs[i];
   names = lik_dist_parmnames(inputs[i]);
   numParms = ncol(names);
   results[i,2] = char(numParms);
   results[i,3:2+numParms] = names;
end;
print results;


QUIT;
