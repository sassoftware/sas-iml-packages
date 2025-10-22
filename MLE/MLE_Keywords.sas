/*********************************************************/
/* MLE_Keywords.sas                                      */
/* Define functions related to distribution keywords.    */
/*********************************************************/


start MLE_Init(y, distname=) global(gMLE_y);   
   gMLE_y = CompleteCases(colvec(y), "extract");
   n = nrow(gMLE_y);
   if n<3 then do;
      run PrintToLog("MLE_INIT: At least 3 nonmissing observations are required.",2);
      return 0;
   end;
   if ^IsSkipped(distname) then do;
      /* check if data are valid for the specified distribution */
      isValid = MLE_LL_ValidData(distname, gMLE_y);
      if isValid=0 then do;
         msg = "MLE_INIT: The data contain values that are invalid for the " + 
                strip(distname) + " distribution.";
         run PrintToLog(msg,2);
         return 0;
      end;
   end;
   return 1;
finish;

start MLE_End(_=0) global(gMLE_y);   
   free gMLE_y;
finish;

/***************************************************/

/* By default, return a unique KEYWORD for any valid distribution name. 
   KEYWORD VALUE:
   1  :  Return the 4-character KEYWORD used in the SAS procedures.
   2  :  Return the SUFFIX used by the functions in the MLE library. 
   3  :  Return the string used by the PDF and CDF functions in SAS
   Return ' ' if the name is not valid. */
start lik_dist_keyword(distname, keyword=1);
   tbl = {
   /* valid inputs | 4-char keywd | std suffix | PDF Name */
      'Beta'         'BETA'         'Beta'      'Beta',
      'Exponential'  'EXPO'         'Expo'      'Expo', 
      'Exp'          'EXPO'         'Expo'      'Expo',
      'Expo'         'EXPO'         'Expo'      'Expo',
      'Gamma'        'GAMM'         'Gamma'     'Gamma',
      'Gumbel'       'GUMB'         'Gumbel'    'Gumbel',
      'IGauss'       'IGAU'         'IGauss'    'IGauss',
      'InvGauss'     'IGAU'         'IGauss'    'IGauss',
      'Wald'         'IGAU'         'IGauss'    'IGauss',
      'Lognormal'    'LN2 '         'LN2'       'Lognormal',
      'LN'           'LN2 '         'LN2'       'Lognormal',
      'LN2'          'LN2 '         'LN2'       'Lognormal',
      'LN3'          'LN3 '         'LN3'       'Lognormal',
      'Normal'       'NORM'         'Normal'    'Normal',
      'Weibull'      'WEI2'         'Weib2'     'Weibull',
      'Weib2'        'WEI2'         'Weib2'     'Weibull',
      'Weib3'        'WEI3'         'Weib3'     'Weibull'
   };
   idx = loc( upcase(distname) = upcase(tbl[,1]) );
   if ncol(idx)=0 then   /* not a built-in distribution */
      return( ' ' );
   if keyword=2 | keyword=3 then 
      return( kstrip(tbl[idx, keyword+1]) );
   else 
      return( kstrip(tbl[idx, 2]) );
finish;

start lik_dist_suffix(distname);
   return lik_dist_keyword(distname, 2);
finish;
start lik_dist_name(distname);
   return lik_dist_keyword(distname, 3);
finish;

/* The built-in functions in the MLE library follow a naming convention:
   MLE_LL_Dist: the log-likelihood function for 'Dist'. For example, MLE_LL_Expo.
   MLE_Mom_Dist: the method-of-moments estimators for 'Dist'
   MLE_ValidData_Dist: function to check if data are valid for 'Dist'
   If you are implementing your own custom distribution, the naming convention is:
   LL_Dist, MoM_Dist, and ValidData_Dist.  For example, LL_MyDist, MoM_MyDist, etc.
*/
start lik_func_name(distname, role);
   if upcase(role)="LL" then do;
      prefix = "lik_LL_";
      user_prefix = "LL_";
   end;
   else if upcase(role)="MOM" then do;
      prefix = "lik_MOM_";
      user_prefix = "MOM_";
   end;
   else if upcase(role)="VALIDDATA" then do;
      prefix = "lik_ValidData_";
      user_prefix = "ValidData_";
   end;
   suffix = lik_dist_suffix(distname);
   if missing(suffix) then do;
      /* User-defined distribution - check if name already has the prefix */
      distname_upper = upcase(strip(distname));
      if substr(distname_upper, 1, length(user_prefix)) = user_prefix then
         func_name = kstrip(distname);  /* Already has prefix, use as-is */
      else
         func_name = kstrip(upcase(role)) + "_" + kstrip(distname);  /* Add prefix by convention */
   end;
   else
      func_name = prefix + kstrip(suffix);    /* built-in func */
   return( func_name );
finish;


/* Given a KEYWORD, return the parameter names of the distribution.
   You can get the number of parameters by NCOL(result).
*/
start lik_dist_parmnames(distName);
   tbl = {
   /* 4-char keyword | Parm1 |  Parm2 |  Parm3 */
      'BETA'           'alpha'  'beta'   ' ',
      'EXPO'           'sigma'  ' '      ' ',
      'GAMM'           'alpha'  'lambda' ' ',
      'GUMB'           'mu'     'sigma'  ' ',
      'IGAU'           'lambda' 'mu'     ' ',
      'LN2 '           'mu'     'sigma'  ' ',
      'LN3 '           'mu'     'sigma'  'theta',
      'NORM'           'mu'     'sigma'  ' ',
      'WEI2'           'c'      'sigma'  ' ',
      'WEI3'           'c'      'sigma'  'theta'
   };
   distKEY = lik_dist_keyword(distName);
   idx = loc( upcase(distKEY) = upcase(tbl[,1]) );
   if ncol(idx)=0 then   /* not a built-in distribution */
      return( ' ' );
   r = tbl[idx, ];
   return( r[ ,2:countn(r)] );
finish;

store module=(
   MLE_Init
   MLE_End
   lik_dist_keyword
   lik_dist_suffix
   lik_dist_name
   lik_func_name
   lik_dist_parmnames
);
