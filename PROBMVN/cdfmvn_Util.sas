%macro DefineCDFMVN;
%if %sysevalf(&sysver = 9.4) %then %do;
   start cdfmvn(b, Sigma, mu=repeat(0,1,ncol(Sigma)));
      if ncol(Sigma)=2 then 
         prob = cdfbvn_mod(b, Sigma, mu);
      else if ncol(Sigma)=3 then 
         prob = cdftvn_mod(b, Sigma, mu);
      else
         prob = cdfmvn_mod(b, Sigma, mu);
      return prob;
   finish;
   store module=(cdfmvn);   
%end;
%mend;

/* Check the SYSVER macro to see if SAS 9.4 is running.
   In SAS Viya, the macro is empty and does nothing.
   In SAS 9.4, the macro defines a function that emulates the PrintToLog call.
   The syntax is as follows:
   call PrintToLog("This is a log message.");
   call PrintToLog("This is a note.", 0);
   call PrintToLog("This is a warning.", 1);
   call PrintToLog("This is an error.", 2);
*/
%macro DefinePrintToLog;
%if %sysevalf(&sysver = 9.4) %then %do;
start PrintToLog(msg,errCode=-1);
   if      errCode=0 then prefix = "NOTE: ";
   else if errCode=1 then prefix = "WARNING: ";
   else if errCode=2 then prefix = "ERROR: ";
   else prefix = "";
   stmt = '%put ' + prefix + msg + ';';
   call execute(stmt);
finish;
store module=(PrintToLog);
%end;
start ErrorToLog(msg);
   run PrintToLog(msg, 2);
finish;
store module=(ErrorToLog);   
%mend;

/* this program runs in SAS 9.4 or in SAS Viya */
proc iml;
%DefineCDFMVN;
%DefinePrintToLog; 

/* validate the arguments for CDFMVN:
   b does not support missing values
   Sigma is SPD
*/
start mvn_IsSym(A);
   if nrow(A) ^= ncol(A) then return(0);    /* A is not square */
   c = max(abs(A));
   sqrteps = constant('SqrtMacEps');
   return( all( abs(A-A`) < c*sqrteps ) );
finish;

start mvn_IsSPD(M);
   if ^mvn_IsSym(M) then return( 0 );
   U = root(M, "NoError");
   if any(U=.) then return( 0 );
   return( 1 );
finish;

start mvn_IsCorr(M);
   if ^mvn_IsSPD(M) then return( 0 );
   if any(vecdiag(M) ^= 1) then return ( 0 );
   return( 1 );
finish;

start mvn_IsValidParmsCDF(b, Sigma, mu);
   if ^mvn_IsSym(Sigma) then do;
      run ErrorToLog( "The Sigma parameter must be symmetric.");
      return( 0 );
   end;
   if ncol(b) ^= ncol(Sigma) then do;
      run ErrorToLog( "The b and Sigma parameters are not compatable dimensions.");
      return( 0 );
   end;
   if ncol(mu) ^= ncol(Sigma) then do;
      run ErrorToLog( "The mu and Sigma parameters are not compatable dimensions.");
      return( 0 );
   end;
   if any(b=.) | any(Sigma=.) | any(mu=.) then do;
      run ErrorToLog( "The parameters cannot contain missing values.");
      return( 0 );
   end;
   if ^mvn_IsSPD(Sigma) then do;
      run ErrorToLog( "The Sigma parameter must be positive definite.");
      return( 0 );
   end;
   return( 1 );
finish;

start mvn_IsValidParmsMVN(b, Sigma, mu);
   isValid = mvn_IsValidParmsCDF(b, Sigma, mu);
   if ^isValid then return( 0 );
   if ncol(b)<2 | ncol(b) > 20 then do;
      run ErrorToLog( "cdfmvn supports problems between 2 and 20 dimensions.");
      return( 0 );
   end;
   return( 1 );
finish;

/* Truncate a value x to be within the bounds ab=[a,b] */
start Clip(x, ab);
   a = ab[1]; b = ab[2];
   return( choose(x=., ., a <> (x >< b)) );
finish;

/* Convert from centered Covariance scale to Correlation scale.
   Assume that we have already verified that Sigma is SPD. Call as
   run StdizeCovToCorr(U, R.
                       b, Sigma, mu);  
   The input arguments are b, Sigma, and mu.
   After the call, U is the standardized upper limits, and 
   R is the correlation matrix. 
*/
start mvn_StdizeCovToCorr(U, R, b, Sigma, mu);
   n = nrow(Sigma);
   D = sqrt(vecdiag(Sigma));
   R = cov2corr(Sigma);
   /* assume b is a row vector or a matrix of row vectors */
   U = (b - mu) / rowvec(D);
finish;

store module=(
   mvn_IsSym 
   mvn_IsSPD 
   mvn_IsCorr 
   mvn_IsValidParmsCDF 
   mvn_IsValidParmsMVN 
   Clip 
   mvn_StdizeCovToCorr
);
QUIT;
