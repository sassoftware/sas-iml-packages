/*****************************************************/
/* Utility functions for PROBMVN package             */
/*****************************************************/

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

/* Sigma is SPD, no missing. mu has no missing values. */
start mvn_IsValidParms(Sigma, mu);
   if ^mvn_IsSym(Sigma) then do;
      run ErrorToLog( "The Sigma parameter must be symmetric.");
      return( 0 );
   end;
   if ncol(mu) ^= ncol(Sigma) then do;
      run ErrorToLog( "The mu and Sigma parameters are not compatable dimensions.");
      return( 0 );
   end;
   if any(Sigma=.) then do;
      run ErrorToLog("The Sigma parameter cannot contain missing values.");
      return( 0 );
   end;
   if any(mu=.) then do;
      run ErrorToLog("The mu parameter cannot contain missing values.");
      return( 0 );
   end;
   if ^mvn_IsSPD(Sigma) then do;
      run ErrorToLog( "The Sigma parameter must be positive definite.");
      return( 0 );
   end;
   return( 1 );
finish;

/* opt[1] is QMC absolute error tolerance; opt[2] is binary flag for error estimate. */
start mvn_IsValidParmsOpt(opt);
   if type(opt) ^= 'N' then do;
      run ErrorToLog( "The opt parameter must be numeric." );
      return( 0 );
   end;
   if opt[1] ^= . then do;
      if opt[1] < 1E-5 | opt[1] > 1E-2 then do;
         run ErrorToLog( "The first element of opt must be missing or in the range [1E-5, 1E-2]." );
         return( 0 );
      end;
   end;
   return( 1 );
finish;

/* dim(U)=dim(Sigma), U is not missing, size for CDFMVN.*/
start mvn_IsValidParmsCDF(U, Sigma, mu);
   isValid = mvn_IsValidParms(Sigma, mu);
   if ^isValid then return( 0 );
   if ncol(U) ^= ncol(Sigma) then do;
      run ErrorToLog( "The U and Sigma parameters are not compatable dimensions.");
      return( 0 );
   end;
   if any(U=.) then do;
      run ErrorToLog("The U parameter cannot contain missing values.");
      return( 0 );
   end;
   if ncol(U)<2 | ncol(U) > 20 then do;
      run ErrorToLog( "CDFMVN supports problems between 2 and 20 dimensions.");
      return( 0 );
   end;
   return( 1 );
finish;

/* dim(L)=dim(U)=dim(Sigma), size for PROBMVN */
start mvn_IsValidParmsProbmvn(L, U, Sigma, mu);
   isValid = mvn_IsValidParms(Sigma, mu);
   if ^isValid then return( 0 );
   if nrow(L) > 1 | nrow(U) > 1 then do;
      run ErrorToLog( "The PROBMVN package does not support multiple limits of integration." );
      return( 0 );
   end;
   if ncol(L) ^= ncol(Sigma) then do;
      run ErrorToLog( "The L and Sigma parameters are not compatable dimensions.");
      return( 0 );
   end;
   if ncol(L) ^= ncol(U) then do;
      run ErrorToLog( "The L and U parameters are not compatable dimensions.");
      return( 0 );
   end;
   if ^IsValidRectLimits(L, U) then do;
      run ErrorToLog( "The L and U parameters must satisfy L[i] < U[i].");
      return( 0 );
   end;
   if ncol(U)<2 | ncol(U) > 100 then do;
      run ErrorToLog( "PROBMVN supports problems between 2 and 100 dimensions.");
      return( 0 );
   end;
   return( 1 );
finish;

/* IsValidRectLimits: validate the lower (L) and upper (U) limit vectors for PROBMVN_MOD.
   Missing values are explicitly allowed: a missing L[i] represents -Infinity
   and a missing U[i] represents +Infinity.
   Returns 1 if valid, 0 otherwise.
   Checks:
   1. L and U have the same number of elements.
   2. For every dimension i where both L[i] and U[i] are non-missing, L[i] < U[i].
*/
start IsValidRectLimits(L, U);
   if ncol(L) ^= ncol(U) then do;
      run ErrorToLog("The L and U parameters must have the same number of columns.");
      return( 0 );
   end;
   if nrow(L) ^= nrow(U) then do;
      run ErrorToLog("The L and U parameters must have the same number of rows.");
      return( 0 );
   end;
   n = ncol(L);
   do j = 1 to nrow(L);
       do i = 1 to n;
          if L[j,i] ^= . & U[j,i] ^= . then do;
             if L[j,i] >= U[j,i] then do;
                run ErrorToLog("L[j,i] must be < U[j,i] for every non-missing pair of limits.");
                return( 0 );
             end;
          end;
       end;
   end;
   return( 1 );
finish;

/*==============================================================*/
/*    Clip an integration limit into [-delta, delta], where     */
/*    delta ~ 8.125 is chosen so that                           */
/*    SDF("Normal", delta) ~ constant("maceps")                 */
/*==============================================================*/
start ClipLimit(x, delta=8.125);
   return( -delta <> (x >< delta) );
finish;

/* Sigma is a kxk covariance matrix and b is a row vector with k elements.
   Convert b to inv(D)*(b-mu), where D is the diagonal matrix of
   standard deviations: sqrt(vecdiag(Sigma))
   Note: b can have multiple rows. b can also have missing values.
*/
start Xform_Limits_Cov2Corr( b, Sigma, mu=j(1,ncol(Sigma),0) );
   D = sqrt(vecdiag(Sigma));
   c = (b - mu) /rowvec(D);
   return c;
finish;

/* Monte Carlo simulation of P( L < X < U | X~MVN(mu,Sigma) )
   where L and U are row vectors and a missing value 
   represents -Infinity in L and represents +Infinity in U.
   NOTE: If N is very large, a straightforward implementation uses a lot of memory.
   Therefore, perform the computations in blocks of size MaxN to
   prevent out-of-memory errors. It actually speeds up the computation, too!
*/
start MC_PROBMVN(N, Lcov, Ucov, Sigma, mu=j(1,ncol(Sigma),0));
   MaxN = 2E5;
   /* standardize the parameters to the correlation scale */
   L = Xform_Limits_Cov2Corr(Lcov, Sigma, mu);
   U = Xform_Limits_Cov2Corr(Ucov, Sigma, mu);
   R = cov2corr(Sigma);
   N_remain = N;
   nIter = ceil(N / MaxN);
   sum = 0;
   do j = 1 to nIter;
      K = min(MaxN, N_remain);
      inRegion = j(K,1,1);
      X = randnormal(K, j(1,ncol(R),0), R);
      /* Note: The '<' operator works correctly with a missing value
         on the left. The expression (a<y & y<b) is correct if a=.;
         However, if b=., you need to use (a<y). */
      do i = 1 to ncol(L);
         if U[i]=. then 
            v = (L[i] < X[,i]);
         else
            v = ((L[i] < X[,i]) & (X[,i] < U[i]));
         inRegion = inRegion & v;
      end;
      N_remain = N_remain - K;
      sum = sum + sum(inRegion);
   end;
   return sum / N;
finish;

/* return a list with the MC est and a 95% CL. The list looks like this:
   [prob, lower95, upper95] */
start MC_PROBMVN_CL(N, L, U, Sigma, mu=j(1,ncol(Sigma),0));
   prob_MC = MC_PROBMVN(N, L, U, Sigma, mu);
   SE_MC = sqrt( prob_MC * (1-prob_MC)/N );
   Lower95 = prob_MC - 1.96*SE_MC;
   Upper95 = prob_MC + 1.96*SE_MC;
   return( [prob_MC, Lower95, Upper95] );
finish;

/* Helper module used by tests. Prints a message is a tet fails.
   A test that passes is silent. */
start check_test(test_name, prob, correct, tol=1E-3);
   maxDiff = max(abs(prob-correct));
   if maxDiff > tol then do;
      msg = cat("--- ",test_name, " --- FAILS ---");
      print msg[L=""], maxDiff prob correct;
   end;
finish;

store module=(
   mvn_IsSym 
   mvn_IsSPD 
   mvn_IsCorr 
   mvn_IsValidParms
   mvn_IsValidParmsOpt
   mvn_IsValidParmsCDF 
   mvn_IsValidParmsProbmvn 
   ClipLimit
   Xform_Limits_Cov2Corr
   IsValidRectLimits
   MC_PROBMVN
   MC_PROBMVN_CL
   check_test
);
QUIT;
