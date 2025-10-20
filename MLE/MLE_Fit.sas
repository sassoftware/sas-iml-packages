/*********************************************************/
/* MLE_MLE.sas                                           */
/* Define functions that return the maximum likelihood   */
/* estimates for parameters in the built-in              */
/* distributions                                         */
/*********************************************************/


/***********************************/
/* MLE and MLE_FIT have the same syntax. The difference is that 
   MLE returns ONLY the MLE estimate, whereas 
   MLE_FIT returns a FItObj. The FitObj contains many statistics, 
   such as LL, grad(LL), Hess(LL),
   and standard errors evaluated at a specified tuple of parameter values.
*/
/***********************************/


/* Helper function to validate and prepare MLE input parameters
   Returns: a list containing validated initial_point, BoundsMatrix, and OptimMethodSetting
   
   This function consolidates the parameter validation logic used by both MLE and MLE_Fit.
   Note: param0_provided, Bounds_provided, OptimMethod_provided are flags (1/0) indicating 
   whether the optional parameters were provided by the caller.
*/
start mle_ValidateInputs(DistName, y, param0, param0_provided, Bounds, Bounds_provided, OptimMethod, OptimMethod_provided);
   /* Get distribution keyword to determine if built-in or user-defined */
   keyword = lik_dist_keyword(DistName);
   
   /* param0: initial guess, MoM for built-in, or user-defined MoM function name */
   if ^param0_provided then do;
      /* No param0 provided - use built-in MoM if available */
      if keyword ^= ' ' then do;  /* Built-in distribution - use MoM */
         initial_point = MLE_MoM(DistName, y);
      end;
      else do;  /* User-defined distribution - require param0 */
         STOP "ERROR: param0 is required for user-defined distribution " + DistName + ". Provide numeric initial values or MoM function name.";
      end;
   end;
   else do;
      /* Check if param0 is a string (function name) or numeric (initial values) */
      if type(param0) = 'C' then do;  /* Character - treat as MoM function name */
         %EVALFUNC1(initial_point, param0, y);
         initial_point = colvec(initial_point);
      end;
      else do;  /* Numeric - use as initial values */
         initial_point = colvec(param0);
         /* Validate param0: must be numeric, non-missing */
         if any(initial_point = .) then 
            STOP "ERROR: param0 contains missing values.";
         
         /* Validate parameter count for built-in distributions only */
         if keyword ^= ' ' then do;  /* Built-in distribution */
            parmNames = lik_dist_parmnames(DistName);
            if nrow(initial_point) ^= ncol(parmNames) then 
               STOP "ERROR: param0 has incorrect number of parameters. Expected " + char(ncol(parmNames),2) + " for " + DistName + " distribution.";
         end;
         /* For user-defined distributions, we cannot validate parameter count */
      end;
   end;

   /* Bounds: user-supplied bounds or default is j(2, nrow(initial_point), .) */
   if ^Bounds_provided then do;
      BoundsMatrix = j(2, nrow(initial_point), .);
   end;
   else do;
      /* Validate Bounds dimensions before transposing */
      if nrow(Bounds) ^= nrow(initial_point) | ncol(Bounds) ^= 2 then do;
         STOP "ERROR: Bounds must be a " + char(nrow(initial_point),2) + " x 2 matrix.";
      end;
      /* user-supplied bounds are stored as a two-column matrix, while NLP routines expect a two-row matrix */
      BoundsMatrix = T(Bounds);
      /* Optionally: check lower <= upper for each parameter */
      if any(BoundsMatrix[1,] > BoundsMatrix[2,]) then do;
         STOP "ERROR: Lower bounds exceed upper bounds in Bounds matrix.";
      end;
   end;

   /* OptimMethod: user-supplied optimization method or default is "NLPQN" */
   if ^OptimMethod_provided then OptimMethodSetting = "NLPQN";
   else OptimMethodSetting = upcase(strip(OptimMethod));

   /* Return validated parameters as a list */
   names = {"initial_point", "BoundsMatrix", "OptimMethodSetting"};
   result = ListCreate(names);
   result$"initial_point" = initial_point;
   result$"BoundsMatrix" = BoundsMatrix;
   result$"OptimMethodSetting" = OptimMethodSetting;
   
   return(result);
finish mle_ValidateInputs;


/* USAGE:
   CALL mle_OPTIM(rc, soln, method, ll_func, initial_point, bounds_matrix);
   
   WHERE:
   rc = return code from optimization
   soln = solution vector (parameter estimates)
   method = optimization method name ("NLPQN", "NLPNRA", "ACTIVE", "IP", "IPDIRECT")
   ll_func = log-likelihood function name
   initial_point = initial parameter guess
   bounds_matrix = 2-column matrix with lower and upper bounds

   Use %IF macro logic to call optimization methods based on SAS version and method name.
   This function handles the differences between SAS 9.4 and SAS Viya optimization routines.
*/
start mle_OPTIM(rc, soln, method, ll_func, initial_point, bounds_matrix);
   rc = -1;
   soln = j(nrow(initial_point), 1, .);
   errmsg = "ERROR: Unknown optimization method: " + kstrip(method) + ". ";
   %if &sysver. = 9.4 %then %do;
      errmsg = errmsg + "Valid methods in SAS 9.4 are NLPQN and NLPNRA.";
      validMethods = {"NLPQN", "NLPNRA"};
   %end;
   %else %do;
      errmsg = errmsg + "Valid methods in SAS Viya are NLPQN, NLPQN, NLPNRA, ACTIVE, IP, and IPDIRECT.";
      validMethods = {"NLPQN", "NLPNRA", "ACTIVE", "IP", "IPDIRECT"};
   %end;
   if type(method)^='C' then do;
      call PrintToLog(errmsg, 2);
      return;
   end;
   if ^element(upcase(method), validMethods) then do;
      call PrintToLog(errmsg, 2);
      return;
   end;
   
   if upcase(method) = "NLPQN" then 
      call NLPQN(rc, soln, ll_func, initial_point) OPT=1 BLC=bounds_matrix;
   else if upcase(method) = "NLPNRA" then 
      call NLPNRA(rc, soln, ll_func, initial_point) OPT=1 BLC=bounds_matrix;
   %if &sysver. ne 9.4 %then %do;
   else do;
      LowerBound = bounds_matrix[1,];
      UpperBound = bounds_matrix[2,];
      rowvec_initial_point = rowvec(initial_point);
      if      upcase(method) = "ACTIVE" then   opt={-1, 0, 0};
      else if upcase(method) = "IP"     then   opt={-1, 0, 1};
      else if upcase(method) = "IPDIRECT" then opt={-1, 0, 2};
      call NLPSOLVE(rc, soln, ll_func, rowvec_initial_point) OPT=opt L=LowerBound U=UpperBound;
   end;
   %end;
finish mle_OPTIM;

/* Compute the MLE estimates that fit DistName model to the Y data:
   SYNTAX:
   Est = MLE(DistName, y, param0=, OptimMethod=, Bounds=);
   WHERE THE INPUT PARAMETERS ARE
   DistName = a string for a built-in function or the name of a user-defined LL function 
   Y = a column vector of numerical values
   Param0 = a (1 x n_p) guess, or a string. 
            Use "MoM" for built-in methods, otherwise it is the name of user-defined func.
            For most common distributions, n_p = 2.
            You can call 
            parm_names = lik_dist_parmnames(distName)
            to get the names and use ncol(parm_names) to get n_p.
            If param0 is missing, then "MoM" is used for built-in distributions.
   OptimMethod = "NLPQN" (or " ") is default, "NLPNRA", [for Viya: "IP", "Active", "IPDirect"]
   Bounds = a two-column matrix. For CUSTOM distribs only. First column is Lower bounds 
       (default is j(n_p,1,.)) and second column is upper bounds on params (default j(n_p,1,.))

   The function returns the MLE estimate for the distribution, given the data.
*/

start MLE(DistName, y, param0=, OptimMethod=, Bounds=);
  /* Initialize the data */
   run MLE_Init(y);

   /* Check which optional parameters were provided */
   param0_provided = ^isskipped(param0);
   Bounds_provided = ^isskipped(Bounds);
   OptimMethod_provided = ^isskipped(OptimMethod);
   
   /* Set default values for skipped parameters (will be ignored if flags are 0) */
   if ^param0_provided then _param0 = .; else _param0 = param0;
   if ^Bounds_provided then _Bounds = .; else _Bounds = Bounds;
   if ^OptimMethod_provided then _OptimMethod = " "; else _OptimMethod = OptimMethod;

   /* Validate and prepare input parameters using helper function */
   validated = MLE_ValidateInputs(DistName, y, _param0, param0_provided, _Bounds, Bounds_provided, _OptimMethod, OptimMethod_provided);
   initial_point = validated$"initial_point";
   BoundsMatrix = validated$"BoundsMatrix";
   OptimMethodSetting = validated$"OptimMethodSetting";

   /* Get the log-likelihood function name */
   ll_func = lik_func_name(DistName, "LL");

   /* Use the CALLOPTIM macro to handle optimization method calls based on SAS version */
   CALL mle_OPTIM(rc, soln, OptimMethodSetting, ll_func, initial_point, BoundsMatrix);

   run MLE_End();
   
   return( soln );
finish MLE;


/* Main function that returns a fit object */
start MLE_Fit(DistName, y,  param0=, OptimMethod=, Bounds=)  global(gMLE_y);

   /* Check which optional parameters were provided */
   param0_provided = ^isskipped(param0);
   Bounds_provided = ^isskipped(Bounds);
   OptimMethod_provided = ^isskipped(OptimMethod);
   
   /* Set default values for skipped parameters (will be ignored if flags are 0) */
   if ^param0_provided then _param0 = .; else _param0 = param0;
   if ^Bounds_provided then _Bounds = .; else _Bounds = Bounds;
   if ^OptimMethod_provided then _OptimMethod = " "; else _OptimMethod = OptimMethod;

   /* Validate and prepare input parameters using helper function */
   validated = MLE_ValidateInputs(DistName, y, _param0, param0_provided, _Bounds, Bounds_provided, _OptimMethod, OptimMethod_provided);
   initial_point = validated$"initial_point";
   BoundsMatrix = validated$"BoundsMatrix";
   OptimMethodSetting = validated$"OptimMethodSetting";

   /* Call MLE with validated parameters - note: we pass the nx2 Bounds format, not the 2xn BoundsMatrix */
   est = MLE(DistName, y, initial_point, OptimMethodSetting, T(BoundsMatrix)); 

   names = {"Dist", "y", "ParmNames", "Estimate", "LL", 
            "Grad", "Hessian", "StdErr", "Crit", "CritNames"};
   FitObj = ListCreate(names);  /* named list */
   FitObj$"Dist" = lik_dist_name(DistName);
   FitObj$"y"    = gMLE_y;
   FitObj$"ParmNames" = lik_dist_parmnames(DistName);
   FitObj$"Estimate"  = est;
   run MLE_LL_and_Deriv(LL, Grad, Hess, DistName, est);
   FitObj$"LL"      = LL;
   FitObj$"Grad"    = Grad;
   FitObj$"Hessian" = Hess;
  /* StdErr: standard errors for estimates */
   cov = inv(-Hess);
   StdErr = rowvec(sqrt(vecdiag(cov)));
   FitObj$"StdErr" = StdErr;
   
   /* Crit, CritNames: LL-based fit criteria */
   n = nrow(y);
   k = nrow(est);
   N2LL = -2*LL;
   AIC = 2*k - 2*LL;
   BIC = k*log(n) - 2*LL;
   AICC = AIC + (2*k*(k+1))/(n-k-1);
   FitObj$"Crit" = N2LL // AIC // AICC // BIC;
   FitObj$"CritNames" = {'-2*LL', 'AIC', 'AICC', 'BIC'};
   return( FitObj );
finish MLE_Fit;

store module=(
      mle_ValidateInputs
      mle_OPTIM
      MLE
      MLE_Fit
);