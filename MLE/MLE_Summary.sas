/***********************************/
/* Define printing functions */
/***********************************/

/* MLE_Summary: Summary of MLE fit results
   
   SYNTAX:
   call MLE_Summary(FitObj, printOpt=, showCI=, alpha=, title=, showZ=, showP=);
   
   PARAMETERS:
   FitObj   = List object returned by MLE_Fit()
   printOpt = Print level (0-2):
              0 = Minimal: Only parameter estimates
              1 = Basic: Estimates + fit criteria (default)
              2 = Extended: + Log-likelihood, gradient, Hessian
   showCI   = 1 to show confidence intervals, 0 otherwise (default=0)
   alpha    = Significance level for CI (default=0.05)
   title    = Custom title for output (default uses distribution name)
   showZ    = 1 to show Wald Z statistics (default=1)
   showP    = 1 to show p-values (default=1)
*/

start MLE_Summary(L, printOpt=1, showCI=0, alpha=0.05, title="", showZ=1, showP=1);
   /* Validate input */
   namesL = ListGetName(L);
   if ncol(namesL)=0 then do;
      call PrintToLog("MLE_Summary: Input is not a valid list.", 2);
      return;
   end;
   
   /* Extract values from list */
   dist = L$"Dist";
   parmNames = L$"ParmNames";
   estimate = L$"Estimate";
   stdErr = L$"StdErr";
   crit = L$"Crit";
   critNames = L$"CritNames";
   
   /* Ensure column-vector shapes for printing */
   estimate = colvec(estimate);
   stdErr = colvec(stdErr);

   /* Validate and setup parameter names */
   k = nrow(estimate);
   if type(parmNames)^='C' | ncol(parmNames)^=k then do;
      /* Fallback: create generic parameter names */
      parmNames = j(1, k, "");
      do i = 1 to k;
         parmNames[i] = "Parm" + strip(char(i, 3));
      end;
   end;

   /* Set title for output */
   if strip(title)^="" then 
      titleText = title;
   else 
      titleText = "MLE Fit Summary: " + dist + " Distribution";
   
   /* Use SAS title statement for GUI display */
   stmt = 'title "' + titleText + '";';
   call execute(stmt);

   /* Build Parameter Estimates table */
   PE = estimate || stdErr;
   cNames = {"Estimate" "StdErr"};

   /* Add Z statistic and p-value (Wald test) */
   if showZ then do;
      z = j(k, 1, .);
      idx = loc(stdErr>0);
      if ncol(idx)>0 then z[idx] = estimate[idx] / stdErr[idx];
      PE = PE || z;
      cNames = cNames || {"Z"};
      
      if showP then do;
         p = j(k, 1, .);
         if ncol(idx)>0 then p[idx] = 2#(1 - cdf("Normal", abs(z[idx]), 0, 1));
         PE = PE || p;
         cNames = cNames || {"Pr>|Z|"};
      end;
   end;

   /* Add confidence intervals to table */
   if showCI then do;
      zCrit = probit(1 - alpha/2);
      lower = estimate - zCrit # stdErr;
      upper = estimate + zCrit # stdErr;
      PE = PE || lower || upper;
      clPct = char(100*(1-alpha), 6, 1);
      lowName = "Lower " + strip(clPct) + "% CL";
      upName  = "Upper " + strip(clPct) + "% CL";
      cNames = cNames || (lowName || upName);
   end;

   /* Print Parameter Estimates table */
   print PE [label="Parameter Estimates" rowname=parmNames colname=cNames format=8.4];

   /* printOpt >= 1: Fit Criteria */
   if printOpt >= 1 then
      print crit [label="Fit Criteria" colname=critNames format=8.4];

   /* printOpt >= 2: Optimization Details */
   if printOpt >= 2 then do;
      LL = L$"LL";
      grad = L$"Grad";
      hessian = L$"Hessian";
      
      print LL [label="Log-Likelihood at MLE" format=12.6];
      print grad [label="Gradient at MLE (should be near zero)" colname=parmNames format=best8.];
      
      grad_norm = sqrt(ssq(grad));
      print grad_norm [label="Gradient Norm" format=best12.];
      
      print hessian [label="Hessian Matrix at MLE" colname=parmNames rowname=parmNames format=10.4];
      
      call eigen(eigval, eigvec, hessian);
      print eigval [label="Hessian Eigenvalues (should be negative)" format=10.4];
      
      if all(eigval < 0) then
         print "NOTE: Hessian is negative definite (local maximum verified)";
      else
         print "WARNING: Hessian is not negative definite - may not be at local maximum";
   end;
   
   /* Reset title */
   call execute('title;');
   
finish;

store module=(
      MLE_Summary
);
