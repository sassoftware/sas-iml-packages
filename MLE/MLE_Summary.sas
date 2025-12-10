/***********************************/
/* Define printing functions */
/***********************************/

/* MLE_Summary: Summary of MLE fit results
   
   SYNTAX:
   call MLE_Summary(FitObj, printOpt=1, showCI=0, alpha=0.05);
   
   PARAMETERS:
   FitObj   = List object returned by MLE_Fit()
   printOpt = Print level (0-2):
              0 = Minimal: Only parameter estimates
              1 = Basic: Estimates + fit criteria (default)
              2 = Extended: + Log-likelihood, gradient, Hessian
   showCI   = 1 to show confidence intervals, 0 otherwise (default=0)
   alpha    = Significance level for CI (default=0.05)
*/

start MLE_Summary(L, printOpt=1, showCI=0, alpha=0.05) global(G_DEBUG);
   IF G_DEBUG THEN  run PrintLoc();
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
   
   /* Ensure column-vector shapes for printing */
   estimate = colvec(estimate);
   stdErr = colvec(stdErr);

   /* Validate and setup parameter names */
   k = nrow(estimate);
   if type(parmNames)^='C' | ncol(parmNames)^=k then
      /* Fallback: create generic parameter names */
      parmNames = "Parm1":("Parm" + strip(char(k, 3)));

   /* Build Parameter Estimates table */
   PE = estimate || stdErr;
   cNames = {"Estimate" "StdErr"};

   /* Add optional columns of statistics and p-value */
   if showCI then do;
      stat = j(k, 1, .);
      idx = loc(stdErr>0);
      if ncol(idx)>0 then stat[idx] = estimate[idx] / stdErr[idx];
      p = j(k, 1, .);

      use_t_statistic = 1;
      if use_t_statistic then do;
         /* Compute degrees of freedom for t-distribution: df = n - k */
         y = L$"y";
         n = nrow(y);
         DF = n - k;
         PE = PE || j(k, 1, DF) ;
         cNames = cNames || {"DF"};
         /* t statistic and p-value*/
         cNames = cNames || {"t Value"};
         if ncol(idx)>0 then p[idx] = 2#(1 - cdf("T", abs(stat[idx]), DF));
         PE = PE || p;
         cNames = cNames || {"Pr>|t|"};
         Crit = quantile("T", 1 - alpha/2, DF);
      end;
      else do;
         /* Wald Z statistic and p-value */
         cNames = cNames || {"Z"};
         if ncol(idx)>0 then p[idx] = 2#(1 - cdf("Normal", abs(z[idx]), 0, 1));
         PE = PE || p;
         cNames = cNames || {"Pr>|Z|"};
         Crit = quantile("Normal", 1 - alpha/2);      
      end;
      /* Add confidence intervals to table */
      lower = estimate - Crit # stdErr;
      upper = estimate + Crit # stdErr;
      PE = PE || lower || upper;
      clPct = char(100*(1-alpha), 6, 1);
      lowName = "Lower " + strip(clPct) + "% CL";
      upName  = "Upper " + strip(clPct) + "% CL";
      cNames = cNames || (lowName || upName);
   end;

   /* Print Parameter Estimates table */
   lablText = "Parameter Estimates for " + dist + " Distribution";
   print PE [label=lablText rowname=parmNames colname=cNames format=8.4];

   /* printOpt >= 1: Fit Criteria */
   if printOpt >= 1 then do;
      crit = L$"Crit";
      critNames = L$"CritNames";
      print crit [label="Fit Criteria" rowname=critNames c="Value" format=8.4];
   end;

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
         call PrintToLog("Hessian is negative definite. Local maximum verified.", 0);
      else
         call PrintToLog("Hessian is not negative definite. Solution is not a local maximum.", 1);
   end;
   
finish;

store module=(
      MLE_Summary
);
