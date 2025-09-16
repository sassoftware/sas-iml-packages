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
   /* TO DO */
   run MLE_Init(y);
   /* TEMP: Return MoM estimate until optimization is implemented */
   est = MLE_MoM(DistName, y);
   return( est );
finish;



start MLE_Fit(DistName, y,  param0=, OptimMethod=, Bounds=);
   /* need to check how these optional parameters are passed into MLE() */
   est = MLE(DistName, y, param0, OptimMethod, Bounds); 
   names = {"Dist", "y", "ParmNames", "Estimate", "LL", 
            "Grad", "Hessian", "StdErr", "Crit", "CritNames"};
   FitObj = ListCreate(names);  /* named list */
   FitObj$"Dist" = lik_dist_keyword(DistName, 0);
   FitObj$"y" = gMLE_y;
   FitObj$"ParmNames" = lik_dist_parmnames(DistName);
   FitObj$"Estimate" = est;
   run MLE_LL_and_Deriv(LL, Grad, Hess, DistName, est);
   FitObj$"LL" = LL;
   FitObj$"Grad" = Grad;
   FitObj$"Hessian" = Hess;
   /* TO DO: compute StdErr and Crit */
   /*
   FitObj$"StdErr" = StdErr;
   FitObj$"Crit" = Crit;
   FitObj$"CritNames" = CritNames;
   */
   FitObj$"StdErr" = j(nrow(est), 1, .);
   FitObj$"CritNames" = {'-2*LL', 'AIC', 'AICC', 'BIC'};
   FitObj$"Crit" = -2*LL // . // . // .;
   return( FitObj );
finish;

store module=(
      MLE
      MLE_Fit
);