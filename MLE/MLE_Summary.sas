/***********************************/
/* Define printing functions */
/***********************************/


start MLE_Summary(L, printOpt=1);
   labl = "Fit for " + L$"Dist" + " distribution";
   T = L$"Estimate" || L$"StdErr";
   print T [label=labl colname=L$"ParmNames" format=8.4 
            rowname={"Estimate" "StdErr"}];

   if printOpt >= 1 then 
      print L$"Crit" [label="Fit Criteria" colname=L$"CritNames" format=8.4];

   if printOpt >= 2 then do;
      print L$"LL" [label="Log-Likelihood at Estimate"];
      print L$"Grad" [label="Gradient at Estimate" colname=L$"ParmNames" rowname="Grad"];
      print L$"Hessian" [label="Hessian at Estimate" colname=L$"ParmNames" rowname=L$"ParmNames"];
   end;  
finish;

store module=(
      MLE_Summary
);