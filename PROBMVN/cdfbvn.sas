/* ------------------------------------------------------------------
   Main Function: cdfbvn_mod
   ------------------------------------------------------------------ */
/* Return the bivariate CDF for MVN(Sigma, mu) at each row of b.
   Validate the parameters, standardize to correlation scale, and call PROBBNRM */
proc iml;
start cdfbvn_mod(b, Sigma, mu={0 0});
   IsValid = mvn_IsValidParmsMVN(b, Sigma, mu);
   if ^IsValid then
      return(.);
   run mvn_StdizeCovToCorr(U, R, b, Sigma, mu);
   prob = probbnrm(U[,1], U[,2], R[1,2]);
   return ( prob );
finish;
store module=(cdfbvn_mod);
quit;
