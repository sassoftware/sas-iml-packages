/*********************************************************/
/* MLE_LL.sas                                            */
/* Let p be a valid set of parameters for a built-in     */
/* distribution. The functions in this file return the   */
/* value of the loglikelihood function at p: LL(p)       */
/*********************************************************/

/* Define LL functions for common distributions */
start lik_LL_Beta(param) global(gMLE_y); 
   alpha = param[1]; 
   beta = param[2];
   LL = sum( logpdf("Beta", gMLE_y, alpha, beta) );
   return( LL );
finish;
start lik_ValidData_Beta(y); 
   return all(y > 0 & y < 1);
finish;

start lik_LL_Expo(param) global(gMLE_y); 
   sigma = param[1]; 
   LL = sum( logpdf("Expo", gMLE_y, sigma) );
   return( LL );
finish;
start lik_ValidData_Expo(y); 
   return all(y >= 0 );
finish;

start lik_LL_Gamma(_param) global(gMLE_y);
   LL = .;
   param = colvec(_param);
   alpha = param[1]; 
   if nrow(param) = 1 then 
      LL = sum( logpdf("Gamma", gMLE_y, alpha) );
   else if nrow(param) = 2 then do;
      lambda = param[2];
      LL = sum( logpdf("Gamma", gMLE_y, alpha, lambda) );
   end;
   return( LL );
finish;
start lik_ValidData_Gamma(y); 
   return all(y >= 0 );
finish;

/* Log likelihood for Gumbel distribution:
       f(x) = (1/sigma) * exp(-(x-mu)/sigma) * exp(-exp(-(x-mu)/sigma))
   log f(x) = -log(sigma) - (x-mu)/sigma - exp(-(x-mu)/sigma)
            = -n*log(sigma) - sum(z) - sum(exp(-z)), where z=(x-mu)/sigma
*/
start lik_LL_Gumbel(_param) global(gMLE_y);
   LL = .;
   param = colvec(_param);
   if nrow(param) = 2 then do;
      mu = param[1];
      sigma = param[2];
      z = (gMLE_y-mu)/sigma;
      LL = -nrow(gMLE_y) * log(sigma) - sum(z) - sum(exp(-z));
   end;
   return( LL );
finish;
start lik_ValidData_Gumbel(y); 
   return 1; /* all real numbers are valid */
finish;

start lik_LL_IGauss(param) global(gMLE_y); /* Inverse Gaussian */
   lambda = param[1];
   mu = param[2];
   LL = sum( logpdf("IGAUSS", gMLE_y, lambda, mu) );
   return( LL );
finish;
start lik_ValidData_IGauss(y); 
   return all(y > 0 );
finish;

start lik_LL_LN2(param) global(gMLE_y); 
   mu = param[1]; 
   sigma = param[2];
   LL = sum( logpdf("Lognormal", gMLE_y, mu, sigma) );
   return( LL );
finish;
start lik_ValidData_LN2(y); 
   return all(y > 0 );
finish;

start lik_LL_Normal(param) global(gMLE_y); 
   mu = param[1]; 
   sigma = param[2];
   LL = sum( logpdf("Normal", gMLE_y, mu, sigma) );
   return( LL );
finish;
start lik_ValidData_Normal(y); 
   return 1; /* all real numbers are valid */
finish;

start lik_LL_Weib2(param) global(gMLE_y); 
   c = param[1];
   lambda = param[2];
   LL = sum( logpdf("Weibull", gMLE_y, c, lambda) );
   return( LL );
finish;
start lik_ValidData_Weib2(y); 
   return all(y >= 0 );
finish;

/*** direct top-level API ***/

start MLE_LL(distname, param);
   /* construct name of function */
   func_name = lik_func_name(distname, "LL");
   if missing(func_name) then 
      func_name = distname;   /* not a built-in function; call directly */
   %EVALFUNC1(LL, func_name, param);
   return( LL );
finish;
start MLE_LL_ValidData(distname, y);
   /* construct name of function */
   func_name = lik_func_name(distname, "ValidData");
   if missing(func_name) then 
      func_name = distname;   /* not a built-in function; call directly */
   %EVALFUNC1(isValid, func_name, y);
   return( isValid );
finish;

start MLE_LL_and_Deriv(f, g, H, distname, param);
   /* construct name of function */
   func_name = lik_func_name(distname, "LL");
   call nlpfdd(f, g, H, func_name, param);
finish;

/* Return PDFs of the various distributions evaluated on a uniform grid of X values:
   INPUT:
   DistNames : a (k x 1) vector of distribution keyword names (eg, 'Normal')
   params    : a (k x 3) matrix of parameters, one row for each distribution
   x         : an (N x 1) vector of values in the support of the distribution. 
               Each PDF is evaluated at the points in x, which is usually on a regular grid.
   OUTPUT:
   Labl: A vector of groups of the form DIST(parm1,parm2)  (eg, 'Normal(1.23,4.56)')
   P   : An (m x k) matrix of PDFs. The i_th column is the PDF at x for the i_th distribution
*/
start MLE_PDF(Labl, P, DistNames, params, x);
   Format = "BEST6.";
   k = nrow(DistNames);
   N = nrow(x);
   P = j(N,k,.);
   Labl = j(k,1,BlankStr(32));
   do i = 1 to k;
      name = lik_dist_name(DistNames[i]);
      parms = params[i,];
      f = putn(parms,Format);
      nParms = ncol(loc(parms^=.));
      if nParms=1 then do;
         labl[i] = cats(name,'(',f[1],')');
         P[,i] = PDF(name, x, parms[1]);
      end;
      else if nParms=2 then do;
         labl[i] = cats(name,'(',f[1],',',f[2],')');
         P[,i] = PDF(name, x, parms[1], parms[2]);
      end;
      else if nParms=3 then do;
         labl[i] = cats(name,'(',f[1],',',f[2],',',f[3],')');
         P[,i] = PDF(name, x, parms[1], parms[2], parms[3]);
      end;
  end;
finish;

store module=(
      lik_LL_Beta
      lik_LL_Expo
      lik_LL_Gamma
      lik_LL_Gumbel
      lik_LL_IGauss
      lik_LL_LN2
      lik_LL_Normal
      lik_LL_Weib2
      MLE_LL
      lik_ValidData_Beta
      lik_ValidData_Expo
      lik_ValidData_Gamma
      lik_ValidData_Gumbel
      lik_ValidData_IGauss
      lik_ValidData_LN2
      lik_ValidData_Normal
      lik_ValidData_Weib2
      MLE_LL_ValidData      
      MLE_LL_and_Deriv
      MLE_PDF
);
