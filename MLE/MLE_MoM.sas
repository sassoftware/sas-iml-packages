/*********************************************************/
/* MLE_MoM.sas                                           */
/* Define functions that return the parameters estimates */
/* for built-in distributions (or some other suitable    */
/* guess for an MLE optimization to use)                 */
/*********************************************************/

/* Method of Moments (MoM) for common distributions */
start lik_MoM_Beta(_y);
   y = colvec(_y);
   m = mean(y);
   s = std(y);
   p = m*(1-m)/s##2 - 1;
   alpha = m*p;
   beta = (1-m)*p;
   return( alpha // beta );
finish;

start lik_MoM_Expo(_y);
   y = colvec(_y);
   scale = mean(y);
   return( scale );
finish;

start lik_MoM_Gamma(_y);
   y = colvec(_y);
   m = mean(y);
   s2 = var(y);
   lambda = s2 / m; 
   alpha = m##2 / s2;
   return( alpha // lambda );
finish;

start lik_MoM_Gumbel(_y);
   pi = constant('pi');
   gamma = constant('Euler');
   y = colvec(_y);
   m = mean(y);
   s = std(y);
   sigma = sqrt(6)/pi * s;
   mu = m - gamma*sigma;
   return( mu // sigma );
finish;

start lik_MoM_IGauss(_y);
   y = colvec(_y);
   n = nrow(y);
   m = mean(y);
   s2_mle = var(y) * ((n-1)/n);  /* sample variance w/o Bessel correction */
   mu = m;
   lambda = m##3 / s2_mle;
   return( lambda // mu );
finish;

start lik_MoM_LN2(_y);
   y = colvec(_y);
   m1 = mean(y);
   m2 = ssq(y) / nrow(y);
   mu = log( m1##2 / sqrt(m2) );
   sigma = sqrt( log(m2 / m1##2) );
   return( mu // sigma );
finish;

start lik_MoM_Normal(_y);
   y = colvec(_y);
   mu = mean(y);
   sigma = std(y);
   return( mu // sigma );
finish;

/* There are several Method of Moments for 2-param Weibull, including 
   some in Johnson, Kotz, Balakrishnan (1994, 2nd Ed)
   "Continuous Univariate Distributions, Volume 1."
   This method follow the blog post by Charles Zaiontz at
   https://real-statistics.com/distribution-fitting/method-of-moments/method-of-moments-weibull/
   which is based on 
   Garcia, O. (1981) Simplified method-of-moments estimation of the Weibull distribution
   https://www.scionresearch.com/__data/assets/pdf_file/0003/59286/NZJFS1131981GARCIA304-306.pdf
   The Weib2 distribution assumes threshold=0 and all Y > 0.

   The names for the parameters vary:
   SOURCE      SHAPE SCALE THRESHOLD
   ------      ----- ----- ---------
   PDF doc       a   lambda    -
 =>UNIVARIATE    c   sigma   theta
   RAND          a   b         -
   Wikipedia     k   lambda  theta
   Zaiontz site beta alpha
*/
/* Here: c is shape parameter and sigma is scale.
   The Gamma function grows rapidly, so instead of finding root of function 
   that involves Gamma functions, find root of the log of the function.
   */
start lik_MoM_Weib2_Root(c) global(gMLE_y);
   m = mean(gMLE_y);
   s2 = var(gMLE_y);
   f = 2*log(m) - log(s2 + m##2) - 2*lgamma(1 + 1/c) + lgamma(1 + 2/c);
   return( f );
finish;

/* TO DO: improve this hard-coded search interval */
start lik_MoM_Weib2(_y);
   y = colvec(_y);
   m = mean(y);
   c = froot("lik_MoM_Weib2_Root", {0.01 10});
   sigma = m / GAMMA(1 + 1/c);
   return( c // sigma );
finish;

/*** direct top-level API ***/

start MLE_MoM(distname, _y) global(G_DEBUG);
   IF G_DEBUG THEN run PrintLoc();   
   isValid = MLE_Init(_y, DistName);           /* creates GLOBAL var */
   if ^isValid then
      return( . );
   /* construct name of function */
   func_name = lik_func_name(distname, "MoM");
   if missing(func_name) then 
      func_name = distname;   /* not a built-in function; call directly */
   y = lik_GetValidatedData(_y);
   %EVALFUNC1( MoM_est, func_name, y);
   return( MoM_est );
finish;

store module=(
      lik_MoM_Beta
      lik_MoM_Expo
      lik_MoM_Gamma
      lik_MoM_Gumbel
      lik_MoM_IGauss
      lik_MoM_LN2
      lik_MoM_Normal
      lik_MoM_Weib2_Root
      lik_MoM_Weib2
      MLE_MoM
);