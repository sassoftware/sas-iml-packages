proc iml;
%include "&MLE_path./MLE_define.sas";
QUIT;


/*******************************************/
/* LL tests */
/*******************************************/

PROC IML;
load module=_ALL_;
print "--- A successfult test prints only 'TEST DONE' ---";

/* Helper function to check for equality between two character matrices */
start NumMatEqual(A, B, tol=1e-3);
   if type(A)^='N' | type(B)^='N' |
      (nrow(A)^=nrow(B)) | (ncol(A)^=ncol(B)) then 
      return(0);
   m = max(abs(A-B));
   *print m;
   return( m < tol);
finish;

Correct_Beta_LL = 1080.6886;
Correct_Expo_LL = -2077.525;
Correct_Gamma_LL = -2708.222;
Correct_Gumbel_LL = -2265.088;
Correct_IGauss_LL = -2255.112;
Correct_Normal_LL = -2128.205;
Correct_lognormal_LL = -3712.441;
Correct_Weibull_LL = -1457.625;

call randseed(12345);
/* Y ~ Beta(alpha, beta) */
N = 1000;
alpha = 3;
beta = 0.5;
Param = alpha // beta;
Y = randfun(N, "Beta", alpha, beta);
isValid = MLE_Init(Y, "beta");
if ^isValid then print "ERROR in MLE_Init for Beta";
ll = MLE_LL("beta", Param);
*print ll[L="Beta LL"];
if ^NumMatEqual(ll, Correct_Beta_LL) then 
   print "ERROR in Beta LL";
else;

/* Y ~ Expo(scale) */
N = 1000;
scale = 3;
Param = scale;
Y = randfun(N, "Expo", scale);
isValid = MLE_Init(Y, "Expo");
if ^isValid then print "ERROR in MLE_Init for Expo";
ll = MLE_LL("Expo", scale);
*print ll[L="Expo LL"];
if ^NumMatEqual(ll, Correct_Expo_LL) then 
   print "ERROR in Expo LL";
else;

/* Y ~ Gamma(alpha, lambda) */
N = 1000;
alpha = 4;
lambda = 2;
Param = alpha // lambda;
Y = randfun(N, "Gamma", alpha, lambda);
isValid = MLE_Init(Y, "Gamma");
if ^isValid then print "ERROR in MLE_Init for Gamma";
ll = MLE_LL("Gamma", Param);
*print ll[L="Gamma LL"];
if ^NumMatEqual(ll, Correct_Gamma_LL) then 
   print "ERROR in Gamma LL";
else;

/* Y ~ Gumbel(mu, sigma) */
N = 1000;
mu = 4;
sigma = 2;
Param = mu // sigma;
Y = randfun(N, "Gumbel", mu, sigma);
isValid = MLE_Init(Y, "Gumbel");
if ^isValid then print "ERROR in MLE_Init for Gumbel";
ll = MLE_LL("Gumbel", Param);
*print ll[L="Gumbel LL"];
if ^NumMatEqual(ll, Correct_Gumbel_LL) then 
   print "ERROR in Gumbel LL";
else;

/* Y ~ IGauss(lambda, mu) */
N = 1000;
mu = 4;
sigma = lambda;
Param = sigma // mu;
Y = randfun(N, "IGauss", lambda, mu);
isValid = MLE_Init(Y, "IGauss");
if ^isValid then print "ERROR in MLE_Init for IGauss";
ll = MLE_LL("IGauss", Param);
*print ll[L="IGauss LL"];
if ^NumMatEqual(ll, Correct_IGauss_LL) then 
   print "ERROR in IGauss LL";
else;

/* Y ~ Normal(mu, sigma) */
N = 1000;
mu = 10;
sigma = 2;
Param = mu // sigma;
Y = randfun(N, "Normal", mu, sigma);
isValid = MLE_Init(Y, "Normal");
if ^isValid then print "ERROR in MLE_Init for Normal";
ll = MLE_LL("Normal", Param);
*print ll[L="Normal LL"];
if ^NumMatEqual(ll, Correct_Normal_LL) then 
   print "ERROR in Normal LL";
else;

/* Y ~ lognormal(mu, sigma) */
N = 1000;
mu = 3;
sigma = 0.5;
Param = mu // sigma;
Y = randfun(N, "lognormal", mu, sigma);
isValid = MLE_Init(Y, "lognormal");
if ^isValid then print "ERROR in MLE_Init for Lognormal";
ll = MLE_LL("lognormal", Param);
*print ll[L="lognormal LL"];
if ^NumMatEqual(ll, Correct_lognormal_LL) then 
   print "ERROR in Lognormal LL";
else;


/* Y ~ Weibull(c, lambda) */
N = 1000;
c = 1.5;
lambda = 2;
Param = c // lambda;
Y = randfun(N, "Weibull", c, lambda);
isValid = MLE_Init(Y, "Weibull");
if ^isValid then print "ERROR in MLE_Init for Weibull";
ll = MLE_LL("Weibull", Param);
*print ll[L="Weibull LL"];
if ^NumMatEqual(ll, Correct_Weibull_LL) then 
   print "ERROR in Weibull LL";
else;

print "--- TEST DONE ---";

*QUIT;
