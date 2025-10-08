proc iml;
%include "&MLE_path./MLE_Util.sas";
%include "&MLE_path./MLE_Keywords.sas";
%include "&MLE_path./MLE_MoM.sas";
QUIT;


/*******************************************/
/* LL tests */
/*******************************************/

PROC IML;
load module=_ALL_;

call randseed(12345);
/* Y ~ Beta(alpha, beta) */
N = 1000;
alpha = 3;
beta = 0.5;
Param = alpha // beta;
Y = randfun(N, "Beta", alpha, beta);
run MLE_Init(Y);
ll = MLE_LL("beta", Param);
print ll[L="Beta LL"];

/* Y ~ Expo(scale) */
N = 1000;
scale = 3;
Param = scale;
Y = randfun(N, "Expo", scale);
run MLE_Init(Y);
ll = MLE_LL("Expo", scale);
print ll[L="Expo LL"];

/* Y ~ Gamma(alpha, lambda) */
N = 1000;
alpha = 4;
lambda = 2;
Param = alpha // lambda;
Y = randfun(N, "Gamma", alpha, lambda);
run MLE_Init(Y);
ll = MLE_LL("Gamma", Param);
print ll[L="Gamma LL"];

/* Y ~ Gumbel(mu, sigma) */
N = 1000;
mu = 4;
sigma = 2;
Param = mu // sigma;
Y = randfun(N, "Gumbel", mu, sigma);
run MLE_Init(Y);
ll = MLE_LL("Gumbel", Param);
print ll[L="Gumbel LL"];

/* Y ~ IGauss(lambda, mu) */
N = 1000;
mu = 4;
sigma = lambda;
Param = sigma // mu;
Y = randfun(N, "IGauss", lambda, mu);
run MLE_Init(Y);
ll = MLE_LL("IGauss", Param);
print ll[L="IGauss LL"];

/* Y ~ Normal(mu, sigma) */
N = 1000;
mu = 10;
sigma = 2;
Param = mu // sigma;
Y = randfun(N, "Normal", mu, sigma);
run MLE_Init(Y);
ll = MLE_LL("Normal", Param);
print ll[L="Normal LL"];

/* Y ~ lognormal(mu, sigma) */
N = 1000;
mu = 3;
sigma = 0.5;
Param = mu // sigma;
Y = randfun(N, "lognormal", mu, sigma);
run MLE_Init(Y);
ll = MLE_LL("lognormal", Param);
print ll[L="lognormal LL"];


/* Y ~ Weibull(c, lambda) */
N = 1000;
c = 1.5;
lambda = 2;
Param = c // lambda;
Y = randfun(N, "Weibull", c, lambda);
run MLE_Init(Y);
ll = MLE_LL("Weibull", Param);
print ll[L="Weibull LL"];



QUIT;
