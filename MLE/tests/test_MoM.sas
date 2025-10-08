proc iml;
%include "&MLE_path./MLE_Util.sas";
%include "&MLE_path./MLE_Keywords.sas";
%include "&MLE_path./MLE_MoM.sas";
QUIT;


/*******************************************/
/* MoM tests */
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
est = MLE_MOM("Beta", Y);
est2 = lik_MOM_Beta(Y);
diff_est = est - est2;
print Param[r={'alpha' 'beta'}] est est2 diff_est;

/* Y ~ Expo(scale) */
N = 1000;
scale = 3;
Param = scale;
Y = randfun(N, "Expo", scale);
est = MLE_MOM("Exp", Y);
est2 = lik_MOM_Expo(Y);
diff_est = est - est2;
print Param[r={'scale'}] est est2 diff_est;

/* Y ~ Gamma(alpha, lambda) */
N = 1000;
alpha = 4;
lambda = 2;
Param = alpha // lambda;
Y = randfun(N, "Gamma", alpha, lambda);
est = MLE_MOM("Gamma", Y);
est2 = lik_MOM_Gamma(Y);
diff_est = est - est2;
print Param[r={'alpha' 'lambda'}] est est2 diff_est;

/* Y ~ Gumbel(mu, sigma) */
N = 1000;
mu = 4;
sigma = 2;
Param = mu // sigma;
Y = randfun(N, "Gumbel", mu, sigma);
est = MLE_MOM("Gumbel", Y);
est2 = lik_MOM_Gumbel(Y);
diff_est = est - est2;
print Param[r={'mu' 'sigma'}] est est2 diff_est;

/* Y ~ IGauss(lambda, mu) */
N = 1000;
mu = 4;
sigma = lambda;
Param = sigma // mu;
Y = randfun(N, "IGauss", lambda, mu);
est = MLE_MOM("IGauss",Y);
est2 = lik_MOM_IGauss(Y);
diff_est = est - est2;
print Param[r={'sigma' 'mu'}] est est2 diff_est;

/* Y ~ Normal(mu, sigma) */
N = 1000;
mu = 10;
sigma = 2;
Param = mu // sigma;
Y = randfun(N, "Normal", mu, sigma);
est = MLE_MOM("Normal", Y);
est2 = lik_MOM_Normal(Y);
diff_est = est - est2;
print Param[r={'mu' 'sigma'}] est est2 diff_est;

/* Y ~ lognormal(mu, sigma) */
N = 1000;
mu = 3;
sigma = 0.5;
Param = mu // sigma;
Y = randfun(N, "lognormal", mu, sigma);
est = MLE_MOM("lognormal", Y);
est2 = lik_MOM_LN2(Y);
diff_est = est - est2;
print Param[r={'mu' 'sigma'}] est est2 diff_est;


/* Y ~ Weibull(c, lambda) */
N = 1000;
c = 1.5;
lambda = 2;
Param = c // lambda;
Y = randfun(N, "Weibull", c, lambda);
est = MLE_MOM("Weibull", Y);
est2 = lik_MOM_Weib2(Y);
diff_est = est - est2;
print Param[r={'c' 'lambda'}] est est2 diff_est;

/* Weibull root equation */
cc = {0.02,0.05,0.1,0.2,0.5,1,2,3,6,10};
s = lik_MoM_Weib2_Root(cc);
print cc s;


QUIT;
