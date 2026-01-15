/* BEFORE running this example, store the modules in the MLE package, as shown in test_Install.sas */

/*******************************************/
/* MoM tests */
/*******************************************/

PROC IML;
load module=_ALL_;

print "--- MLE_MoM: A successfult test prints multiple tables and the text 'TEST DONE' ---";
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
results = Param || est || est2 || diff_est;
cnames = {'Param' 'Estimate' 'Est2' 'diff'};
print results[r={'alpha' 'beta'} c = cnames L="Beta Estimates"];

correct = {
3 2.8009079 2.8009079 0 ,
0.5 0.4894146 0.4894146 0 
};
if ^all(max(abs(results-correct)<1e-6)) then 
   print "ERROR: Beta results are incorrect";

/* Y ~ Expo(scale) */
N = 1000;
scale = 3;
Param = scale;
Y = randfun(N, "Expo", scale);
est = MLE_MOM("Exp", Y);
est2 = lik_MOM_Expo(Y);
diff_est = est - est2;
results = Param || est || est2 || diff_est;
cnames = {'Param' 'Estimate' 'Est2' 'diff'};
print results[r={'scale'} c=cnames L="Expo Estimates"];

correct = {
3 2.9367369 2.9367369 0 
};
if ^all(max(abs(results-correct)<1e-6)) then 
   print "ERROR: Expo results are incorrect";


/* Y ~ Gamma(alpha, lambda) */
N = 1000;
alpha = 4;
lambda = 2;
Param = alpha // lambda;
Y = randfun(N, "Gamma", alpha, lambda);
est = MLE_MOM("Gamma", Y);
est2 = lik_MOM_Gamma(Y);
diff_est = est - est2;
results = Param || est || est2 || diff_est;
cnames = {'Param' 'Estimate' 'Est2' 'diff'};
print results[r={'alpha' 'lambda'} c=cnames L="Gamma Estimates"];

correct = {
4 4.0594573 4.0594573 0 ,
2 1.9570362 1.9570362 0 
};
if ^all(max(abs(results-correct)<1e-6)) then 
   print "ERROR: Gamma results are incorrect";

/* Y ~ Gumbel(mu, sigma) */
N = 1000;
mu = 3;
sigma = 2.2;
Param = mu // sigma;
Y = randfun(N, "Gumbel", mu, sigma);
est = MLE_MOM("Gumbel", Y);
est2 = lik_MOM_Gumbel(Y);
diff_est = est - est2;
results = Param || est || est2 || diff_est;
cnames = {'Param' 'Estimate' 'Est2' 'diff'};
print results[r={'mu' 'sigma'} c=cnames L="Gumbel Estimates"];

correct = {
3 3.0683966 3.0683966 0 ,
2.2 2.1508534 2.1508534 0 
};
if ^all(max(abs(results-correct)<1e-6)) then 
   print "ERROR: Gumbel results are incorrect";


/* Y ~ IGauss(lambda, mu) */
N = 1000;
mu = 4;
sigma = 1.8;
Param = sigma // mu;
Y = randfun(N, "IGauss", lambda, mu);
est = MLE_MOM("IGauss",Y);
est2 = lik_MOM_IGauss(Y);
diff_est = est - est2;
results = Param || est || est2 || diff_est;
cnames = {'Param' 'Estimate' 'Est2' 'diff'};
print results[r={'mu' 'sigma'} c=cnames L="IGauss Estimates"];

correct = {
1.8 2.0879949 2.0879949 0 ,
4 4.0226402 4.0226402 0 
};
if ^all(max(abs(results-correct)<1e-6)) then 
   print "ERROR: IGauss results are incorrect";


/* Y ~ Normal(mu, sigma) */
N = 1000;
mu = 10;
sigma = 2;
Param = mu // sigma;
Y = randfun(N, "Normal", mu, sigma);
est = MLE_MOM("Normal", Y);
est2 = lik_MOM_Normal(Y);
diff_est = est - est2;
results = Param || est || est2 || diff_est;
cnames = {'Param' 'Estimate' 'Est2' 'diff'};
print results[r={'mu' 'sigma'} c=cnames L="Normal Estimates"];

correct = {
10 9.9244503 9.9244503 0 ,
2 2.0315931 2.0315931 0 
};
if ^all(max(abs(results-correct)<1e-6)) then 
   print "ERROR: Normal results are incorrect";


/* Y ~ lognormal(mu, sigma) */
N = 1000;
mu = 3;
sigma = 0.5;
Param = mu // sigma;
Y = randfun(N, "lognormal", mu, sigma);
est = MLE_MOM("lognormal", Y);
est2 = lik_MOM_LN2(Y);
diff_est = est - est2;
results = Param || est || est2 || diff_est;
cnames = {'Param' 'Estimate' 'Est2' 'diff'};
print results[r={'mu' 'sigma'} c=cnames L="Lognormal Estimates"];

correct = {
3 3.0187233 3.0187233 0 ,
0.5 0.4769941 0.4769941 0 
};
if ^all(max(abs(results-correct)<1e-6)) then 
   print "ERROR: Lognormal results are incorrect";


/* Y ~ Weibull(c, lambda) */
N = 1000;
c = 1.5;
lambda = 2;
Param = c // lambda;
Y = randfun(N, "Weibull", c, lambda);
est = MLE_MOM("Weibull", Y);
est2 = lik_MOM_Weib2(Y);
diff_est = est - est2;
results = Param || est || est2 || diff_est;
cnames = {'Param' 'Estimate' 'Est2' 'diff'};
print results[r={'c' 'lambda'} c=cnames L="Weibull Estimates"];

correct = {
1.5 1.5181732 1.5181732 0, 
2 1.9761111 1.9761111 0 
};
if ^all(max(abs(results-correct)<1e-6)) then 
   print "ERROR: Weibull results are incorrect";

/* Weibull root equation */
c = {0.02,0.05,0.1,0.2,0.5,1,2,3,6,10};
s = lik_MoM_Weib2_Root(c);

correct = {
66.411667 ,
25.277233 ,
11.754617 ,
5.1572548 ,
1.4195852 ,
0.3209729 ,
-0.13061 ,
-0.248106 ,
-0.335314 ,
-0.357803 
};
if ^all(max(abs(s-correct)<1e-6)) then 
   print "ERROR: Weib2_root results are incorrect";

print "--- TEST DONE ---";

QUIT;
