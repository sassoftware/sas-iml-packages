/* BEFORE running this example, store the modules in the package, as shown in Install_Pkg.sas */

/* Examples in the documentation for the MLE package */
proc iml;
load module=_all_;     /* load the MLE library */

/* Example 1: Define limits and covariance matrix */
b = {1 4 2};
Sigma = {1.0 0.6 0.3333333333,
         0.6 1.0 0.7333333333,
         0.3333333333 0.7333333333 1.0 };
prob = cdfmvn(b, Sigma);
print prob;

/* Example 2: 3-D example where X1 is independent of (X2, X3) */
b = {0.5 1.0 1.5};
Sigma = {1.0 0.0 0.0,
         0.0 1.0 0.5,
         0.0 0.5 1.0};
prob = cdfmvn(b, Sigma);
/* Validation: Phi(0.5) * Phi2(1.0, 1.5, 0.5) */
check = cdf("Normal", 0.5) * probbnrm(1.0, 1.5, 0.5);
print prob check;

/* Example 3: 5-D example with non-zero mean and non-diagonal covariance. */
Sigma = {1 1 1 1 1, 
        1 2 2 2 2, 
        1 2 3 3 3, 
        1 2 3 4 4, 
        1 2 3 4 5};
b = 0:4;
prob_centered = cdfmvn(b, Sigma);

mu = {2 2 2 2 2};
prob_mu = cdfmvn(b, Sigma, mu);
print prob_centered prob_mu;

QUIT;

