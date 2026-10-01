# The PROBMVN repo

## Description

This project is a library of SAS IML functions for estimating high-dimensional 
probabilities for multivariate normal distributions. 
In SAS Viya, these functions are implemented as built-in functions. These IML modules are provided for customers who are using SAS 9.4.

## Documentation

The SAS IML functions are described in the documentation for the  package. The file probmvn.pdf is a PDF file 
that describes the syntax of each public function. The documentation shows how to call each top-level public 
function and provides examples of each function's output.

## Main functions

The following high-level functions are designed to be called directly:

- **CDFMVN_MOD**: The main function for estimating the CDF of a multivariate normal random variable, 
X ~ MVN(mu, Sigma), where mu is a k-dimensional row vector, and Sigma is a kxk covariance matrix.
If U = (U_1, U_2, ..., U_k) is the upper limit of integration, the function returns the probability that the random variable is in the left-tailed region {X_1 < U_1 & X_2 < U_2 & ... & X_k < U_k}.
In this implementation, 2 <= k <= 20.

- **PROBMVN_MOD**: The main function for estimating the probability of a multivariate normal random variable on a rectangular domain. The random variable is 
X ~ MVN(mu, Sigma), where mu is a k-dimensional row vector, and Sigma is a kxk covariance matrix.
If L = (L_1, L_2, ..., L_k) is the lower limit of integration, and
U = (U_1, U_2, ..., U_k) is the upper limit of integration, 
the function returns the probability that the random variable is in the rectangular region {L_1 < X_1 < U_1 & L_2 < X_2 < U_2 & ... & L_k < X_k < U_k}. If you want a lower limit to be -Infinity, use the special missing value .M in the L vector. If you want an upper limit to be +Infinity, 
use the special missing value .I in the U vector. 
In this implementation, 2 <= k <= 100.

## Example

```sas
proc iml;
load module=_all_;     /* load the library */

/* Example 1: Define limits and covariance matrix */
U = {1 4 2};
Sigma = {1.0 0.6 0.3333333333,
         0.6 1.0 0.7333333333,
         0.3333333333 0.7333333333 1.0 };
CDF = cdfmvn_mod(U, Sigma);
print CDF;

L = {-1 0 -2};
prob_rect = probmvn_mod(L, U, Sigma);
print prob_rect;
QUIT;
```