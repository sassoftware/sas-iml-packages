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

- **CDFMVN**: The main function for estimating the CDF of a multivariate normal random variable, 
X ~ MVN(mu, Sigma), where mu is a k-dimensional row vector and Sigma is a kxk covariance matrix.
If b = (b_1, b_2, ..., b_k) is the upper limit of integration, the function returns the probability that the random variable is in the left-tailed region {X_1 < b_1 & X_2 < b_2 & ... & X_k < b_k}.

In future releases, the package will support additional functions for probabilities of muiltivariate distributions.

## Example

```sas
proc iml;
load module=_all_;     /* load the library */

/* Example 1: Define limits and covariance matrix */
b = {1 4 2};
Sigma = {1.0 0.6 0.3333333333,
         0.6 1.0 0.7333333333,
         0.3333333333 0.7333333333 1.0 };
prob = cdfmvn(b, Sigma);
print prob;
QUIT;
```