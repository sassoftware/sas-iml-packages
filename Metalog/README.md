# The Metalog Package

## Overview

The metalog family of distributions (Keelin 2016)
is a system of flexible distributions that can model a
wide range of univariate data distributions.  The **Metalog** package provides a library of SAS® IML functions that can model univariate data by using the metalog distributions.

The metalog system and the SAS IML functions are described in the documentation of the **Metalog** package.  To download the files for the package, see the [README file](../README.md) for the sas-iml-packages repository. Instructions are also provided in Appendix A of the **Metalog** package documentation. 


## Getting Started

Before you can use the SAS IML functions that implement the Metalog distribution, you must read the definitions into the current SAS® session. In PROC IML, you can use the %INCLUDE statement to read the files. For example, the program that follows performs these steps in PROC IML:

1. Use %INCLUDE statements to define the function in the **Metalog** package.
2. Fit a five-term unbounded metalog model to the data.
3. Visualize the model and the empirical distribution of the data.
4. Display the 5th, 25th, 75th, and 95th percentiles of the model.
5. Generate a random sample of size 10 from the model.

```sas
%let repoPath = C:/sas-iml-packages;          /* set the path */
proc iml;
/* define the metalog functions */
%include "&repoPath./Metalog/ML_proc.sas";    /* required for PROC IML */
%include "&repoPath./Metalog/ML_define.sas";  /* required for procedure or action */

/* fit a five-term unbounded metalog model to data */
x  = {14, 17.3, 20, 22, 24, 27, 30, 33, 37.5};
order = 5;
MLObj = ML_CreateFromData(x, order);
/* visualize the model and the empirical distribution of the data */
title "Metalog Model and ECDF";
run ML_PlotECDF(MLObj);
/* display the 5th, 25th, 75th, and 95th percentiles of the model */
p = {0.05, 0.25, 0.75, 0.95};
q = ML_Quantile(MLObj, p);
print p[F=PERCENT9.] q[F=6.2]; 
/* generate a random sample of size 10 from the model */
call randseed(12345);
Z = ML_Rand(MLObj, 10);
print Z;
quit;
```

When you run a SAS IML program that includes the definition files, the function definitions are stored to the active library. If you want to use the function again in the same SAS session, use the LOAD statement to reload the definitions, as follows:

```sas
proc iml;
/* define the metalog functions */
load module= _all_;
x  = {14, 17.3, 20, 22, 24, 27, 30, 33, 37.5};
MLObj = ML_CreateFromData(x);
run ML_Summary(MLObj);
quit;
```

## Additional Resources

- Keelin, T. (2016). "The Metalog Distributions." *Decision Analysis* 13:243–277. https://doi.org/10.1287/deca.2016.0338.
- Keelin, T. (2016). "The Metalog Distributions." http://www.metalogdistributions.com.

