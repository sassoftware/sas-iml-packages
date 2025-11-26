# The MLE repo

# THIS PACKAGE IS UNDER DEVELOPMENT. DO NOT USE.

## Description

This project is a library of SAS IML functions for fitting univariate distributions to data by using maximum likelihood estimation (MLE).
The functions serve several purposes:
- Show programmers how to use the optimization routines in SAS IML to perform MLE
- Demonstrate best practices for MLE
- Compare the older NLPNRA and NLPQN optimization routines with the newer NLPSOLVE routine. The NLPSOLVE routine is only available in SAS Viya.
- Provide visualization routines for MLE
- Provide standard errors for parameter estimates

## Documentation

The SAS IML functions are described in the documentation for the MLE package. The file MLE_Doc.docx is a Word file that describes the syntax of each public function. The documentation shows how to call the functions and provides examples of output.

## Main functions

The following high-level functions are designed to be called directly. These functions all start with the 'MLE' prefix (all caps). They are documented and fully tested.
Although you can call the lower-level helper functions as well, the lower-level functions are documented only in the source code.

- **MLE**: The main function for fitting distributions to data. The function returns a (n_p x 1) vector of parameter estimates, where n_p
is the number of parameters in the distribution. You can get the names of the parameters by calling the 
**lik_dist_parmnames** function.
- **MLE_Fit**: Similar to the MLE function, but the MLE_Fit function returns a list of results, which is called the "MLE object." This list can be printed by passing it to the MLE_Summary function.
It can be used to create a graph by passing it to the MLE_Plot function. 
- **MLE_Summary**: Displays tables that summarize the results of the MLE_Fit function.
- **MLE_Plot**: Creates a histogram of the data and overlays the MLE model. You pass in the MLE object that is returned by the **MLE_Fit** routine.
- **MLE_MoM**: Returns the parameter estimates by using the method of moments (MoM). The MoM estimate can be used as an initial guess for MLE optimization.
- **MLE_LL**: Evaluate the loglikelihood function for a distribution at a specified value of the parameters.
- **MLE_Init** and **MLE_End**: Initializes and deletes (respectively) global variables needed for low-level routines.

