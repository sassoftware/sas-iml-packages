# The MLE repo

# THIS PACKAGE IS UNDER DEVELOPMENT. DO NOT USE.

## Description

This project is a set of helper routines for fitting univariate distributions to data by using maximum likelihood estimation (MLE).
The functions are written in the SAS IML language. The functions serve several purposes:
- Show programmers how to use the optimization routines in SAS IML to perform MLE
- Demonstrate best practices for MLE
- Compare the older NLPNRA and NLPQN optimization routines with the newer NLPSOLVE routine. The NLPSOLVE routine is onlyu available in SAS Viya.
- Provide visualization routines for MLE
- Provide standard errors for parameter estimates

The SAS IML functions are described in the documentation for the MLE package.

## Main functions

The following high-level functions are designed to be called directly. These functions all start with the 'MLE' prefix. They are documented and fully tested.
Although you can call the lower-level helper functions as well, the lower-level functions might be less robust and less documented.

- **MLE**: The main function for fitting distributions to data. The function returns a (1 x n_p) vector of parameter estimates, where n_p
is the number of parameters in the distribution. 
- **MLE_Fit**: Similar to the MLE function, but the MLE_Fit function returns a list of results. This list can be printed by passing it to the MLE_Summary function.
It can be used to create a graph by passing it to the MLE_Plot function. 
- **MLE_Plot**: Creates a histogram of the data and overlays the MLE model.
- **MLE_Summary**: Displays tables that summarize the results of the MLE_Fit function.
- **MLE_MoM**: Returns the parameter estimates by using the method of moments (MoM). The MoM estimate can be used as an initial guess for MLE.
- **MLE_LL**: Evaluate the loglikelihood function for a distribution at the specified value of the parameters.
- **MLE_Init** and **MLE_End**: Initializes and deletes (respectively) global variables needed for low-level routines.

## Documentation

The file MLE_Doc.docx is a Word file that describes the syntax of each public function. The documentation shows how to call the functions and provides examples of output.