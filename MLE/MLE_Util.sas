/*********************************************************/
/* MLE_Util.sas                                          */
/* Define initialization functions and helper functions. */
/*********************************************************/

/* define functions for SAS 9.4 that are built-in for Viya */

/* Check the SYSVER macro to see if SAS 9.4 is running.
   In SAS Viya, the macro is empty and does nothing.
   In SAS 9.4, the macro defines a function that emulates the 
   COMPLETECASES function in SAS Viya.
*/
%macro DefineSAS9UtilFuncs;
%IF %sysevalf(&sysver = 9.4) %THEN %DO;

/* Return the rows of a matrix M that are nonmissing. 
   Default: return a vector of the row numbers (indices) that are complete.
   If method="EXTRACT", return the data values for the complete rows.

   EXAMPLE:
   M = {1 2, 3 ., . 5, . ., 9 10};
   nonMissIdx = CompleteCases(M);
   nonMissData = CompleteCases(M, "Extract");
*/
start CompleteCases(M, method="ROW");
   idx = loc( countn(M, "row")=ncol(M) );
   if ncol(idx)=0 then return(idx);
   if upcase(method)="ROW" then 
      return( idx );
   return( M[idx, ] );
finish;

/* Check the SYSVER macro to see if SAS 9.4 is running.
   In SAS Viya, the macro is empty and does nothing.
   In SAS 9.4, the macro defines a function that emulates the PrintToLog call.
   The syntax is as follows:
   call PrintToLog("This is a log message");
   call PrintToLog("This is a note", 0);
   call PrintToLog("This is a warning", 1);
   call PrintToLog("This is an error", 2);
*/
start PrintToLog(msg,errCode=-1);
   if      errCode=0 then prefix = "NOTE: ";
   else if errCode=1 then prefix = "WARNING: ";
   else if errCode=2 then prefix = "ERROR: ";
   else prefix = "";
   stmt = '%put ' + prefix + msg + ';';
   call execute(stmt);
finish;
store module=(
      CompleteCases
      PrintToLog
      );
%end;
%mend DefineSAS9UtilFuncs;

%DefineSAS9UtilFuncs;

/* Indirect calling of a function that has 1 arg. In SAS 9
   use CALL EXECUTE. In Viya, use FUNCEVAL.
   EXAMPLE:
    start foo(x); 
       return(x+1);
    finish;

    %EVALFUNC1(result1, foo, 7);
    print result1;
    funcname = "foo";
    %EVALFUNC1(result2, funcname, 7);
    print result2;
*/
%macro EVALFUNC1(rOut, fname, arg);
%if %sysevalf(&sysver = 9.4) %then %do;
   if type(&fname.)='C' then 
      cmd = "&rOut. = " + value("&fname.") + "( &arg. );";
   else 
      cmd = "&rOut. = &fname.( &arg. );";
   call execute( cmd );
%end;
%else %do;
   if type(&fname.)='C' then 
      &rOut. = funceval( &fname., &arg. );
   else if type(&fname.)='U' then
      &rOut. = funceval( "&fname.", &arg. );
   else 
      call PrintToLog('Incorrect call to %EVALFUNC1 macro', 1);
%end;
%mend EVALFUNC1; 

/* Syntax:
     %EmulateHistogramSetup(dsIn=DATASET, varIn=VARIABLE)
   where 
     DATASET = name of a SAS data set 
     VARIABLE= name of variable in data set whose distribution you want to model

   The macro does the following:
      1. Writes a data set called _HistBins that contains variables 
         _MIDPT_ : Centers of histogram bins
         _COUNT_ : Frequency count in each bin
         _OBSPCT_: Percentage of observations in each bin
         _ZERO_  : The constant value 0, which is the lower edge of the high-low plot
      2. Creates the following macro variables:
         &_BINSTART : the value of the center of the first bin
         &_BINEND   : the value of the center of the last bin
         &_BINWIDTH : the width of the bins
         &_VARNAME  : the name of the variable whose distribution is modeled
   You can emulate a histogram by using the HIGHLOW stmt in PROC SGPLOT:
   proc sgplot data=_HistBins;
      highlow x=_midpt_ low=_zero_ high=_obspct_ / type=bar barwidth=1;
      yaxis min=0 offsetmin=0 grid;
      xaxis values=(&_binStart to &_binEnd by &_binWidth) valueshint;
   run;
*/
%macro EmulateHistogramSetup(dsIn=, varIn=);
%global _binStart _binEnd _binWidth _NObs;
proc univariate data=&dsIn noprint;
   var &varIn;
   histogram &varIn / outhist=_HistBins(rename=(_OBSPCT_=_PCT_)) noplot;
   output out=_HistOut n=_NOBS_;
run;
data _null_;
   set _HistBins end=EOF;
   if _N_=1 then 
      call symputx("_binStart", _MIDPT_);
   h = dif(_MIDPT_);
   if EOF then do;
      call symputx("_binEnd", _MIDPT_);
      call symputx("_binWidth", h);
      call symputx("_varName", &varIn);
   end;
run;
data _null_;
   set _HistOut;
   call symputx("_NOBS", _NOBS_);
run;

data _HistBins;
set _HistBins;
_ZERO_ = 0;        /* add baseline for histogram */
label _MIDPT_ = &varIn
      _PCT_   = "Percent"
      _COUNT_ = "Count";
run;
/*
%PUT The following macro variables have been defined:;
%PUT &=_binStart &=_binEnd &=_binWidth &=_NObs;
*/
%mend EmulateHistogramSetup;
