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

start IsSAS9(_=0);
%IF (&sysver = 9.4) %THEN %DO;
   return(1);
%END;
%ELSE %DO;
   return(0);
%END; 
finish;
store module=(
   IsSAS9
   );

%IF (&sysver = 9.4) %THEN %DO;

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
%if (&sysver = 9.4) %then %do;
   if type(&fname.)='C' then 
      _cmd = "&rOut. = " + value("&fname.") + "( &arg. );";
   else 
      _cmd = "&rOut. = &fname.( &arg. );";
   call execute( _cmd );
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
