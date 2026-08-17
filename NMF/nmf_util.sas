/* nmf_util.sas: Utility routines for the NMF package */


/*
   Check the SYSVER macro to see if SAS 9.4 is running.
   In SAS Viya, the macro is empty and does nothing.
   In SAS 9.4, the macro defines a function that emulates the PrintToLog call.
   The syntax is as follows:
      call PrintToLog("This is a log message.");
      call PrintToLog("This is a note.", 0);
      call PrintToLog("This is a warning.", 1);
      call PrintToLog("This is an error.", 2);
*/
%macro DefinePrintToLog;
   %if %sysevalf(&sysver = 9.4) %then %do;
      start PrintToLog(msg,errCode=-1);
         if      errCode=0 then prefix = "NOTE: ";
         else if errCode=1 then prefix = "WARNING: ";
         else if errCode=2 then prefix = "ERROR: ";
         else prefix = "";
         stmt = '%put ' + prefix + msg + ';';
         call execute(stmt);
      finish;
      store module=(PrintToLog);

      start ErrorToLog(msg);
         run PrintToLog(msg, 2);
      finish;
      store module=(ErrorToLog);
   %end;
%mend;

%macro DefineProcNMF;
   %if %substr(&sysver, 1, 4) ^= V.04 %then %do;
      start proc_nmf(w_proc_nmf, h_proc_nmf, A, k, if_stdize=0, seed=12345, options='', timing=);
         msg = 'PROC NMF is not supported in this version of SAS. Please use SAS VIYA 4.0 or later.';
         call PrintToLog(msg, 2);
         stop;
      finish;
      
      /*Overwrite proc_nmf function in SAS V9.4 since proc nmf is shipped with Viya*/
      store module=(proc_nmf); 
   %end;
%mend;

proc iml;
%DefinePrintToLog;
%DefineProcNMF;
QUIT;