proc iml;
%include "&MLE_path./MLE_util.sas";
QUIT;

/* test calls to the top-level subroutine */
PROC IML;
%IF (&sysver = 9.4) %THEN %DO;
load module=(
     CompleteCases
     PrintToLog
     );
%END;

/* CompleteCases */
M = {1 2, 3 ., . 5, . ., 9 10};
nonMissIdx = CompleteCases(M);
nonMissData = CompleteCases(M, "Extract");
print nonMissIdx nonMissData;

/* PrintToLog */
call PrintToLog("This is a log message");
call PrintToLog("This is a note", 0);
call PrintToLog("This is a warning", 1);
call PrintToLog("This is an error", 2);

/* %EVALFUNC1 */
start foo(x); 
   return(x+1);
finish;

%EVALFUNC1(result1, foo, 7);
print result1;
funcname = "foo";
%EVALFUNC1(result2, funcname, 7);
print result2;

QUIT;
