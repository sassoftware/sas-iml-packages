/* Master %INCLUDE for all module definitions. */
/**** You MUST define the PKG_path macro before including this file. ****/

/*For example:
  %let repoPath = u:\gitpp\DEV\sas-iml-packages;
  proc iml;
  %include "&repoPath\PROBBVN\PROBMVN_define.sas";
  quit;

  For Windows vs Linux, you might need to use backslashes instead of forward slashes.
*/

%include "&repoPath\PROBMVN\probmvn_Util.sas";
%include "&repoPath\PROBMVN\probbvn.sas";
%include "&repoPath\PROBMVN\cdftvn.sas";
%include "&repoPath\PROBMVN\cdfmvn.sas";
%include "&repoPath\PROBMVN\probmvn.sas";

/* Alternatively, you can use the DLGCDIR to set the current working directory,
   and then refer to the files without specifying a path. Here is a Linux example:

%let rc = %sysfunc(dlgcdir('U:\gitpp\DEV\sas-iml-packages\PROBMVN'));
%include "probmvn_Util.sas";
%include "probbvn.sas";
%include "cdftvn.sas";
%include "cdfmvn.sas";
%include "probmvn.sas";

*/
