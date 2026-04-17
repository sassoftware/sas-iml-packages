/* Master %INCLUDE for all module definitions. */
/**** You MUST define the PKG_path macro before including this file. ****/

/*For example:
  %let repoPath = u:\gitpp\DEV\sas-iml-packages;
  proc iml;
  %include "&repoPath\PROBBVN\probmvn_Define.sas";
  quit;

  For Windows vs Linux, you might need to use backslashes instead of forward slashes.
*/

%include "&repoPath\PROBMVN\cdfmvn_Util.sas";
%include "&repoPath\PROBMVN\cdfbvn.sas";
%include "&repoPath\PROBMVN\cdftvn.sas";
%include "&repoPath\PROBMVN\cdfmvn.sas";
