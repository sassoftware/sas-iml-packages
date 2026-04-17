/* Master %INCLUDE for all module definitions. */
/**** You MUST define the PKG_path macro before including this file. ****/

/*For example:
  %let PKG_Path = u:\gitpp\DEV\sas-iml-packages;
  proc iml;
  %include "&PKG_path\PROBBVN\probmvn_Define.sas";
  quit;

  For Windows vs Linux, you might need to use backslashes instead of forward slashes.
*/

%include "&PKG_path\PROBMVN\cdfmvn_Util.sas";
%include "&PKG_path\PROBMVN\cdfbnv.sas";
%include "&PKG_path\PROBMVN\cdftvn.sas";
%include "&PKG_path\PROBMVN\cdfmvn.sas";
