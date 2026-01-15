/* Master %INCLUDE for all module definitions. */
/**** You MUST define the MLE_path macro before including this file. ****/

/*For example:
  %let MLE_Path = u:\gitpp\DEV\sas-iml-packages\MLE;
  proc iml;
  %include "&MLE_path/MLE_Define.sas";
  quit;

  For Windows vs Linux, you might need to use backslashes instead of forward slashes.
*/

%include "&MLE_path/MLE_Util.sas";
%include "&MLE_path/MLE_Fit.sas";
%include "&MLE_path/MLE_Keywords.sas";
%include "&MLE_path/MLE_LL.sas";
%include "&MLE_path/MLE_MoM.sas";
%include "&MLE_path/MLE_Plot.sas";
%include "&MLE_path/MLE_Summary.sas";
