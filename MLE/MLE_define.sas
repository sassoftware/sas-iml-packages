/* Master %INCLUDE for all module definitions. */
%let MLE_path = u:\gitpp\DEV\sas-iml-packages\MLE;

proc iml;
%include "&MLE_path/MLE_Util.sas";
%include "&MLE_path/MLE_Fit.sas";
%include "&MLE_path/MLE_Keywords.sas";
%include "&MLE_path/MLE_LL.sas";
%include "&MLE_path/MLE_MoM.sas";
%include "&MLE_path/MLE_Plot.sas";
%include "&MLE_path/MLE_Summary.sas";
quit;

