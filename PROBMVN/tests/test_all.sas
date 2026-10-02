%let rc = %sysfunc(dlgcdir('U:\gitpp\DEV\sas-iml-packages\PROBMVN\tests'));

%include "test_cdfbvn.sas";
%include "test_cdftvn.sas";
%include "test_cdfmvn.sas";
%include "test_probmvn_Doc.sas";
%include "test_probmvn_opt.sas";
%include "test_probmvn_2D.sas";
%include "test_probmvn_big.sas";
%include "test_probmvn.sas";
