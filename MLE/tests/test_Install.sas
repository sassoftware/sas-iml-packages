/* Use Git to fetch the package from GitHub into a temporary directory such as WORK.
   To learn more about using git commands in SAS to read packages, See  
   https://blogs.sas.com/content/iml/2023/02/06/git-share-sas-programs.html
*/
options dlcreatedir;
%let repoPath = %sysfunc(getoption(WORK))/sas-iml-packages;  
%let MLE_path = &repoPath/MLE; 

data _null_;
if fileexist("&repoPath.") then 
   put 'Repository already exists; skipping the clone operation'; 
else do;
   put "Cloning repository 'sas-iml-packages'";
   rc = gitfn_clone("https://github.com/sassoftware/sas-iml-packages/", "&repoPath." ); 
end;
run;

/* Use %INCLUDE to read source code and STORE functions to IML storage library */
proc iml;
%include "&MLE_path/MLE_Define.sas";  /* one file with all modules */
quit;