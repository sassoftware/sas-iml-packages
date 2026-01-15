# SAS® IML Packages

## Overview

A *package* is a collection of SAS IML function modules that share a common purpose. Programmers can store the modules and use them in their SAS IML programs. 

Packages in SAS Viya® are managed by using a repository on GitHub. You can clone the repository to obtain a copy of the functions. When you execute the source files, the functions are stored in a module library. You can then load the functions and use them in your programs.

This method enables you to use libraries of function modules in the IML procedure and in the IML action. The same process also works for PROC IML in SAS® 9.4.

## List of Packages

* The **Compgeom** package is supported only in SAS Viya. The package shows how to use the CONVEXHULL, DELAUNAY, and VORONOI subroutines (released in Viya 2024.10) to solve problems in computational geometry. The package contains visualization subroutines that you can run in PROC IML. The routines visualize two-dimensional convex hulls, Delaunay triangulations, and Voronoi diagrams.
* The **Metalog** package is supported for SAS 9.4 and for SAS Viya. The package supports using the metalog system of distributions.
* The **MLE** package is supported for SAS 9.4 and for SAS Viya. The package demonstrates how to fit and visualize common univariate probability distributions.

### Prerequisites

* The graphical subroutines must be run from PROC IML.
* The Compgeom package requires SAS Viya 2024.10.
* The Metalog package works in every release of SAS Viya. It also works in SAS/IML® in SAS 9.4 Maintenance 7 and beyond.
* The MLE package works in every release of SAS Viya. It also works in SAS/IML® in SAS 9.4 Maintenance 7 and beyond.

## Getting Started

To use the packages, you must download the source code and include the source files in your SAS IML program.

### Installing the Packages

When you download the repository, you get all packages. Each package is contained in its own directory. Each package contains documentation that shows how to use the package. The installation steps depend on the client: SAS, R, Python, and so forth. 

On the SAS client, you can use the GITFN_CLONE function in the DATA step to clone the repository to a local directory or to a libref such as WORK. You can then use a %INCLUDE statement to run the statements that define and store the function modules. For example, in PROC IML, you can use the following statements to clone the sas-iml-packages repository to the WORK libref:

```sas
%let gitURL = https://github.com/sassoftware/sas-iml-packages/;
options dlcreatedir;
%let repoPath = %sysfunc(getoption(WORK))/sas-iml-packages/;

/* clone repository; if repository exists, skip download */
data _null_;
if fileexist("&repoPath.") then 
   put 'Repository already exists; skipping the clone operation'; 
else do;
   put 'Cloning repository sas-iml-packages';
   rc = gitfn_clone("&gitURL", "&repoPath." ); 
end;
run;
```

If you want the package files to remain after the SAS session ends, you can download the files to a permanent location. 

### Running a Package

You need to clone the repository only one time. Whenever you want to use a package, use the %INCLUDE statement on the SAS client to read and define the functions in the package. For example, if you want to use the **Metalog** package, use the following statements:

```sas
proc iml;
%include "&repoPath./Metalog/ML_define.sas";
%include "&repoPath./Metalog/ML_proc.sas";
/* use the functions here */
quit;
```
If you want to use the **Compgeom** package, use the following statements:

```sas
proc iml;
%include "&repoPath./Compgeom/CG_define.sas";
%include "&repoPath./Compgeom/CG_proc.sas";
/* use the functions here */
quit;
```

If you want to use the **MLE** package, use the following statements:

```sas
proc iml;
%include "&repoPath./MLE/MLE_define.sas";
quit;
```

Each file ends with a STORE statement, which means that the functions are stored to the active library. Consequently, you can use the LOAD statement to read the functions into subsequent calls to PROC IML in the same SAS session:

```sas
proc iml;
LOAD MODULE= _ALL_;
/* use the functions here */
quit;
```

See the individual packages for documentation and examples of usage.

## Contributing

We welcome your contributions. Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details about submitting contributions to this project. 

## License

This project is licensed under the [Apache 2.0 License](LICENSE).

## Additional Resources

The directory for each package contains references and additional resources for the package.
