/*---------------------------------------------------------------
 |  exnmf.sas: Example of Nonnegative Matrix Factorization (NMF)
 |
 |  Purpose:
 |     Demonstrate the CALL NMF subroutine (SAS/IML implementation)
 |     and validate its result against PROC NMF.
 |
 |  Description:
 |     1. Build a small 5x5 nonnegative test matrix A.
 |     2. Factorize A as A ~ W*H of rank k using CALL NMF with the
 |        Alternating Least Squares ('ALS') method.
 |     3. Compute the root-mean-square residual (dn) of the fit.
 |     4. Repeat the factorization with PROC NMF (via proc_nmf) and
 |        compute its residual (dn_proc_nmf).
 |     5. The test passes if CALL NMF achieves a residual no larger
 |        than PROC NMF; the resulting W and H matrices are printed
 |        side by side for comparison.
 |
 |  Requires:
 |     The NMF modules stored by nmf.sas (run nmf.sas first so that
 |     LOAD MODULE=_ALL_ can find them). PROC NMF requires SAS Viya.
 ---------------------------------------------------------------*/
proc iml;
load module=_all_;

A = diag({4,2,3,6,9});  /* 5x5 example */
k = 3; 

verbose=0; maxiter=200; xtol=1e-6; normtol=1e-6;
opt = verbose || maxiter || . || normtol || xtol;
call nmf(w,h,a,k) method='ALS' Opt=opt; /* opt = verbose || maxiter || maxtry || normtol || xtol; */

d =(a-w*h);
dn= sqrt( sum(d#d) / ncol(a)/nrow(a) );

run proc_nmf(w_proc_nmf, h_proc_nmf, A, k, 1);
d =(a - w_proc_nmf * h_proc_nmf);
dn_proc_nmf= sqrt( sum(d#d) / ncol(a)/nrow(a) );

if dn<=dn_proc_nmf then print 'NMF test passed!';
else print dn, dn_proc_nmf ;
print (round(w,1e-4))[l='w'] (round(w_proc_nmf,1e-4))[l='w_proc_nmf'],
      (round(h,1e-4))[l='h'] (round(h_proc_nmf,1e-4))[l='h_proc_nmf'];
quit;
