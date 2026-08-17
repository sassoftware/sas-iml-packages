resetline;proc iml;

/*---------------------------------------------------------------
 |  nmf: Nonnegative Matrix Factorization (SAS/IML)
 |  Inputs:
 |     a        -- input matrix
 |     k        -- approximation rank
 |     method   -- 'ALS': Alternating Least Square   or
 |                 'MUP': Multiplicative updates or
 |                 'APG': Alternating Proximal Gradient
 |     Opt      -- 1x5 numerical vector of options default:{0 200 10 1e-5 1e-5}
 |       Opt[1] = verbose  -- 0=none, 1=final, 2=all iterations
 |       Opt[2] = maxiter  -- maximum iterations
 |       Opt[3] = maxtry   -- Maximum number of initial guesses to be tried
 |       Opt[4] = normtol  -- convergence tolerance for objective function
 |       Opt[5] = xtol     -- convergence tolerance for factor change
 |     w0      -- initial factors:      "SVD", "RAND", or numerical  mxk matrix
 |     h0      -- initial coefficients: "SVD", "RAND", or numerical  kxn matrix
 |
 |  Outputs:
 |     w, h     -- factor matrices
 ---------------------------------------------------------------*/


/* Standardize W and H matrices, sort by column norms, and check for zero rows/columns */
start stdizeWH(w,h,status,tol=1e-4, input_zero_cols=);
   /* standardize the rows of H to have ssq of 1., scale W accordingly */
   hscsq=sqrt(h[,##]);
   low_rank=any(hscsq=0);            /* all zero row */
   hscsq = choose(hscsq=0, 1,hscsq); /* avoid divide-by-zero for zero-norm rows */
   w=w#hscsq`;
   h=h/hscsq;
   
   /* sort W columns based on their norm in descending order  */
   call sortndx(ndx, w[##,]`, ,1);   
   w=w[,ndx];
   h=h[ndx,];

   /* report if zero row or column is seen and optionally ignore zero_cols from input*/
   h_zero_cols=h[##,]<tol;
   if ^isskipped(input_zero_cols) then do;
      non_input_zero_cols = setdif(loc(h_zero_cols),input_zero_cols);
      status = low_rank>0 | ^isempty(non_input_zero_cols);
   end;
   else
      status = low_rank>0 | any(h_zero_cols);
finish;

/* Initialize W and H using the SVD of A. 
   This is a common initialization strategy for NMF. 
*/
start init_wh_svd(w0,h0,A,k);
   eps = 1E-2;
   call svd(U, Q, V, A);
   q=sqrt(Q[1:k]);
   w0=(U[,1:k] # q`) <>  eps; 
   h0=(q # (V[,1:k])`) <> eps;
finish;


/* This routine performs a single NMF pass given a w0 and h0. 
   The NMF routine by default tries up to 10 repetitions of try_nmf calls. 
*/
start try_nmf(w, h, rmsr, status, tryi, a, k, method, maxiter, normtol, xtol, verbose, w0, h0, zero_cols);
   status=0;
   m=nrow(a); n=ncol(a);
   mn = m*n;
   eps=constant("MACEPS");
   sqrteps = sqrt(eps);

   if (verbose>0 & tryi=1) then create __nmfdata__ var {"try", "j", "rms", "delta"};

   do j = 1 to maxiter;

      if method='ALS' then do;
         /*--- ALS: Alternating Least Squares ---*/
         CALL APPCORT (h, lindep, w0, a);       /* least-squares solve for H given W0 in W0*H ~ A */
         h = choose(h>0, h, 0);

         CALL APPCORT (w, lindep, h`, a`);w=w`; /* least-squares solve for W given H in W*H ~ A */
         w = choose(w>0, w, 0);
      end;
      else if method='APG' then do;
         /*--- APG: Proximal Gradient Descent ---*/
         gradH = w0` * (w0*h0 - a);
         Lh = norm(w0`*w0, 2);     /* Lipschitz constant */
         h = h0 - (1/Lh) * gradH;  /* gradient step */
         h = choose(h>0, h, 0);    /* projection */
         
         /* ---- Update W using APG ---- */
         gradW = (w0*h - a) * h`;
         Lw = norm(h*h`, 2);       /* Lipschitz constant */
         w = w0 - (1/Lw) * gradW;  /* gradient step */
         w = choose(w>0, w, 0);    /* projection */
      end;
      else do;
         /*--- Multiplicative Updates ---*/
         numer = w0` * a;
         denom = (w0`*w0)*h0 + eps;
         h = h0 # (numer / denom);
         h = choose(h>0, h, 0);

         numer = a * h`;
         denom = w0*(h*h`) + eps;
         w = w0 # (numer / denom);
         w = choose(w>0, w, 0);
      end;

      /*--- Compute norms and convergence ---*/
      d = a - w*h;
      rmsr = sqrt( sum(d#d) / mn );

      dw = max(abs(w-w0)) / (sqrteps + max(abs(w0)));
      dh = max(abs(h-h0)) / (sqrteps + max(abs(h0)));
      delta = max(dw, dh);
      
      if j>1 then do;
         if delta <= xtol then leave;
         delta_rmsr = rmsr0 - rmsr;
         if delta_rmsr >= 0 & delta_rmsr <= normtol*max(1, rmsr0) then leave;
      end;

      if verbose=2 then do;
         dbg_data= tryi || j || rmsr || delta;
         append from dbg_data;
      end;

      rmsr0 = rmsr;
      w0 = w;
      h0 = h;
   end;

   /* standardize the rows of H and columns of W */
   run stdizeWH(w,h,status, ,zero_cols);

   converged = (j<=maxiter);                /* LEAVE keeps j; normal completion gives maxiter+1 */
   if verbose>0 then do;
      iter = min(j, maxiter);
      if (verbose=1 | converged) then do;   /* don't duplicate the final MaxIter row */
         dbg_data = tryi || iter || rmsr || delta;
         append from dbg_data;
      end;
   end;
finish;   


start nmf(w, h, input, k, method='ALS', opt={0 . 10 1e-5 1e-5}, w0='SVD', h0='SVD');

   /* Validate argument types before doing any numeric comparison */
   if (type(opt)^='N' | type(k)^='N' | type(input)^='N') then do;
      call PrintToLog('The input matrix, the rank k, and Opt must be numeric.', 2);
      stop;
   end;

   if any(input=.) then do;
      msg = 'Input matrix contains missing values. Please impute or remove missing values and rerun.';
      call PrintToLog(msg, 2);
      stop;
   end;
   A = input;

   /* Normalize and validate the method up front so the defaults below match */
   if (type(method)^='C') then do;
      call PrintToLog('The Method argument must be a character matrix.', 2);
      stop;
   end;
   method = strip(upcase(method));
   if (method^='ALS' & method^='MUP' & method^='APG') then do;
      call PrintToLog('Invalid method. Use ALS, MUP, or APG.', 2);
      stop;
   end;

   /* Check and verify the options, set defaults if not specified, and assign to named variables for readability */
   options=rowvec(opt);
   if ncol(options)<5 then 
     options = options || j(1, 5-ncol(options), .);

   if options[1]=. then verbose=0;    else verbose=options[1];
   if options[3]=. then maxtry =10;   else maxtry =options[3]; /* for ALS method, maxtry=1 */
   if options[4]=. then normtol=1e-5; else normtol=options[4];
   if options[5]=. then xtol   =1e-5; else xtol   =options[5];
   if options[2]=. then do;
      if      method='ALS' then maxiter=200;
      else if method='MUP' then maxiter=20000;
      else if method='APG' then maxiter=20000;
   end;
   else maxiter=options[2];


   SVD_INIT=0;
   if type(w0)='C' then do;
      w0=strip(upcase(w0));
      if w0='SVD' then SVD_INIT=1;
   end;
   if type(h0)='C' then do;
      h0=strip(upcase(h0));
      if h0='SVD' then SVD_INIT=1;
   end;
   m=nrow(a); n=ncol(a);

   if (type(w0)^='N' & type(w0)^='C') | 
      (type(w0)='N' & nrow(w0)^=m ) then do;
      call PrintToLog('Invalid W0 specified.', 2); 
      stop; 
   end;

   if (type(w0)='C') then do;
      if (w0 ^='SVD' & w0 ^='RAND') then do;
         call PrintToLog('Invalid W0 specified.', 2); 
         stop; 
      end;
   end;

   if (type(h0)^='N' & type(h0)^='C') | 
      (type(h0)='N' & ncol(h0)^=n ) then do;
      call PrintToLog('Invalid H0 specified.', 2); 
      stop; 
   end;

   if (type(h0)='C') then do;
      if (h0^='SVD' & h0^='RAND') then do;
         call PrintToLog('Invalid H0 specified.', 2); 
         stop; 
      end;
   end;
   
   w0_len=ncol(w0)*nrow(w0);
   h0_len=ncol(h0)*nrow(h0);
   if (w0_len>1 & h0_len=1 & type(w0)=type(h0)) | ((w0_len=1 & h0_len>1 & type(w0)=type(h0))) then do;
      call PrintToLog('Incompatible H0 and W0 identified.', 2); 
      stop; 
   end;

   if (ncol(a)=0) then do; call PrintToLog('Input matrix cannot be empty.', 2); stop; end;
   if (any(a<0)) then do;  call PrintToLog('Input matrix must be non-negative.', 2); stop; end;
   if (k<1 | k>n | k^=round(k)) then do; call PrintToLog('Invalid rank.', 2); stop; end;
   if (maxiter<1 | normtol<0 | xtol<0 |verbose<0 | maxtry<1) then do; 
                           call PrintToLog('Invalid input.', 2); stop; end;
   
   if type(w0)='N'  then do;
      if any(w0<0) then do; call PrintToLog('W0 must be non-negative.', 2); stop; end;
   end;
   
   if type(h0)='N' then do;
      if any(h0<0) then do; call PrintToLog('H0 must be non-negative.', 2); stop; end;
   end;

   w0t=w0; h0t=h0; /* save the user-specified w0 and h0 */
   if (SVD_INIT=1 & ncol(A)<=5000) then do;
      run init_wh_svd(w0t,h0t,a,k);
      if options[3]=. then maxtry =1; /* if SVD is used, only one try is needed */
    end;
   else do;
      /* For 'RAND', or for 'SVD' on a large matrix, fall back to random factors */
      if (type(w0)='C') then w0t=randfun(m||k,"Uniform");
      if (type(h0)='C') then h0t=randfun(k||n,"Uniform");
   end;

   tmp=a[<>,];
   zero_cols=loc(tmp[##,]<1e-6);
   free tmp;

   status =1; 
   rmsr=1e16; err=1;
   
   do tryi=1 to maxtry while(status>0 | err>1e-4);
      
      run try_nmf(w_i, h_i, rmsr_i, status, tryi, a, k, method, maxiter, normtol, xtol, verbose, w0t, h0t, zero_cols);
      err = abs(rmsr-rmsr_i);
      
      if (tryi=1 | rmsr_i<rmsr) then do;
         w=w_i; 
         h=h_i; 
         rmsr=rmsr_i;
      end;
      
      w0t=randfun(m||k,"Uniform");
      h0t=randfun(k||n,"Uniform");
   end;

   if verbose>0 then do;
      close __nmfdata__;
      submit verbose;
         title "CALL NMF Iteration Steps";
         proc print data=work.__nmfdata__ label noobs;
            label try="Trial Number"
                  j  ="Iteration Number"
                 rms ="RMS Residual"
               delta ="Max Relative Change";
         run;
         
         title;
      endsubmit;
      call delete("work", "__nmfdata__");
   end;

finish;


/* 
Helper IML subroutine, named proc_nmf, which calls PROC NMF to 
compute the NMF factorization.
Because PROC NMF is only available in SAS Viya 4.0 and later, this 
function is conditionally defined based on the SYSVER macro.
PROC NMF on the input matrix A with rank k, and return W and H. 
If if_stdize is nonzero, also standardize the output W and H. 
The inputs are:
   A: the input matrix to factorize
   k: the target rank for the factorization
   if_stdize: whether to standardize W and H (default 0)
   seed: the random seed for PROC NMF (default 12345)
   options: additional options to pass to PROC NMF (default '')
   timing: if provided, will be set to the time taken by PROC NMF
The outputs are:
   w_proc_nmf: the W matrix from PROC NMF
   h_proc_nmf: the H matrix from PROC NMF
*/
start proc_nmf(w_proc_nmf, h_proc_nmf, A, k, if_stdize=0, seed=12345, options='', timing=);
   create mydata from A;
   append from A;
   close mydata;

   submit k seed options;
      ods exclude all; 
      ods output Timing=NMFTimeData;
      proc nmf data=mydata nthreads=64 rank=&k seed=&seed outh=H &options;
         var col:;
         output out=W;
         *display / excludeall;
      run;
      ods output close;
      ods exclude none;
   endsubmit;
   use NMFTimeData; read all var {Time} into proc_time_matrix; close NMFTimeData;
   if ^isskipped(timing) then do;
      timing=proc_time_matrix[nrow(proc_time_matrix)];
   end;

   use work.H(keep=col:); read all into h_proc_nmf; close work.H;
   h_proc_nmf=h_proc_nmf[,2:ncol(h_proc_nmf)]; /* remove the first column which is the ID variable */ 
   use work.W; read all into w_proc_nmf; close work.W;

   status = 0;   /* default when the convergence check is skipped */
   if (if_stdize^=0) then run stdizeWH(w_proc_nmf, h_proc_nmf, status);
   if (status) then call PrintToLog('PROC NMF did not converge.', 1);
finish;

store module=(stdizeWH init_wh_svd try_nmf nmf proc_nmf);
QUIT;
