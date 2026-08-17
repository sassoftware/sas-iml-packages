/*---------------------------------------------------------------
 |  imagenmf.sas: Image approximation with Nonnegative Matrix
 |                Factorization (NMF)
 |
 |  Purpose:
 |     Illustrate how increasing the approximation rank k improves
 |     the NMF reconstruction of a simple binary image.
 |
 |  Description:
 |     1. Build a 7x19 binary matrix A whose nonzero entries spell
 |        the letters "NMF".
 |     2. Display A as a grayscale heatmap ("Original Image").
 |     3. For each rank k = 1, ..., 5, factorize A as A ~ W*H using
 |        CALL NMF with the Alternating Least Squares ('ALS') method.
 |     4. Reconstruct X = W*H, rescale it so that max(X) = 1, and
 |        display X as a heatmap titled "Rank = k".
 |     As k increases, the reconstruction sharpens until the "NMF"
 |        text is clearly legible.
 |
 |  Requires:
 |     The NMF modules stored by nmf.sas (run nmf.sas first so that
 |     LOAD MODULE=_ALL_ can find them). Available in SAS 9.4 and later.
 ---------------------------------------------------------------*/
proc iml;
load module=_all_;

NMF = {0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0,
       0 1 0 0 0 1 0 1 0 0 0 1 0 1 1 1 1 1 0,
       0 1 1 0 0 1 0 1 1 0 1 1 0 1 0 0 0 0 0,
       0 1 0 1 0 1 0 1 0 1 0 1 0 1 1 1 1 0 0,
       0 1 0 0 1 1 0 1 0 0 0 1 0 1 0 0 0 0 0,
       0 1 0 0 0 1 0 1 0 0 0 1 0 1 0 0 0 0 0,
       0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 };
A = NMF;
ramp = {CXFFFFFF CXF0F0F0 CXBDBDBD CX636363 CX000000};
call heatmapcont(A) colorramp=ramp showlegend=0 title="Original Image" range={0,1};
m=nrow(A); n=ncol(A);
mn = m*n;

call randseed(1234,1);
/* opt = verbose || maxiter || maxtry || normtol || xtol; */
maxtry=10; verbose=0; xtol=1e-4; normtol=1e-4;
opt = verbose || . || maxtry || normtol || xtol;
do k=1 to 5;
    call nmf(w_mult,h_mult,A,k) method='ALS' Opt=opt;
    s = "Rank = " + char(k,1);
    X = w_mult*h_mult;
    X = X / max(X); /* scale X so that max(X)=1 */
    call heatmapcont(X) colorramp=ramp showlegend=0 title=s range={0,1};
    *print X[F=4.2];
end;
quit;