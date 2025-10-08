/***********************************/
/* Define plotting functions for PROC IML */
/* See                                    */
/***********************************/

/* call PROC UNIVARIATE to output histogram bins.
   Return the k x 2 matrix, M, where 
      M[,1] = midpoints of bins
      M[,2] = percentages of obs in each bin */
start mle_Plot_HistBins(Y, varName=parentname("Y"));
    create _HistData from Y[c=varName];
    append from Y;
    close;
    submit varName;
      %let _currNotes = %sysfunc(getoption(NOTES)); /* save state */
      OPTIONS NONOTES;
      %EmulateHistogramSetup(dsIn=_HistData, varIn=&varName);
      OPTIONS &_currNotes;                          /* restore state */
    endsubmit;
    use _HistBins;
    read all var {'_MIDPT_' '_PCT_'} into M;
    close;
    return( M );
finish;

/* Input: midPts is the _MIDPTS_ column that contains the midpoints of histogram bins
   Output:
      h        = width of bins
      startBin = left edge of first bin
      endBin   = right edge of last bin 
*/
start mle_Plot_BinsAndWidth(h, startBin, endBin, midPts);
   nBins = nrow(midPts) * ncol(midPts);
   h = midPts[2] - midPts[1];
   startBin = midPts[1] - h/2;
   endBin = midPts[nBins] + h/2;
finish;

/* Given a variable X and a set of distributions with parameters,
    overlay the histogram of X with the fitted distributions.

   INPUT:
      X         : (n x 1) vector of data values
      DistNames : (k x 1) vector of distribution names (eg, 'Normal')
      params    : (k x 3) matrix of parameters, one row for each distribution
   SIDE EFFECTS:
      Writes the data set _HistBins that contains the histogram bin info
      Writes the data set _Model that contains the model fit info
   OUTPUT:
      A plot is produced that overlays a histogram of X with the fitted distributions.
*/
start mle_Plot_Overlay(X, DistNames, params, varName=parentname("Y"));
   M = mle_Plot_HistBins(X, varName);   /* side effect: writes _HistBins */
   call mle_Plot_BinsAndWidth(h, b1, bn, M[,1]);
   t = T( do( b1, bn, (bn-b1)/49) );
   call MLE_PDF(Labl, PDF, DistNames, params, t);
   /* scale the density curves to the percentage scale. 
      See https://blogs.sas.com/content/iml/2024/06/19/scale-density-curve-histogram.html */
   Pct = h*100*PDF;
   /* convert from wide to long format; use Labl as group values */
   call widetolong(LongX, LongY, Group, Pct, t, Labl);
   /* Write model Pct to second data set */
   create _Model from LongX LongY Group[c=(varName || {"Pred_Pct" "Distrib"})];
   append from LongX LongY Group;
   close;   

   /* create a label statement for variables in the _MODEL data set */
   labelStmt = 'label Pred_Pct="Percent" Distrib="Distribution";';
   /* overlay a high-low plot (emulate a histogram) and the model */
   submit varName LabelStmt;
   data _HistOverlay;
   set _HistBins _Model;
   &LabelStmt;
   run;
   proc sgplot data=_HistOverlay;
      highlow x=_midpt_ low=_zero_ high=_pct_ / type=bar barwidth=1;
      series x=&varName y=Pred_Pct / group=Distrib nomissinggroup;
      yaxis min=0 offsetmin=0 grid;
      xaxis values=(&_binStart to &_binEnd by &_binWidth) valueshint;
   run;
   endsubmit;
finish;

/* The list L contains the following fields:
   FitObj$"Dist" = a keyword for the distribution
   FitObj$"y" = the raw data
   FitObj$"ParmNames" = the parameter names
   FitObj$"Estimate" = the parameter estimates

   Extract the relevant information and call the lower-level routine to create the plot.
*/
start MLE_Plot(L, L2=, L3=, L4=, L5=, L6=, L7=, L8=, L9=, L10=);
    y = L$"y";
    Dist = L$"Dist";
    params = rowvec(L$"Estimate");
    if ^isSkipped(L2) then do;   Dist=Dist//L2$"Dist";  params=params//rowvec(L2$"Estimate");  end;
    if ^isSkipped(L3) then do;   Dist=Dist//L3$"Dist";  params=params//rowvec(L3$"Estimate");  end;
    if ^isSkipped(L4) then do;   Dist=Dist//L4$"Dist";  params=params//rowvec(L4$"Estimate");  end;
    if ^isSkipped(L5) then do;   Dist=Dist//L5$"Dist";  params=params//rowvec(L5$"Estimate");  end;
    if ^isSkipped(L6) then do;   Dist=Dist//L6$"Dist";  params=params//rowvec(L6$"Estimate");  end;
    if ^isSkipped(L7) then do;   Dist=Dist//L7$"Dist";  params=params//rowvec(L7$"Estimate");  end;
    if ^isSkipped(L8) then do;   Dist=Dist//L8$"Dist";  params=params//rowvec(L8$"Estimate");  end;
    if ^isSkipped(L9) then do;   Dist=Dist//L9$"Dist";  params=params//rowvec(L9$"Estimate");  end;
    if ^isSkipped(L10) then do;  Dist=Dist//L10$"Dist"; params=params//rowvec(L10$"Estimate"); end;
    call mle_Plot_Overlay(y, Dist, params);
    return;
finish;


store module=(
   mle_Plot_HistBins
   mle_Plot_BinsAndWidth
   mle_Plot_Overlay
   MLE_Plot
);