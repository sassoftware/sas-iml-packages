/*
   vdc: Generate the van der Corput sequence in the specified base.
   Inputs:
      base      - Integer radix >= 2.
      numPts - Positive integer number of points.
   Output: numPts x 1 column vector in [0,1), containing the radical
           inverses of integers 0 through numPts-1; the first point is 0.
   Depends: cbe()
*/
start vdc(base, numPts);
   if numPts = 1 then return( 0 );
   seq = T(0:(numPts-1));
   /* cbe() keeps the most significant coefficients on the left. The MATLAB
      vdc() from the Handbook of Monte Carlo Methods reverses its digit vector
      before weighting; here the weights are reversed instead, preserving the
      same radical inverses without reversing the digit columns.
   */
   digits = cbe(seq, base);
   if cmiss(digits[1]) then return(.);
   numDigits = ncol(digits);
   pow = -(numDigits:1);   /* reverse the weights instead of the digits */
   weights = base##pow;
   points = (digits # weights)[,+];   /* sum of weighted columns */
   return(points);
finish;
store module=(vdc);