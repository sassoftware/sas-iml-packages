/*
   vdc: Generate the van der Corput sequence in the specified base.
   Inputs:
      base      - Integer radix >= 2.
      numPts - Positive integer number of points.
   Output: numPts x 1 column vector in [0,1), containing the radical
           inverses of integers 0 through numPts-1; the first point is 0.
   Depends: ConvertFromBase10(), numDigBase()
*/
start vdc(base, numPts);
   if numPts = 1 then return( 0 );
   seq = T(0:(numPts-1));
   digits = ConvertFromBase10(seq, base);
   numDigits = numDigBase(numPts-1, base);
   pow = -(numDigits:1);   /* reverse the weights instead of the digits */
   weights = base##pow;
   points = (digits # weights)[,+];   /* sum of weighted columns */
   return(points);
finish;

/*
start vdc_orig(base, numPts);
   if numPts = 1 then return( 0 );
   points = j(numPts, 1, 0);
   numDigits = numDigBase(numPts, base);
   pow = -(1:numDigits);
   weights = base##pow;
   digits = 0;
   points[1] = 0;
   do pointIndex = 2 to numPts;
      digits = nbe(digits, base);
      reversedDigits = digits[ , ncol(digits):1];
      points[pointIndex] = sum(reversedDigits # weights[,1:ncol(digits)]);
   end;
   return(points);
finish;
*/

store module=(vdc);