/*
   vdc: Generate the van der Corput sequence in the specified base.
   Inputs:
      base      - Integer radix >= 2.
      numPoints - Positive integer number of points.
   Output: numPoints x 1 column vector in [0,1), containing the radical
           inverses of integers 0 through numPoints-1; the first point is 0.
*/
start vdc(base, numPoints);
   if numPoints = 1 then return( 0 );
   seq = T(0:(numPoints-1));
   digits = ConvertFromBase10(seq, base);
   numDigits = numDigBase(numPoints-1, base);
   pow = -(numDigits:1);   /* reverse the weights instead of the digits */
   weights = base##pow;
   points = (digits # weights)[,+];   /* sum of weighted columns */
   return(points);
finish;

start vdc_orig(base, numPoints);
   if numPoints = 1 then return( 0 );
   points = j(numPoints, 1, 0);
   numDigits = numDigBase(numPoints, base);
   pow = -(1:numDigits);
   weights = base##pow;
   digits = 0;
   points[1] = 0;
   do pointIndex = 2 to numPoints;
      digits = nbe(digits, base);
      reversedDigits = digits[ , ncol(digits):1];
      points[pointIndex] = sum(reversedDigits # weights[,1:ncol(digits)]);
   end;
   return(points);
finish;

store module=(vdc);