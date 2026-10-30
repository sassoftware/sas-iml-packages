/*
   vdc: Generate the van der Corput sequence in the specified base.
   Inputs:
      base      - Integer radix >= 2.
      numPoints - Positive integer number of points.
   Output: numPoints x 1 column vector in [0,1), containing the radical
           inverses of integers 0 through numPoints-1; the first point is 0.
*/
start vdc(base, numPoints);
   points = j(numPoints, 1, 0);
   if numPoints = 1 then return(points);
   numDigits = ceil(log(numPoints) / log(base));
   weights = 1 / (base##(1:numDigits));
   digits = {};
   points[1] = 0;
   do pointIndex = 2 to numPoints;
      digits = nbe(digits, base);
      reversedDigits = digits[ , ncol(digits):1];
      points[pointIndex] = sum(reversedDigits # weights[1:ncol(digits)]);
   end;
   return(points);
finish;
store module=(vdc);