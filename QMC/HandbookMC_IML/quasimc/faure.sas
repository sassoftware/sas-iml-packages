/*
   faure: Generate a Faure low-discrepancy sequence using Pascal transforms.
   Inputs:
      base      - Prime integer >= 2 and >= dimension.
      dimension - Positive integer number of coordinates.
      numPoints - Positive integer final sequence index, not the row count.
   Output: (numPoints+1) x dimension matrix in [0,1), containing sequence
           indices 0 through numPoints, including the initial all-zero row.
*/
start faure(base, dimension, numPoints);
   numDigits = floor(log(numPoints) / log(base)) + 1;
   weights = repeat(1 / (base##(1:numDigits)), numPoints + 1, 1);
   points = j(numPoints + 1, dimension, 0);
   transform = j(numDigits, numDigits, 0);
   do columnIndex = 1 to numDigits;
      do rowIndex = 1 to columnIndex;
         transform[rowIndex, columnIndex] =
            mod(comb(columnIndex - 1, rowIndex - 1), base);
      end;
   end;
   transform = T(transform);
   digits = repeat(T(0:numPoints), 1, numDigits);
   do digitIndex = 1 to numDigits - 1;
      digits[ , digitIndex] = mod(digits[ , digitIndex], base);
      digits[ , (digitIndex + 1):numDigits] =
         floor(digits[ , (digitIndex + 1):numDigits] / base);
   end;
   weightedDigits = weights # digits;
   points[ , 1] = weightedDigits[ , +];
   do dimIndex = 2 to dimension;
      digits = mod(digits * transform, base);
      weightedDigits = weights # digits;
      points[ , dimIndex] = weightedDigits[ , +];
   end;
   return(points);
finish;
store module=(faure);