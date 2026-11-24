/*
   faure: Generate a Faure low-discrepancy sequence using Pascal transforms.
   Inputs:
      base      - Prime integer >= 2 and >= dim.
      dim - Positive integer number of coordinates.
      numPts - Positive integer final sequence index, not the row count.
   Output: (numPts+1) x dim matrix in [0,1), containing sequence
           indices 0 through numPts, including the initial all-zero row.
   Depends: NONE
*/
/* vectorized version */
start faure(base, dim, numPts);
   numDigits = numDigBase(numPts, base);
   /* 1. Lower-triangular Pascal matrix mod base */
   idx = j(numDigits, numDigits);
   r = row(idx);
   c = col(idx);
   lowerMask = (r >= c);
   transform = lowerMask # mod(comb(r - 1, lowerMask # (c - 1)), base);

   /* 2. Base-b digit matrix  */
   y = floor(repeat(T(0:numPts), 1, numDigits) / (base ## (0:(numDigits - 1))));
   digits = mod(y, base);

   /* 3. Column vector of weights */
   weights = T(1 / (base ## (1:numDigits)));

   points = j(numPts + 1, dim, 0);
   points[ , 1] = digits * weights;
   do dimIdx = 2 to dim;
      digits = mod(digits * transform, base);
      points[ , dimIdx] = digits * weights;
   end;

   return(points);
finish;

/*
start faure_orig(base, dim, numPts);
   numDigits = numDigBase(numPts, base);
   transform = j(numDigits, numDigits, 0);
   do columnIndex = 1 to numDigits;
      do rowIndex = 1 to columnIndex;
         transform[rowIndex, columnIndex] =
            mod(comb(columnIndex - 1, rowIndex - 1), base);
      end;
   end;
   transform = T(transform);
   digits = repeat(T(0:numPts), 1, numDigits);
   do digitIdx = 1 to numDigits - 1;
      digits[ , digitIdx] = mod(digits[ , digitIdx], base);
      digits[ , (digitIdx + 1):numDigits] =
         floor(digits[ , (digitIdx + 1):numDigits] / base);
   end;
   weights = repeat(1 / (base##(1:numDigits)), numPts + 1, 1);
   weightedDigits = weights # digits;
   points = j(numPts + 1, dim, 0);
   points[ , 1] = weightedDigits[ , +];
   do dimIdx = 2 to dim;
      digits = mod(digits * transform, base);
      weightedDigits = weights # digits;
      points[ , dimIdx] = weightedDigits[ , +];
   end;
   return(points);
finish;
*/
store module=(faure);