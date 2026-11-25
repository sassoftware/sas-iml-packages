/*
   sobol: Generate a base-2 Sobol sequence in binary-index order.
   Inputs:
      dim - Positive integer number of coordinates, at most 20
                  as supported by sobmat.
      numPts - Positive integer final sequence index, not the row count.
   Output: (numPts+1) x dim matrix in [0,1), containing sequence
           indices 0 through numPts, including the initial all-zero row.
*/
/*
   sobol: Generate a base-2 Sobol sequence in binary-index order.
   Inputs:
      dim - Positive integer number of coordinates, at most 20
                  as supported by sobmat.
      numPts - Positive integer final sequence index, not the row count.
   Output: (numPts+1) x dim matrix in [0,1), containing sequence
           indices 0 through numPts, including the initial all-zero row.
   Depends: cbe(), sobmat()
*/
start sobol(dim, numPts);
   base = 2;

   /* 1. Base-2 digit matrix.
      cbe() puts the most significant coefficients on the left, unlike the
      least-significant-first internal digits in the MATLAB sobol() from
      the Handbook of Monte Carlo Methods. Direction-block input columns
      below are reordered; output weights and sequence points are unchanged.
   */
   digits = cbe(T(0:numPts), base);
   if cmiss(digits[1]) then return(.);
   numDigits = ncol(digits);
   reverse = numDigits:1;

   /* 2. Column vector of weights (replaces repeat, elementwise #, and [, +]) */
   weights = base ## -T(1:numDigits);

   matrixBlocks = sobmat(dim, numDigits);
   points = j(numPts + 1, dim, 0);

   do dimIdx = 1 to dim;
      columns = (dimIdx - 1) * numDigits + (1:numDigits);
      directions = matrixBlocks[ , columns];
      directions = directions[ , reverse];
      points[ , dimIdx] = mod(digits * T(directions), base) * weights;
   end;
   return(points);
finish;

/*
start sobol_orig(dim, numPts);
   base = 2;
   numDigits = numDigBase(numPts, base);
   if cmiss(numDigits) then return(.);
   numDigits = numDigBase(numPts, base);
   weights = repeat(base##-(1:numDigits), numPts + 1, 1);
   points = j(numPts + 1, dim, 0);
   matrixBlocks = sobmat(dim, numDigits);
   if cmiss(matrixBlocks[1]) then return(.);
   digits = repeat(T(0:numPts), 1, numDigits);
   do digitIdx = 1 to numDigits;
      digits[ , digitIdx] = mod(digits[ , digitIdx], base);
      if digitIdx < numDigits then
         digits[ , (digitIdx + 1):numDigits] = floor(digits[ , (digitIdx + 1):numDigits] / base);
   end;
   do dimIdx = 1 to dim;
      columns = (dimIdx - 1) * numDigits + (1:numDigits);
      directions = matrixBlocks[ , columns];
      weightedDigits = weights # mod(digits * T(directions), base);
      points[ , dimIdx] = weightedDigits[ , +];
   end;
   return(points);
finish;
*/
store module=(sobol);