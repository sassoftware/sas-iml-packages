/*
   sobol: Generate a base-2 Sobol sequence in binary-index order.
   Inputs:
      dimension - Positive integer number of coordinates, at most 20
                  as supported by sobmat.
      numPoints - Positive integer final sequence index, not the row count.
   Output: (numPoints+1) x dimension matrix in [0,1), containing sequence
           indices 0 through numPoints, including the initial all-zero row.
*/
start sobol(dimension, numPoints);
   base = 2;
   numDigits = floor(log(numPoints) / log(base)) + 1;
   weights = repeat(1 / (base##(1:numDigits)), numPoints + 1, 1);
   points = j(numPoints + 1, dimension, 0);
   matrixBlocks = sobmat(dimension, numDigits);
   digits = repeat(T(0:numPoints), 1, numDigits);
   do digitIndex = 1 to numDigits;
      digits[ , digitIndex] = mod(digits[ , digitIndex], base);
      if digitIndex < numDigits then
         digits[ , (digitIndex + 1):numDigits] =
            floor(digits[ , (digitIndex + 1):numDigits] / base);
   end;
   do dimIndex = 1 to dimension;
      columns = (dimIndex - 1) * numDigits + (1:numDigits);
      directions = matrixBlocks[ , columns];
      weightedDigits = weights # mod(digits * T(directions), base);
      points[ , dimIndex] = weightedDigits[ , +];
   end;
   return(points);
finish;
store module=(sobol);