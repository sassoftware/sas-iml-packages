/*
   halton: Generate a multidimensional Halton sequence.
   Inputs:
      bases     - Nonempty vector of integer bases >= 2, one per coordinate;
                  use pairwise coprime bases, typically distinct primes.
      numPoints - Positive integer number of points.
   Output: numPoints x d matrix in [0,1), where d is the number of bases.
           Each column is a van der Corput sequence; the first row is zero.
*/
start halton(bases, numPoints);
   dimension = nrow(bases) * ncol(bases);
   points = j(numPoints, dimension, 0);
   do dimIndex = 1 to dimension;
      points[ , dimIndex] = vdc(bases[dimIndex], numPoints);
   end;
   return(points);
finish;
store module=(halton);