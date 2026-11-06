/*
   halton: Generate a multidimensional Halton sequence.
   Inputs:
      bases     - Nonempty vector of integer bases >= 2;
                  Use pairwise coprime bases, typically the first d prime numbers.
      numPoints - Positive integer number of points.
   Output: numPoints x d matrix in [0,1), where d is the number of bases.
           Each column is a van der Corput sequence; the first row is zero.
*/
start halton(bases, numPoints);
   dim = nrow(bases) * ncol(bases);
   points = j(numPoints, dim, 0);
   do i = 1 to dim;
      points[ ,i] = vdc(bases[i], numPoints);
   end;
   return(points);
finish;
store module=(halton);

