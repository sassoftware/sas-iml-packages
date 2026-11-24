/*
   halton: Generate a multidimensional Halton sequence.
   Inputs:
      bases     - Nonempty vector of integer bases >= 2;
                  Use pairwise coprime bases, typically the first d prime numbers.
      numPts - Positive integer number of points.
   Output: numPts x d matrix in [0,1), where d is the number of bases.
           Each column is a van der Corput sequence; the first row is zero.
   Depends: vdc()
*/
start halton(bases, numPts);
   dim = nrow(bases) * ncol(bases);
   points = j(numPts, dim, 0);
   do i = 1 to dim;
      points[ ,i] = vdc(bases[i], numPts);
   end;
   return(points);
finish;
store module=(halton);

