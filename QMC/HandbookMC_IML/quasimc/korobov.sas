/*
   korobov: Generate a rank-1 Korobov lattice point set.
   Inputs:
      multiplier - Positive integer lattice multiplier, preferably coprime
                   to numPts.
      dim  - Positive integer number of coordinates.
      numPts  - Positive integer number of points.
   Output: numPts x dim matrix in [0,1), using lattice indices
           1 through numPts, so the all-zero point is the final row.
*/
start korobov(multiplier, dim, numPts);
   generator = j(1, dim, 0);
   generator[1] = 1;
   do dimIdx = 2 to dim;
      generator[dimIdx] = mod(generator[dimIdx - 1] * multiplier, numPts);
   end;
   generators = repeat(generator, numPts, 1);
   indices = repeat(T(1:numPts), 1, dim);
   points = mod(indices # generators / numPts, 1);
   return(points);
finish;
store module=(korobov);