/*
   korobov: Generate a rank-1 Korobov lattice point set.
   Inputs:
      multiplier - Positive integer lattice multiplier, preferably coprime
                   to numPoints.
      dimension  - Positive integer number of coordinates.
      numPoints  - Positive integer number of points.
   Output: numPoints x dimension matrix in [0,1), using lattice indices
           1 through numPoints, so the all-zero point is the final row.
*/
start korobov(multiplier, dimension, numPoints);
   generator = j(1, dimension, 0);
   generator[1] = 1;
   do dimIndex = 2 to dimension;
      generator[dimIndex] = mod(generator[dimIndex - 1] * multiplier, numPoints);
   end;
   generators = repeat(generator, numPoints, 1);
   indices = repeat(T(1:numPoints), 1, dimension);
   points = mod(indices # generators / numPoints, 1);
   return(points);
finish;
store module=(korobov);