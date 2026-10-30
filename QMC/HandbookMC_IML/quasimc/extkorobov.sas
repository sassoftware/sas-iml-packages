/*
   extkorobov: Generate an extended Korobov point set in base-2 sequence order.
   Inputs:
      multiplier - Positive integer lattice multiplier; use an odd value
                   for dyadic lattices.
      dimension  - Positive integer number of coordinates.
      numPoints  - Positive integer number of points; powers of 2 give
                   dyadic Korobov lattice point sets.
   Output: numPoints x dimension matrix in [0,1), ordered by the base-2
           van der Corput sequence, with the all-zero point first.
*/
start extkorobov(multiplier, dimension, numPoints);
   base = 2;
   generator = j(1, dimension, 0);
   generator[1] = 1;
   do dimIndex = 2 to dimension;
      generator[dimIndex] = mod(generator[dimIndex - 1] * multiplier, numPoints);
   end;
   generators = repeat(generator, numPoints, 1);
   sequence = vdc(base, numPoints);
   indices = repeat(sequence, 1, dimension);
   points = mod(indices # generators, 1);
   return(points);
finish;
store module=(extkorobov);