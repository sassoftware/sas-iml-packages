/*
   extkorobov: Generate an extended Korobov point set in base-2 sequence order.
   Inputs:
      multiplier - Positive integer lattice multiplier; use an odd value
                   for dyadic lattices.
      dim  - Positive integer number of coordinates.
      numPts  - Positive integer number of points; powers of 2 give
                   dyadic Korobov lattice point sets.
   Output: numPts x dim matrix in [0,1), ordered by the base-2
           van der Corput sequence, with the all-zero point first.
*/
start extkorobov(multiplier, dim, numPts);
   base = 2;
   generator = j(1, dim, 0);
   generator[1] = 1;
   do dimIdx = 2 to dim;
      generator[dimIdx] = mod(generator[dimIdx - 1] * multiplier, numPts);
   end;
   generators = repeat(generator, numPts, 1);
   sequence = vdc(base, numPts);
   indices = repeat(sequence, 1, dim);
   points = mod(indices # generators, 1);
   return(points);
finish;
store module=(extkorobov);