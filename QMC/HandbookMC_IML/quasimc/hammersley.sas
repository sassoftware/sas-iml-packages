/*
   hammersley: Generate a Hammersley point set with the source point ordering.
   Inputs:
      bases     - Nonempty vector of integer bases >= 2; its length is the
                  dimension d. Only the first d-1 bases are used and should
                  be pairwise coprime; the final base is ignored.
      numPoints - Positive integer number of points.
   Output: numPoints x d matrix in [0,1). The first coordinate is
           k/numPoints for k=0,...,numPoints-1; the remaining coordinates
           prepend a zero row to the first numPoints-1 Halton rows.
*/
start hammersley(bases, numPoints);
   dimension = nrow(bases) * ncol(bases);
   points = j(numPoints, dimension, 0);
   if numPoints > 1 & dimension > 1 then
      points[2:numPoints, 2:dimension] =
         halton(bases[1:(dimension - 1)], numPoints - 1);
   points[ , 1] = T(0:(numPoints - 1)) / numPoints;
   return(points);
finish;
store module=(hammersley);