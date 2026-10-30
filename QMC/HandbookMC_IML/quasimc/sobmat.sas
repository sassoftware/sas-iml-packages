/*
   sobmat: Construct binary generating matrices for the Sobol sequence.
   Inputs:
      dimension - Positive integer number of coordinates, at most 20
                  because of the stored initialization table.
      numDigits - Positive integer number of binary digits per coordinate.
   Output: numDigits x (numDigits*dimension) binary matrix of horizontally
           concatenated square generating matrices, one per coordinate;
           the first block is the identity matrix.
*/
start sobmat(dimension, numDigits);
   polynomialNumbers = {1 3 7 11 13 19 25 37 59 47 61 55 41 67 97 91 109 103
                        115 131 193 137 145 143 241 157 185 167 229 171 213 191
                        253 203 211 239 247 285 369 299 425 301 361 333 357 351
                        501 355 397 391 451 463 487};
   initialValues = {1 1 1 1 1 1 1 1,
                   1 3 5 15 17 51 85 255,
                   1 1 7 11 13 61 67 79,
                   1 3 7 5 7 43 49 147,
                   1 1 5 3 15 51 125 141,
                   1 3 1 1 9 59 25 89,
                   1 1 3 7 31 47 109 173,
                   1 3 3 9 9 57 43 43,
                   1 3 7 13 3 35 89 9,
                   1 1 5 11 27 53 69 25,
                   1 3 5 1 15 19 113 115,
                   1 1 7 3 29 51 47 97,
                   1 3 7 7 21 61 55 19,
                   1 1 1 9 23 39 97 97,
                   1 3 3 5 19 33 3 197,
                   1 1 3 13 11 7 37 101,
                   1 1 7 13 25 5 83 255,
                   1 3 5 11 7 11 103 29,
                   1 1 1 3 13 39 27 203,
                   1 3 1 15 17 63 13 65};
   matrices = j(numDigits, numDigits * dimension, 0);
   matrices[ , 1:numDigits] = I(numDigits);
   do dimIndex = 2 to dimension;
      polynomial = polynomialNumbers[dimIndex];
      directionNumbers = j(1, max(8, numDigits), 0);
      directionNumbers[1:8] = initialValues[dimIndex, ];
      coefficients = cbe(polynomial, 2);
      coefficients = coefficients[2:ncol(coefficients)];
      degree = ncol(coefficients);
      do digitIndex = 9 to numDigits;
         xorTotal = 0;
         do coefficientIndex = 1 to degree;
            xorTotal = bxor(xorTotal, (2##coefficientIndex) *
               coefficients[coefficientIndex] *
               directionNumbers[digitIndex - coefficientIndex]);
         end;
         directionNumbers[digitIndex] =
            bxor(xorTotal, directionNumbers[digitIndex - degree]);
      end;
      do digitIndex = 1 to numDigits;
         binaryDigits = cbe(directionNumbers[digitIndex], 2);
         digitCount = ncol(binaryDigits);
         rows = (digitIndex - digitCount + 1):digitIndex;
         column = (dimIndex - 1) * numDigits + digitIndex;
         matrices[rows, column] = T(binaryDigits);
      end;
   end;
   return(matrices);
finish;
store module=(sobmat);