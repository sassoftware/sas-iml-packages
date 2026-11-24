/*
   sobmat: Construct binary generating matrices for the Sobol sequence.
   Inputs:
      dim - Positive integer number of coordinates, at most 20
                  because of the stored initialization table.
      numDigits - Positive integer number of binary digits per coordinate.
   Output: numDigits x (numDigits*dim) binary matrix of horizontally
           concatenated square generating matrices, one per coordinate.
           The first block is the identity matrix.
   Depends: cbe()
   Note: BXOR is a built-in SAS function for the logical operation "Binary eXclusive OR"
         It uses the binary representation of two integers in the range
         [0, 2^32 - 1] to return a new integer whose binary representation has
         the i_th bit equal to 1 if the i_th bit of x or y is 1, but not both.
*/
start sobmat(dim, numDigits);
   polyNumbers = {1 3 7 11 13 19 25 37 59 47 61 55 41 67 97 91 109 103
                  115 131 193 137 145 143 241 157 185 167 229 171 213 191
                  253 203 211 239 247 285 369 299 425 301 361 333 357 351
                  501 355 397 391 451 463 487};
   initVals ={1 1 1 1 1 1 1 1,
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
   matrices = j(numDigits, numDigits * dim, 0);
   matrices[ , 1:numDigits] = I(numDigits);
   /* Note: All calls to cbe() use Base=2 */
   do dimIdx = 2 to dim;
      polynomial = polyNumbers[dimIdx];
      dirNumbers = j(1, max(8, numDigits), 0);
      dirNumbers[1:8] = initVals[dimIdx, ];
      coeffs = cbe(polynomial, 2);
      coeffs = coeffs[, 2:ncol(coeffs)];
      degree = ncol(coeffs);
      do digitIdx = 9 to numDigits;
         xorTotal = 0;
         do coeffIdx = 1 to degree;
            xorTotal = bxor(xorTotal, (2##coeffIdx) * coeffs[coeffIdx] * dirNumbers[digitIdx - coeffIdx]);
         end;
         dirNumbers[digitIdx] = bxor(xorTotal, dirNumbers[digitIdx - degree]);
      end;
      do digitIdx = 1 to numDigits;
         binaryDigits = cbe(dirNumbers[digitIdx], 2);
         digitCount = ncol(binaryDigits);
         rows = (digitIdx - digitCount + 1):digitIdx;
         column = (dimIdx - 1) * numDigits + digitIdx;
         matrices[rows, column] = T(binaryDigits);
      end;
   end;
   return(matrices);
finish;
store module=(sobmat);