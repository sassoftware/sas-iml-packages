/* cbe: Convert nonnegative integer(s) to positional base-b digit expansion.
   Inputs:
      value - Scalar or vector of nonnegative integers (<= 1E9).
      base  - Integer radix >= 2.
   Output: Matrix of base-b digits (most significant first), with one row per
           element of value and columns padded on the left with zeros.
           Scalar input returns a 1 x k row vector; 0 returns {0}.
   Depends: numDigBase
*/
start cbe(value, base);
   v = colvec(value);
   numDigits = max(numDigBase(v, base));
   if cmiss(numDigits) then return(.);

   /* Place values from most significant (base^(k-1)) to least significant (1) */
   Numer = repeat(v, 1, numDigits);
   Denom = base ## ((numDigits - 1):0);
   digits = mod(floor(Numer / Denom), base);
   return(digits);
finish;

/*
start cbe_orig(value, base);
   if value = 0 then return(0);
   numDigits = max(0, floor(log(value) / log(base))) + 1;
   digits = j(1, numDigits, 0);
   placeValue = base##(numDigits - 1);
   remainder = value;
   do digitIdx = 1 to numDigits;
      digits[digitIdx] = floor(remainder / placeValue);
      remainder = remainder - placeValue * digits[digitIdx];
      placeValue = placeValue / base;
   end;
   return(digits);
finish;
*/
store module=(cbe);