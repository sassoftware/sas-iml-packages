/*
   cbe: Convert a nonnegative integer to its positional digit expansion.
   Inputs:
      value - Nonnegative integer scalar to expand.
      base  - Integer radix >= 2.
   Output: Row vector of digits in the specified base, most significant first;
           zero returns {0}.
*/
start cbe(value, base);
   if value = 0 then return(0);
   numDigits = max(0, floor(log(value) / log(base))) + 1;
   digits = j(1, numDigits, 0);
   placeValue = base##(numDigits - 1);
   remainder = value;
   do digitIndex = 1 to numDigits;
      digits[digitIndex] = floor(remainder / placeValue);
      remainder = remainder - placeValue * digits[digitIndex];
      placeValue = placeValue / base;
   end;
   return(digits);
finish;
store module=(cbe);