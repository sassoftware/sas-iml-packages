/*
   nbe: Increment an integer represented by a vector in base b.
   Inputs:
   x - Row vector of integer digits in base b
       The digits must be 0..base-1, with the most significant
       digit first.
   base - Integer base >= 2.
   Output: Row vector representing x+1 in base b.
   Depends: NONE
*/
start nbe(x, base);
   if any(x >= base) | any(x < 0) then return( . ); /* invalid input */
   numDigits = ncol(x);
   xp1 = x;
   carry = 1;
   do digitIdx = numDigits to 1 by -1;
      if carry then do;
         if x[digitIdx] ^= base - 1 then
            carry = 0;
         xp1[digitIdx] = mod(x[digitIdx] + 1, base);
      end;
   end;
   if carry then xp1 = 1 || xp1;
   return(xp1);
finish;
store module=(nbe);

