/*
   nbe: Increment an integer represented by a vector in base b.
   Inputs:
   x - Row vector of integer digits in base b
       The digits must be 0..base-1, with the most significant
       digit first.
   base - Integer base >= 2.
   Output: Row vector representing x+1 in base b.
*/
start nbe(x, base);
   if any(x >= base) | any(x < 0) then return( . ); /* invalid input */
   numDigits = ncol(x);
   xp1 = x;
   carry = 1;
   do digitIndex = numDigits to 1 by -1;
      if carry then do;
         if x[digitIndex] ^= base - 1 then
            carry = 0;
         xp1[digitIndex] = mod(x[digitIndex] + 1, base);
      end;
   end;
   if carry then xp1 = 1 || xp1;
   return(xp1);
finish;
store module=(nbe);
QUIT;
