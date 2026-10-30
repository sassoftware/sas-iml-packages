/* Utility functions to represent an integer in various bases b >= 2.
   https://blogs.sas.com/content/iml/2025/11/03/convert-integers-bases.html
*/

/* Compute LOG_b(x) for a vector of x values and for an integer base, b > 0.
   If an input value is not positive, this function returns a missing value.
   Example: 
     t = {2,3,4,5,8,9,15,16,48};
     log3 = logbase(t, 3);
*/
start logbase(x, base=10);
   if all(x>0) then 
      return log2(x) / log2(base);
   /* otherwise, return missing values for x <= 0 */
   y = j(nrow(x), ncol(x), .);
   idx = loc(x>0);
   if ncol(idx)>0 then 
      y[idx] = log2(x[idx]) / log2(base);
   return y;
finish;

/* This function returns the number of digits in the base-b representation of an integer, x.
   If n >= 0 is a base-10 integer, it has
      k = ceil(log_b(n+1))
   digits when represented in base b. See
   https://blogs.sas.com/content/iml/2015/08/31/digits-in-integer.html
   If any input exceeds safe_cutoff (default 1E9), return a scalar missing value.
   Example:
      t = {2,3,4,5,8,9,15,16,48};
      numDig3 = numDigBase(t, 3);
*/
start numDigBase(x, base=10);
   safe_cutoff = 1E9;
   if any(x > safe_cutoff) then return(.);
   y = j(nrow(x), ncol(x), 1);
   n = round(x);              /* ensure argument is an integer */
   idx = loc(n>0);
   if ncol(idx)>0 then 
      y[idx] = ceil( logbase(n[idx]+1,base) );
   return y;
finish;
 
/* Convert integer x > 0 to a row vector in base b.
   If x is a column vector, return a matrix where each row is the base-b representation of x[i].
   The most significant bit is to the left; the least significant bit is to the right.
   For example, n=15 and base=3 gives (120)_3 because 15 = 1*3##2 + 2*3##1 + 0*3##0.
   If any input exceeds safe_cutoff (default 1E9), return a scalar missing value.
   Example:
      t = {2,3,4,5,8,9,15,16,48};
      base3 = convertFromBase10(t, 3);
*/
start convertFromBase10(x, base);
   safe_cutoff = 1E9;
   if any(x > safe_cutoff) then return(.);
   n = round(x);     /* ensure inputs are integers */
   numDig = numDigBase(max(n), base);
   /* For an explanation, see https://blogs.sas.com/content/iml/2015/08/31/digits-in-integer.html */
   c = j(nrow(n), numDig, 0);
   do i = 1 to numDig;
      a = mod(n, base);
      n = floor( n/base );
      c[ , numDig - i + 1] = a;
   end;
   return c;
finish;

/* Convert a row vector from any base to a number in base 10.
   If c is a kxm matrix, then return a base-10 number fo each row.
   For example, in base 3, define:
   c = {0 1 2 1,   
        1 2 1 0 };
   The first row represents 16. The second row represents 48.
   See https://blogs.sas.com/content/iml/2011/11/16/converting-from-base-2-to-base-10.html 
   If any result exceeds safe_cutoff (default 1E9), return a scalar missing value.
   Example: 
      b10 = convertToBase10(c, 3);
*/
start convertToBase10(c, base);
   safe_cutoff = 1E9;
   pow = (ncol(c)-1):0;
   factor = base##pow;
   base10 = (c # factor)[ ,+];  /* c[k-1]*b##(k-1) + ... + c[1]*b##1 + c[0]*b##0 */ 
   if any(base10 > safe_cutoff) then return(.);
   return base10;
finish;

store module=(logbase numDigBase convertFromBase10 convertToBase10);
