proc iml;
%include "check_test.sas";
%include "../convertbase.sas";
%include "../cbe.sas";

actual = cbe(0, 2);
correct = {0};
run check_test("cbe: zero", actual, correct, 0);

value = 13;
base = 2;
actual = cbe(value, base);
correct = {1 1 0 1};
run check_test("cbe: binary expansion", actual, correct, 0);

actual = cbe(T(0:7), base);  /* 8x3 matrix of binary expansion */
correct = {0 0 0,
           0 0 1,
           0 1 0,
           0 1 1,
           1 0 0,
           1 0 1,
           1 1 0,
           1 1 1};
run check_test("cbe: binary expansion", actual, correct, 0);

actual = cbe({0, 8, 9}, 3);  /* 3x3 matrix of ternary expansion */
correct = {0 0 0,
           0 2 2,
           1 0 0};
run check_test("cbe: ternary expansion", actual, correct, 0);

actual = cbe(26, 3);
correct = {2 2 2};
run check_test("cbe: just below a base power", actual, correct, 0);

actual = cbe(27, 3);
correct = {1 0 0 0};
run check_test("cbe: exact base power", actual, correct, 0);

actual = cbe(5, 7);
correct = {5};
run check_test("cbe: one digit", actual, correct, 0);

actual = cbe(100203, 10);
correct = {1 0 0 2 0 3};
run check_test("cbe: decimal expansion with zeros", actual, correct, 0);

do base = 2 to 7;
   do value = 0 to 64;
      digits = cbe(value, base);
      weights = base##((ncol(digits) - 1):0);
      actual = sum(digits # weights);
      run check_test("cbe: reconstruct integer", actual, value, 0);
      actual = all(digits >= 0 & digits < base & digits = floor(digits));
      correct = 1;
      run check_test("cbe: valid base digits", actual, correct, 0);
   end;
end;
quit;