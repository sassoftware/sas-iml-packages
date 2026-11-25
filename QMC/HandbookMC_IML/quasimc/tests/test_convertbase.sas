proc iml;
%include "check_test.sas";
%include "../convertbase.sas";
%include "../cbe.sas";

exampleValues = {2, 3, 4, 5, 8, 9, 15, 16, 48};
originalValues = exampleValues;

actual = logbase(exampleValues, 3);
correct = log(exampleValues) / log(3);
run check_test("logbase: documented base-3 example", actual, correct);

actual = logbase({1 10 100 1000});
correct = {0 1 2 3};
run check_test("logbase: default base and row vector", actual, correct);

actual = logbase({1 2 4 8}, 2);
correct = {0 1 2 3};
run check_test("logbase: binary powers", actual, correct);

mixedValues = {-8 0 1 8 .};
mixedLogs = logbase(mixedValues, 2);
actual = (mixedLogs = .);
correct = {1 1 0 0 1};
run check_test("logbase: nonpositive and missing inputs", actual, correct, 0);

actual = mixedLogs[ , 3:4];
correct = {0 3};
run check_test("logbase: valid entries in mixed input", actual, correct);

actual = numDigBase(exampleValues, 3);
correct = {1, 2, 2, 2, 2, 3, 3, 3, 4};
run check_test("numDigBase: documented base-3 example", actual, correct, 0);

actual = numDigBase({1 9 10 99 100});
correct = {1 1 2 2 3};
run check_test("numDigBase: default base and powers", actual, correct, 0);

actual = numDigBase({7 8 9 26 27 28}, 3);
correct = {2 2 3 3 4 4};
run check_test("numDigBase: neighbors of base powers", actual, correct, 0);

actual = numDigBase({2.4 2.6 8.6}, 3);
correct = {1 2 3};
run check_test("numDigBase: rounded inputs", actual, correct, 0);

actual = convertFromBase10(exampleValues, 3);
correct = {0 0 0 2,
           0 0 1 0,
           0 0 1 1,
           0 0 1 2,
           0 0 2 2,
           0 1 0 0,
           0 1 2 0,
           0 1 2 1,
           1 2 1 0};
run check_test("convertFromBase10: documented example", actual, correct, 0);

actual = exampleValues;
correct = originalValues;
run check_test("convertbase: value inputs unchanged", actual, correct, 0);

actual = convertFromBase10(13, 2);
correct = {1 1 0 1};
run check_test("convertFromBase10: scalar binary expansion", actual, correct, 0);

actual = convertFromBase10(100203, 10);
correct = {1 0 0 2 0 3};
run check_test("convertFromBase10: decimal digits and zeros", actual, correct, 0);

actual = convertFromBase10({2.4, 2.6, 8.6}, 3);
correct = {0 0 2, 0 1 0, 1 0 0};
run check_test("convertFromBase10: rounded column inputs", actual, correct, 0);

exampleDigits = {0 1 2 1,
                 1 2 1 0};
originalDigits = exampleDigits;
actual = convertToBase10(exampleDigits, 3);
correct = {16, 48};
run check_test("convertToBase10: documented base-3 example", actual, correct, 0);

actual = exampleDigits;
correct = originalDigits;
run check_test("convertToBase10: digit input unchanged", actual, correct, 0);

actual = convertToBase10({1 1 0 1}, 2);
correct = 13;
run check_test("convertToBase10: binary row vector", actual, correct, 0);

actual = convertToBase10({15 15}, 16);
correct = 255;
run check_test("convertToBase10: hexadecimal digits", actual, correct, 0);

actual = convertToBase10({0 0 0, 0 0 1}, 7);
correct = {0, 1};
run check_test("convertToBase10: zero and leading zeros", actual, correct, 0);

safe_cutoff = 1E9;
values = T(0:64);
do base = 2 to 16;
   digits = convertFromBase10(values, base);
   actual = convertToBase10(digits, base);
   correct = values;
   run check_test("convertbase: round trip for bases 2 through 16", actual, correct, 0);

   actual = all(digits >= 0 & digits < base & digits = floor(digits));
   correct = 1;
   run check_test("convertFromBase10: valid integer digits", actual, correct, 0);

   power = base;
   powerDigits = 2;
   do while (power <= safe_cutoff);
      boundaryValues = (power - 1) // power;
      boundaryCounts = (powerDigits - 1) // powerDigits;
      if power < safe_cutoff then do;
         boundaryValues = boundaryValues // (power + 1);
         boundaryCounts = boundaryCounts // powerDigits;
      end;
      actual = numDigBase(boundaryValues, base);
      run check_test("numDigBase: base-power boundaries below cutoff", actual, boundaryCounts, 0);
      power = power * base;
      powerDigits = powerDigits + 1;
   end;
end;

actual = numDigBase(0, 3);
correct = 1;
run check_test("numDigBase: zero has one digit", actual, correct, 0);

actual = numDigBase((safe_cutoff - 1) || safe_cutoff, 10);
correct = {9 10};
run check_test("numDigBase: cutoff and predecessor", actual, correct, 0);

cutoffValues = (safe_cutoff - 1) // safe_cutoff;
cutoffDigits = convertFromBase10(cutoffValues, 10);
correct = {0 9 9 9 9 9 9 9 9 9,
           1 0 0 0 0 0 0 0 0 0};
run check_test("convertFromBase10: cutoff is included", cutoffDigits, correct, 0);

actual = convertToBase10(cutoffDigits, 10);
correct = cutoffValues;
run check_test("convertbase: cutoff boundary round trip", actual, correct, 0);

overLimitValues = 1 // (safe_cutoff + 1);
actual = numDigBase(overLimitValues, 10);
actual = (actual = .);
correct = 1;
run check_test("numDigBase: reject whole over-limit call", actual, correct, 0);

actual = convertFromBase10(overLimitValues, 10);
actual = (actual = .);
run check_test("convertFromBase10: reject whole over-limit call", actual, correct, 0);

overLimitDigits = {0 0 0 0 0 0 0 0 0 1,
                  1 0 0 0 0 0 0 0 0 1};
actual = convertToBase10(overLimitDigits, 10);
actual = (actual = .);
run check_test("convertToBase10: reject whole over-limit call", actual, correct, 0);

bigint = constant('exactint');
exactBits = round(log2(bigint));
actual = numDigBase(bigint, 2);
actual = (actual = .);
run check_test("numDigBase: reject EXACTINT", actual, correct, 0);

actual = convertFromBase10(bigint, 2);
actual = (actual = .);
run check_test("convertFromBase10: reject EXACTINT", actual, correct, 0);

actual = convertToBase10(1 || j(1, exactBits, 0), 2);
actual = (actual = .);
run check_test("convertToBase10: reject EXACTINT", actual, correct, 0);

actual = numDigBase(bigint - 1, 2);
actual = (actual = .);
run check_test("numDigBase: reject EXACTINT predecessor", actual, correct, 0);

actual = convertFromBase10(bigint - 1, 2);
actual = (actual = .);
run check_test("convertFromBase10: reject EXACTINT predecessor", actual, correct, 0);

actual = convertToBase10(j(1, exactBits, 1), 2);
actual = (actual = .);
run check_test("convertToBase10: reject EXACTINT predecessor", actual, correct, 0);

negativeLogs = logbase({-8 -1}, 2);
actual = (negativeLogs = .);
correct = {1 1};
run check_test("logbase: negative-only input returns missing", actual, correct, 0);
quit;