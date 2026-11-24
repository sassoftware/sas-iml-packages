proc iml;
%include "check_test.sas";
%include "../convertbase.sas";
%include "../cbe.sas";
%include "../sobmat.sas";

matrices = sobmat(3, 3);
correct = {1 0 0 1 1 1 1 0 1,
           0 1 0 0 1 0 0 1 1,
           0 0 1 0 0 1 0 0 1};
run check_test("sobmat: three packed direction matrices", matrices, correct, 0);

actual = sobmat(1, 5);
correct = I(5);
run check_test("sobmat: first dimension is identity", actual, correct, 0);

actual = sobmat(3, 1);
correct = {1 1 1};
run check_test("sobmat: single-digit blocks", actual, correct, 0);

matrices = sobmat(3, 9);
actual = matrices[ , 18];
correct = {1, 0, 0, 0, 0, 0, 0, 0, 1};
run check_test("sobmat: second-dimension recurrence", actual, correct, 0);

actual = matrices[ , 27];
correct = {1, 1, 1, 0, 1, 0, 0, 0, 1};
run check_test("sobmat: third-dimension recurrence", actual, correct, 0);

matrices = sobmat(20, 9);
actual = matrices[ , 180];
correct = {1, 1, 1, 0, 0, 0, 0, 1, 1};
run check_test("sobmat: twentieth-dimension recurrence", actual, correct, 0);

actual = all(matrices = 0 | matrices = 1);
correct = 1;
run check_test("sobmat: binary coefficients", actual, correct, 0);
quit;