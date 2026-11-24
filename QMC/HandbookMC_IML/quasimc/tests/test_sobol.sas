proc iml;
%include "check_test.sas";
%include "../convertbase.sas";
%include "../cbe.sas";
%include "../sobmat.sas";
%include "../sobol.sas";

points = sobol(3, 7);
correct = {0 0 0,
           4 4 4,
           2 6 2,
           6 2 6,
           1 5 7,
           5 1 3,
           3 3 5,
           7 7 1} / 8;
run check_test("sobol: first eight source-ordered points", points, correct, 0);

actual = nrow(points) || ncol(points);
correct = {8 3};
run check_test("sobol: source returns N plus one points", actual, correct, 0);

actual = sobol(1, 1);
correct = {0, 1} / 2;
run check_test("sobol: one dimension and one digit", actual, correct, 0);

shortPoints = sobol(3, 5);
correct = points[1:6, ];
run check_test("sobol: non-power-of-two prefix", shortPoints, correct, 0);

longPoints = sobol(3, 15);
actual = longPoints[1:8, ];
run check_test("sobol: prefix across digit lengths", actual, points, 0);

do dimIndex = 1 to 3;
   actual = floor(16 * longPoints[ , dimIndex]);
   call sort(actual, 1);
   correct = T(0:15);
   run check_test("sobol: one-dimensional binary balance", actual, correct, 0);
end;

longPoints = sobol(3, 256);
actual = longPoints[257, ];
correct = {1 257 465} / 512;
run check_test("sobol: direction recurrence beyond eight digits", actual, correct, 0);

/* Test 8: 2D (0, 4, 2)-net balance on a 4 x 4 dyadic grid for N = 16 points */
pts16 = sobol(3, 15);
cellX = floor(4 * pts16[ , 1]); /* 0..3 */
cellY = floor(4 * pts16[ , 2]); /* 0..3 */
actual = 4 * cellX + cellY;     /* Unique cell index in 0..15 */
call sort(actual, 1);
correct = T(0:15);
run check_test("sobol: 2D (0,4,2)-net 4x4 grid balance", actual, correct, 0);

quit;