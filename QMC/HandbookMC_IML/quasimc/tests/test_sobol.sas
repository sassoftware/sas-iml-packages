proc iml;
%include "QMC/HandbookMC_IML/quasimc/tests/check_test.sas";
%include "QMC/HandbookMC_IML/quasimc/cbe.sas";
%include "QMC/HandbookMC_IML/quasimc/sobmat.sas";
%include "QMC/HandbookMC_IML/quasimc/sobol.sas";

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

actual = all(longPoints >= 0 & longPoints < 1);
correct = 1;
run check_test("sobol: unit hypercube", actual, correct, 0);
quit;