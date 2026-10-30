proc iml;
%include "QMC/HandbookMC_IML/quasimc/tests/check_test.sas";
%include "QMC/HandbookMC_IML/quasimc/nbe.sas";
%include "QMC/HandbookMC_IML/quasimc/vdc.sas";

actual = vdc(2, 9);
correct = {0, 8, 4, 12, 2, 10, 6, 14, 1} / 16;
run check_test("vdc: first nine binary points", actual, correct, 0);

actual = vdc(3, 10);
correct = {0, 9, 18, 3, 12, 21, 6, 15, 24, 1} / 27;
run check_test("vdc: first ten ternary points", actual, correct);

actual = vdc(2, 1);
correct = {0};
run check_test("vdc: single point", actual, correct, 0);

points = vdc(2, 16);
actual = nrow(points) || ncol(points);
correct = {16 1};
run check_test("vdc: output dimensions", actual, correct, 0);

actual = floor(16 * points);
call sort(actual, 1);
correct = T(0:15);
run check_test("vdc: binary interval balance", actual, correct, 0);

actual = vdc(2, 8);
correct = points[1:8, ];
run check_test("vdc: sequence prefix", actual, correct, 0);
quit;