proc iml;
%include "QMC/HandbookMC_IML/quasimc/tests/check_test.sas";
%include "QMC/HandbookMC_IML/quasimc/faure.sas";

actual = faure(2, 2, 7);
correct = {0 0,
           4 4,
           2 6,
           6 2,
           1 5,
           5 1,
           3 3,
           7 7} / 8;
run check_test("faure: first eight binary points", actual, correct, 0);

points = faure(3, 3, 8);
correct = {0 0 0,
           3 3 3,
           6 6 6,
           1 4 7,
           4 7 1,
           7 1 4,
           2 8 5,
           5 2 8,
           8 5 2} / 9;
run check_test("faure: first nine ternary points", points, correct);

actual = faure(2, 1, 4);
correct = {0, 4, 2, 6, 1} / 8;
run check_test("faure: one dimension and a power boundary", actual, correct, 0);

actual = faure(2, 2, 1);
correct = {0 0, 1 1} / 2;
run check_test("faure: one-digit expansion", actual, correct, 0);

actual = all(points >= 0 & points < 1);
correct = 1;
run check_test("faure: unit hypercube", actual, correct, 0);
