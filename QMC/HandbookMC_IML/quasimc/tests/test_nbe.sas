proc iml;
%include "QMC/HandbookMC_IML/quasimc/tests/check_test.sas";
%include "QMC/HandbookMC_IML/quasimc/nbe.sas";

digits = {0};
actual = nbe(digits, 2);
correct = {1};
run check_test("nbe: increment zero", actual, correct, 0);

digits = {1 0 1 1};
original = digits;
actual = nbe(digits, 2);
correct = {1 1 0 0};
run check_test("nbe: binary carry", actual, correct, 0);
run check_test("nbe: input unchanged", digits, original, 0);

digits = {1 1 1};
actual = nbe(digits, 2);
correct = {1 0 0 0};
run check_test("nbe: binary overflow", actual, correct, 0);

digits = {1 2 8};
actual = nbe(digits, 10);
correct = {1 2 9};
run check_test("nbe: decimal increment", actual, correct, 0);

digits = {9 9};
actual = nbe(digits, 10);
correct = {1 0 0};
run check_test("nbe: decimal overflow", actual, correct, 0);

digits = {1 0 2 2};
actual = nbe(digits, 3);
correct = {1 1 0 0};
run check_test("nbe: ternary operation", actual, correct, 0);
quit;