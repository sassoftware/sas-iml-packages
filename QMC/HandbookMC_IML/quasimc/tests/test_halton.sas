proc iml;
%include "QMC/HandbookMC_IML/quasimc/tests/check_test.sas";
%include "QMC/HandbookMC_IML/quasimc/nbe.sas";
%include "QMC/HandbookMC_IML/quasimc/vdc.sas";
%include "QMC/HandbookMC_IML/quasimc/halton.sas";

bases = {2 3 5};
points = halton(bases, 5);
correct = ({0, 4, 2, 6, 1} / 8) ||
          ({0, 3, 6, 1, 4} / 9) ||
          ({0, 1, 2, 3, 4} / 5);
run check_test("halton: first five three-dimensional points", points, correct);

actual = nrow(points) || ncol(points);
correct = {5 3};
run check_test("halton: output dimensions", actual, correct, 0);

actual = halton(T(bases), 5);
run check_test("halton: column-vector bases", actual, points);

actual = halton(bases, 1);
correct = j(1, 3, 0);
run check_test("halton: single point", actual, correct, 0);

actual = halton(2, 5);
correct = {0, 4, 2, 6, 1} / 8;
run check_test("halton: one dimension", actual, correct, 0);

actual = all(points >= 0 & points < 1);
correct = 1;
run check_test("halton: unit hypercube", actual, correct, 0);
quit;