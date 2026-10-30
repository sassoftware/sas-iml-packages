proc iml;
%include "QMC/HandbookMC_IML/quasimc/tests/check_test.sas";
%include "QMC/HandbookMC_IML/quasimc/nbe.sas";
%include "QMC/HandbookMC_IML/quasimc/vdc.sas";
%include "QMC/HandbookMC_IML/quasimc/halton.sas";
%include "QMC/HandbookMC_IML/quasimc/hammersley.sas";

bases = {2 3 5};
points = hammersley(bases, 4);
correct = ({0, 1, 2, 3} / 4) ||
          ({0, 0, 2, 1} / 4) ||
          ({0, 0, 1, 2} / 3);
run check_test("hammersley: source point ordering", points, correct);

actual = nrow(points) || ncol(points);
correct = {4 3};
run check_test("hammersley: output dimensions", actual, correct, 0);

actual = hammersley(2, 4);
correct = {0, 1, 2, 3} / 4;
run check_test("hammersley: one dimension", actual, correct, 0);

actual = hammersley(bases, 1);
correct = j(1, 3, 0);
run check_test("hammersley: single point", actual, correct, 0);

actual = hammersley({2 3 7}, 4);
run check_test("hammersley: last base unused in source", actual, points);

actual = all(points >= 0 & points < 1);
correct = 1;
run check_test("hammersley: unit hypercube", actual, correct, 0);
quit;