proc iml;
%include "QMC/HandbookMC_IML/quasimc/tests/check_test.sas";
%include "QMC/HandbookMC_IML/quasimc/nbe.sas";
%include "QMC/HandbookMC_IML/quasimc/vdc.sas";
%include "QMC/HandbookMC_IML/quasimc/korobov.sas";
%include "QMC/HandbookMC_IML/quasimc/extkorobov.sas";

points = extkorobov(3, 3, 8);
correct = {0 0 0,
           4 4 4,
           2 6 2,
           6 2 6,
           1 3 1,
           5 7 5,
           3 1 3,
           7 5 7} / 8;
run check_test("extkorobov: van der Corput point order", points, correct, 0);

actual = nrow(points) || ncol(points);
correct = {8 3};
run check_test("extkorobov: output dimensions", actual, correct, 0);

actual = points[1, ];
correct = j(1, 3, 0);
run check_test("extkorobov: origin is first", actual, correct, 0);

actual = points;
call sort(actual, 1);
correct = korobov(3, 3, 8);
call sort(correct, 1);
run check_test("extkorobov: same dyadic lattice point set", actual, correct, 0);

longPoints = extkorobov(3, 3, 16);
actual = longPoints[1:8, ];
run check_test("extkorobov: dyadic extension preserves prefix", actual, points, 0);

actual = extkorobov(2, 2, 5);
correct = {0 0, 4 0, 2 4, 6 4, 1 2} / 8;
run check_test("extkorobov: non-dyadic source behavior", actual, correct, 0);

actual = extkorobov(3, 1, 5);
correct = {0, 4, 2, 6, 1} / 8;
run check_test("extkorobov: one dimension", actual, correct, 0);

actual = extkorobov(3, 3, 1);
correct = j(1, 3, 0);
run check_test("extkorobov: single point", actual, correct, 0);

actual = all(points >= 0 & points < 1);
correct = 1;
run check_test("extkorobov: unit hypercube", actual, correct, 0);
quit;