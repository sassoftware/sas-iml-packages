proc iml;
%include "check_test.sas";
%include "../korobov.sas";

points = korobov(2, 3, 5);
correct = {1 2 4,
           2 4 3,
           3 1 2,
           4 3 1,
           0 0 0} / 5;
run check_test("korobov: five-point lattice", points, correct);

actual = nrow(points) || ncol(points);
correct = {5 3};
run check_test("korobov: output dimensions", actual, correct, 0);

actual = points[5, ];
correct = j(1, 3, 0);
run check_test("korobov: origin is last", actual, correct, 0);

actual = points[ , 1];
call sort(actual, 1);
correct = T(0:4) / 5;
run check_test("korobov: first-coordinate lattice mesh", actual, correct);

actual = korobov(2, 1, 4);
correct = {1, 2, 3, 0} / 4;
run check_test("korobov: one dimension", actual, correct, 0);

actual = korobov(2, 3, 1);
correct = j(1, 3, 0);
run check_test("korobov: single point", actual, correct, 0);

actual = korobov(1, 2, 4);
correct = {1 1, 2 2, 3 3, 0 0} / 4;
run check_test("korobov: unit multiplier", actual, correct, 0);

actual = all(points >= 0 & points < 1);
correct = 1;
run check_test("korobov: unit hypercube", actual, correct, 0);
quit;