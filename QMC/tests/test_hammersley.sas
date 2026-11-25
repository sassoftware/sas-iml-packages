proc iml;
%include "check_test.sas";
%include "../convertbase.sas";
%include "../cbe.sas";
%include "../vdc.sas";
%include "../halton.sas";
%include "../hammersley.sas";

bases = {2 3 5};
points = hammersley(bases, 4);
correct = {0 0 0} //
          {0.25 0 0} //
          {0.5 0.5 0.333333333333} //
          {0.75 0.25 0.666666666667 };
run check_test("hammersley: source point ordering", points, correct);


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

numPoints = 1000;
bases = FirstPrimes(10);
h = hammersley(bases, numPoints);
actual = h[,1];
correct = T(0:(numPoints - 1)) / numPoints;
run check_test("hammersley: First column", actual, correct);

actual = h[1,];
correct = j(1, ncol(bases`), 0);
run check_test("hammersley: First column", actual, correct);

actual = h[3,];
correct = (2/numPoints) || (1/ rowvec(bases[1:nrow(bases)-1]));
run check_test("hammersley: First non-zero point", actual, correct, 0);

actual = h[100,];
correct = 99/1000 || 35/128 || 226/243 || 98/125 || 2/343 || 118/121 || 98/169 || 226/289 || 62/361 || 142/529 ;
run check_test("hammersley: 100th point", actual, correct);
quit;