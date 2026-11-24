proc iml;
%include "check_test.sas";
%include "../convertbase.sas";
%include "../nbe.sas";
%include "../vdc.sas";
%include "../halton.sas";

bases = {2 3 5};
points = halton(bases, 5);
correct = ({0, 4, 2, 6, 1} / 8) ||
          ({0, 3, 6, 1, 4} / 9) ||
          ({0, 1, 2, 3, 4} / 5);
run check_test("halton: first five three-dimensional points", points, correct);

actual = nrow(points) || ncol(points);
correct = {5 3};
run check_test("halton: output dimensions", actual, correct, 0);


actual = halton(bases, 1);
correct = j(1, 3, 0);
run check_test("halton: single point", actual, correct, 0);

actual = halton(2, 5);
correct = {0, 4, 2, 6, 1} / 8;
run check_test("halton: one dimension", actual, correct, 0);

actual = all(points >= 0 & points < 1);
correct = 1;
run check_test("halton: unit hypercube", actual, correct, 0);

bases = FirstPrimes(10);
h = halton(bases, 1000);
actual = h[2,];
correct = 1 / bases`;
run check_test("halton: First non-zero point", actual, correct, 0);

actual = h[100,];
correct = 99/128 || 19/243 || 123/125 || 51/343 || 9/121 || 111/169 || 243/289 || 81/361 || 165/529 || 351/841;
run check_test("halton: 100th point", actual, correct, 1E-14);

/*
title "Halton Sequence: bases = (2,3)";
call scatter(h[,1], h[,2]) grid={x y};

title "Halton Sequence: bases = (23,29)";
call scatter(h[,9], h[,10]) grid={x y};
*/
