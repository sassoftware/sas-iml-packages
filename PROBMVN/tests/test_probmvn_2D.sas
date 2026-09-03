/*
IML tests that call the PROBMVN for the following situations:
1. 2-D problems. Include finite and semi-infinite domains of integration.
2. Degenerate problem. Create 3-D and 4-D examples that have a coordinate where the domain of integration is (-Infty, Infty), and thus it reduces to a 0-D, 1-D, or 2-D problems. Again, include both finite and semi-infinite limits for the variables that remain.
3. Degenerate problems that reduce to a diagonal matrix when the (-Infty, Infity) coordinates are excluded.
*/

/* Bivariate normal probabilities on rectangular domains. */
options ps=32000 nodate nonumber;

proc iml;
load module = _all_;

/*==============================*/
/* 1. 2-D test cases            */
/*==============================*/
print "Tests for 2-D PROBNRM and Extra Degenerate Dimensions";

/* 1a. 2-D finite rectangle, independent variables */
testName = "Test 1a: 2-D finite rectangle, R=I";
R = I(2);
L = {-1.0 -0.5};
U = { 0.8  1.2};
prob = probmvn_mod(L, U, R);
correct = (cdf("Normal", U[1]) - cdf("Normal", L[1])) #
          (cdf("Normal", U[2]) - cdf("Normal", L[2]));
run check_test(testName, prob, correct, 1E-12);

/* 1b. 2-D semi-infinite, positive orthant formula
   P(X>0, Y>0) = 1/4 + asin(rho)/(2*pi) */
testName = "Test 1b: 2-D negative orthant";
rho = 0.6;
R = I(2);
R[{2 3}] = rho;
U = {0 0};
L = {. .};
prob = probmvn_mod(L, U, R);
correct = 1/4 + arsin(rho)/(2*constant("pi"));
run check_test(testName, prob, correct, 1E-12);

/* 1c. 2-D semi-infinite, mixed orthant formula
   P(X<0, Y>0) = 1/4 - asin(rho)/(2*pi) */
testName = "Test 1c: 2-D mixed orthant";
rho = 0.4;
R = I(2);
R[{2 3}] = rho;
L = {. 0};
U = {0 .};
prob = probmvn_mod(L, U, R);
correct = 1/4 - arsin(rho)/(2*constant("pi"));
run check_test(testName, prob, correct, 1E-12);

/* 1d. 2-D semi-infinite, positive orthant formula
   P(X>0, Y>0) = 1/4 + asin(rho)/(2*pi) */
testName = "Test 1d: 2-D positive orthant";
rho = 0.6;
R = I(2);
R[{2 3}] = rho;
L = {0 0};
U = {. .};
prob = probmvn_mod(L, U, R);
correct = 1/4 + arsin(rho)/(2*constant("pi"));
run check_test(testName, prob, correct, 1E-12);


/*============================================*/
/* 2. Degenerate 3-D/4-D -> 0-D/1-D/2-D cases */
/*============================================*/

/* 2a. 3-D all coordinates unbounded -> 0-D problem */
testName = "Test 2a: 3-D all unbounded -> 0-D";
R = {1 0.3 -0.2,
     0.3 1 0.4,
    -0.2 0.4 1};
L = {. . .};
U = {. . .};
prob = probmvn_mod(L, U, R);
correct = 1;
run check_test(testName, prob, correct, 1E-12);

/* 2b. 3-D with two unbounded coordinates -> 1-D finite interval */
testName = "Test 2b: 3-D -> 1-D finite";
R = {1 0.3 -0.2,
     0.3 1 0.4,
    -0.2 0.4 1};
L = {. -0.75 .};
U = {.  1.10 .};
prob = probmvn_mod(L, U, R);
correct = cdf("Normal", U[2]) - cdf("Normal", L[2]);
run check_test(testName, prob, correct, 1E-12);

/* 2c. 3-D with two unbounded coordinates -> 1-D semi-infinite interval */
testName = "Test 2c: 3-D -> 1-D semi-infinite";
R = {1 0.3 -0.2,
     0.3 1 0.4,
    -0.2 0.4 1};
L = {. 0.20 .};
U = {. . .};
prob = probmvn_mod(L, U, R);
correct = 1 - cdf("Normal", L[2]);
run check_test(testName, prob, correct, 1E-12);

/* 2d. 4-D with two unbounded coordinates -> 2-D finite correlated rectangle */
testName = "Test 2d: 4-D -> 2-D finite correlated";
R = {1    0.2  0.35  0.0,
     0.2  1    0.1   0.4,
     0.35 0.1  1    -0.2,
     0.0  0.4 -0.2   1};
L = {-1.0 . -0.5 .};
U = { 0.7 .  1.2 .};
prob = probmvn_mod(L, U, R);
rho = R[1,3];
correct = probbnrm(U[1], U[3], rho) - probbnrm(L[1], U[3], rho)
        - probbnrm(U[1], L[3], rho) + probbnrm(L[1], L[3], rho);
run check_test(testName, prob, correct, 1E-12);

/* 2e. 4-D with two unbounded coordinates -> 2-D mixed semi-infinite */
testName = "Test 2e: 4-D -> 2-D mixed semi-infinite";
R = {1    0.2  0.35  0.0,
     0.2  1    0.1   0.4,
     0.35 0.1  1    -0.2,
     0.0  0.4 -0.2   1};
L = {0.0 . -1.0 .};
U = { .  .  0.8 .};
prob = probmvn_mod(L, U, R);
rho = R[1,3];
/* P(X1>0, -1<X3<0.8) = P(X3<0.8)-P(X1<0,X3<0.8)-P(X3<-1)+P(X1<0,X3<-1) */
correct = cdf("Normal", U[3]) - probbnrm(0, U[3], rho)
        - cdf("Normal", L[3]) + probbnrm(0, L[3], rho);
run check_test(testName, prob, correct, 1E-12);


/*==============================================================*/
/* 3. Degenerate -> reduced correlation matrix is diagonal      */
/*==============================================================*/

/* Matrix has nonzero correlations only with coordinate 4.
   If coord 4 is unbounded, reduced 3x3 matrix is identity. */
Sigma = {9    0     0    0.30,
         0    4     0    0.20,
         0    0     1   -0.25,
         0.30 0.20 -0.25 1};
mu = {2 -1 1 0};             /* mean of each effective var */
sd = rowvec(sqrt(vecdiag(Sigma)));   /* stddev of each effective var */
/* 3a. Reduced diagonal case with finite bounds on retained coordinates */
testName = "Test 3a: Reduced diagonal finite";
L = {-1.0 -0.5 -0.25 .};
U = { 0.8  0.7  1.10 .};
prob = probmvn_mod(L, U, Sigma, mu);
stdU = (U - mu) / sd;
stdL = (L - mu) / sd;
correct = (cdf("Normal", stdU[1]) - cdf("Normal", stdL[1])) #
          (cdf("Normal", stdU[2]) - cdf("Normal", stdL[2])) #
          (cdf("Normal", stdU[3]) - cdf("Normal", stdL[3]));
run check_test(testName, prob, correct, 1E-12);

/* 3b. Reduced diagonal case with mixed finite/semi-infinite retained bounds */
testName = "Test 3b: Reduced diagonal mixed bounds";
L = {.   0.1 -0.4 .};
U = {0.2  .   0.6 .};
prob = probmvn_mod(L, U, Sigma, mu);
stdU = (U - mu) / sd;
stdL = (L - mu) / sd;
correct = cdf("Normal", stdU[1])       #
          (1 - cdf("Normal", stdL[2])) #
          (cdf("Normal", stdU[3]) - cdf("Normal", stdL[3]));
run check_test(testName, prob, correct, 1E-12);

/*==============================================================*/
/* 4. Test finite truncation to (-delta, detla) for numerical   */
/*    stability, where delta ~ 8.125 is chosen so that          */
/*    SDF("Normal", delta) ~ constant("maceps")                 */
/*==============================================================*/
testName = "Test 4a: Right truncation onto (a, delta)";
L = {0  4};
U = {5  9};
Sigma = {1 -0.3,
         -0.3 1};
prob = probmvn_mod(L, U, Sigma);
clip_L = ClipInterval(L);
clip_U = ClipInterval(U);
correct = probmvn_mod(clip_L, clip_U, Sigma);
run check_test(testName, prob, correct, 1E-12);

testName = "Test 4b: Double right truncation collapse onto (delta, delta)";
L = {0  9};
U = {5  10};
Sigma = {1 -0.3,
         -0.3 1};
prob = probmvn_mod(L, U, Sigma);
correct = 0.0;
run check_test(testName, prob, correct, 1E-12);

testName = "Test 4c: Left truncation onto (-delta, b)";
L = {-9 0};
U = { 4 5};
Sigma = {1 -0.6,
         -0.6 1};
prob = probmvn_mod(L, U, Sigma);
clip_L = ClipInterval(L);
clip_U = ClipInterval(U);
correct = probmvn_mod(clip_L, clip_U, Sigma);
run check_test(testName, prob, correct, 1E-12);

testName = "Test 4d: Double left truncation collapse onto (-delta, -delta)";
L = {-10  -10};
U = {-9    -9};
Sigma = {1 -0.9,
         -0.9 1};
prob = probmvn_mod(L, U, Sigma);
correct = 0.0;
run check_test(testName, prob, correct, 1E-12);

testName = "Test 4e: Infinitely wide; infinitely thin";
L = {-9   9};
U = {10  10};
Sigma = {1 -0.9,
         -0.9 1};
prob = probmvn_mod(L, U, Sigma);
correct = 0.0;
run check_test(testName, prob, correct, 1E-12);


/***********************************************************/
/* Test 5: 2-D Test Cases: All 9 tic-tac-toe regions 
   See https://blogs.sas.com/content/iml/2023/12/04/bivariate-normal-probability-sas.html
*/
testName = "Test 5: 2-D Quadrant Test";
call randseed(123);
rho = 0.4;
R = (1 || rho) //
    (rho || 1);
a = -0.5; b =  0.5;  
c = -1.2; d =  0.8;
tol = 1e-3;
region = {"SW", "S", "SE", "W", "C", "E", "NW", "N", "NE"};

L = ( .M||.M ) //
    (  a||.M ) //
    (  b||.M ) //
    ( .M|| c ) //
    (  a|| c ) //
    (  b|| c ) //
    ( .M|| d ) //
    (  a|| d ) //
    (  b|| d );

U = (  a|| c ) //
    (  b|| c ) //
    ( .I|| c ) //
    (  a|| d ) //
    (  b|| d ) //
    ( .I|| d ) //
    (  a||.I ) //
    (  b||.I ) //
    ( .I||.I );                                        

Correct = j(9, 1, .);
Correct[1] = ProbBVN(.M,  a, .M,  c, rho);   /* SW */
Correct[2] = ProbBVN( a,  b, .M,  c, rho);   /* S  */
Correct[3] = ProbBVN( b, .P, .M,  c, rho);   /* SE */
Correct[4] = ProbBVN(.M,  a,  c,  d, rho);   /* W  */
Correct[5] = ProbBVN( a,  b,  c,  d, rho);   /* C  */
Correct[6] = ProbBVN( b, .P,  c,  d, rho);   /* E  */
Correct[7] = ProbBVN(.M,  a,  d, .P, rho);   /* NW */
Correct[8] = ProbBVN( a,  b,  d, .P, rho);   /* N  */
Correct[9] = ProbBVN( b, .P,  d, .P, rho);   /* NE */

Prob = j(nrow(region), 1, .);
do i = 1 to nrow(region);
   Prob[i] = probmvn_mod( L[i,], U[i,], R );
   *      print (region[i])[L='Region'] (L[i,])[L='L'] (U[i,])[L='U'] (Prob[i])[L='p'] (Correct[i])[L='Correct'];
   run check_test(testName, Prob[i], Correct[i], tol);
end;


quit;

