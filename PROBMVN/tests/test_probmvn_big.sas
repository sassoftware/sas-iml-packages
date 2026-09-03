options ps=max nodate nonumber;

proc iml;
reset wide;
load module=_all_;

/* validation test for high dimensions */
call randseed(12345);

/* Big 1: Test A large matrix.
   Use an NxN AR(1) correlation matrix, 
   which is a banded structure where the off-diagonal elements are 
   rho, rho^2, rho^3, etc.
   Each variables has limits of integration (-Infinity, 2).
   You might need to use -MEMSIZE 12G to run the dim=100 case.
*/
BaseTestName = "Big 1: AR(1) Corr. MC Validataion";
rho = 0.90;
NN = {32, 64, 100};
do i = 1 to nrow(NN);
   N = NN[i];
   TestName = cat(BaseTestName, " (dim=" + char(N,3) + ")");
   /* Limits: (-Infinity, 2) for all variables */
   L = j(1, N, .M);
   U = j(1, N, 2);
   /* Construct AR(1) Covariance Matrix */
   Sigma = CreateAR1(N, rho);
   prob = probmvn_mod(L, U, Sigma);

   /* Compute regular Monte Carlo estimate for comparison. Get a 95% confidence interval. */
   MC_list = MC_PROBMVN_CL(5E5, L, U, Sigma);
   MC_est  = MC_list$1;
   lower95 = MC_list$2;
   upper95 = MC_list$3;
   if prob < lower95 | prob > upper95 then do;
      run check_test(TestName, prob, MC_est);
      if max(abs(prob-MC_est)) > 1E-3 then 
         print prob MC_est lower95 upper95;
   end;
   else
      print (cat(TestName, " passes ---"));
end;


/* Big 2: Use orthant probabilities and an NxN AR(1) correlation matrix, 
   which is a banded structure where the off-diagonal elements are 
   rho, rho^2, rho^3, etc.
   Each variables has limits of integration [0, Infinity).
*/
TestName = "Big 2: AR(1) Corr, orthant prob";
rho = 0.90;
NN = {32, 64, 100};
correct = { 0.0293829, 0.0020252, 0.0001022};
do i = 1 to nrow(NN);
   N = NN[i];
   /* Orthant limits: [0, Infinity) for all variables */
   L = j(1, N, 0);
   U = j(1, N, .I);
   /* Construct AR(1) Covariance Matrix */
   Sigma = CreateAR1(N, rho);
   p = probmvn_mod(L, U, Sigma);
   run check_test(TestName + " (dim=" + char(N,3) + ")", p, correct[i], 0.0002);
end;

/* Big 3: Construct a dim=32 problem as the block diagonal of four dim=8 problems. 
   The permute the matrix so that the blocks are not contiguous.
   Use PROBMVN to solve each dim=8 problem. Compare the product of the four probabilities to the probability 
   for the 32x32 matrix.
*/
TestName = "Big 3: Block Corr with four 8x8 blocks, dim=32 ";
prob = j(4,1,.);
blockSize = 8;
/* Block 1: 8x8 AR(1) correlation with rho=0.5 */
rho1 = 0.5;
L1 = j(1, blockSize, .M);
U1 = j(1, blockSize, 2); 
R1 = CreateAR1(blockSize, rho1); /* AR(1) */
prob[1] = probmvn_mod(L1, U1, R1);
/* Block 2: 8x8 AR(1) correlation with negative rho=-0.6 */
rho2 = -0.6;
L2 = {-3 -2 -1 0 -1 -2 -3 0};
U2 = { 3  2  1 2  1  2  3 3};
R2 = CreateAR1(blockSize, rho2); /* AR(1) */ 
prob[2] = probmvn_mod(L2, U2, R2);
/* Block 3: Rank-1 correlation of the form vv' + D, where D is diagonal and D[i] = 1-v[i]**2 */
L3 = {-2 -1 0 -1 -1 -2 -3 0};
U3 = j(1, blockSize, .I); 
v = {-0.9 0.8 -0.7 0.6 -0.5 0.5 -0.6 0.7};
D = diag(1 - v##2);
R3 = v`*v + D;               /* PD correlation matrix */
prob[3] = probmvn_mod(L3, U3, R3);
/* Block 4: 8x8 AR(1) correlation from 
   https://blogs.sas.com/content/iml/2015/09/23/large-spd-matrix.html */
L4 = j(1, blockSize, -2);
U4 = j(1, blockSize,  2); 
h = 2/blockSize;           /* stepsize for decreasing sequence */
v = do(1, -1+h, -h);       /* {1, 0.75, 0.5, ..., -0.75} */
R4 = toeplitz(v);       /* PD correlation matrix */
prob[4] = probmvn_mod(L4, U4, R4);
correct = prod(prob);
/* Construct the block diagonal matrix */
R = block(R1, R2, R3, R4);
/* Permute the matrix so that the blocks are not contiguous */
perm = ranperm(4*blockSize);
R_perm = R[perm, perm];
L_perm = L1 || L2 || L3 || L4;
L_perm = L_perm[,perm];
U_perm = U1 || U2 || U3 || U4;   
U_perm = U_perm[,perm];
prob = probmvn_mod(L_perm, U_perm, R_perm);
run check_test(TestName, prob, correct);


