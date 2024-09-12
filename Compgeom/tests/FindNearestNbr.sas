proc iml;
%include "CG_define.sas";
%include "CG_proc.sas";
QUIT;

/***************************************/
/* test the function CG_FindNearestNbr */
/***************************************/

PROC IML;
load module=(CG_FindNearestNbr);

sites = {0  3.1, 0 2, 3.3 2.5, 1 1, 0.5 0.5, 0.45 -1, -2 1};
NN = CG_FindNearestNbr( sites );
print NN[c={'NNIndex' 'NNDist'}];

/* documentation example */
sites = {0  3.1, 0 2, 3.3 2.5, 1 1, 0.5 0.5, 0.45 -1, -2 1};
*run CG_DrawDelaunay(sites, {2,0,1}); /* label only the sites */

NN = CG_FindNearestNbr( sites );
Index = T(1:nrow(sites));
NNIndex = NN[,1]; 
NNDist = NN[,2];
print Index sites[c={'x' 'y'} L=""] NNIndex NNDist[F=5.3]; 

/*******************************************************************/
/* 2-D points: NN for many points uniformly at random in the plane */
/*******************************************************************/
call randseed(123);
N = 5000;
X = j(N, 2);
call randgen(X, "Uniform", 0, 10);

t0 = time();
NN = CG_FindNearestNbr(X);
time = time()-t0;
print time[L="Time to find nearest neighbors of 5,000 points in 2D"];

/*******************************************************************/
/* 3-D points: NN for many points uniformly at random in 3-D */
N = 5000;
X = j(N, 3);
call randgen(X, "Uniform", 0, 10);

t0 = time();
NN = CG_FindNearestNbr(X);
time = time()-t0;
print time[L="Time to find nearest neighbors of 5,000 points in 3D"];

N = T(do(5000, 15000, 2500));
times = j(nrow(N), 1);
do i = 1 to nrow(N);
   X = j(N[i], 2);
   call randgen(X, "Uniform", 0, 10);
   t0 = time();
   NN = CG_FindNearestNbr(X);
   times[i] = time()-t0;
end;

print N times;

ods graphics / width=400px height=300px;
title "Time to Find Nearest Neighbors for N Points";
call series(N, times) grid={x y};

QUIT;
