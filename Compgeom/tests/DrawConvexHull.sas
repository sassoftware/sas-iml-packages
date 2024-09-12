proc iml;
%include "CG_define.sas";
%include "CG_proc.sas";
QUIT;

/* test calls to the CG_DrawConvexHull subroutine */
PROC IML;
load module=(CG_DrawConvexHull);

sites = {0  2, 0.5 2, 1 2, 0.5 1, 0 0, 0.5 0, 1  0,
         2 -1,   2 0, 2 1,   3 0, 4 1,   4 0, 4 -1,
         5  2,   5 1, 5 0,   6 0 };

title "Convex Hull with Default Options";
run CG_DrawConvexHull(sites);

/*  test 'short' options vector */
title "Convex Hull with Opt={0}";
opt = {0};
run CG_DrawConvexHull(sites, opt); 
QUIT;

PROC IML;
load module=(CG_DrawConvexHull);

sites = {0  2, 0.5 2, 1 2, 0.5 1, 0 0, 0.5 0, 1  0,
         2 -1,   2 0, 2 1,   3 0, 4 1,   4 0, 4 -1,
         5  2,   5 1, 5 0,   6 0 };

title "Convex Hull with Default Options";
run CG_DrawConvexHull(sites);
/*  test 'short' options vector */
title "Convex Hull with Opt={0}";
opt = {0};
run CG_DrawConvexHull(sites, opt);

title "Convex Hull with Opt={2, 1}";
opt = {2, 1};
run CG_DrawConvexHull(sites, opt);

/*  a Plackett-Burman design uses 9 experiments */
title "Convex Hull with Opt={0,0,0}";
opt = {0,0,0};
run CG_DrawConvexHull(sites, opt);

title "Convex Hull with Opt={1,1,1}";
opt = {1,1,1};
run CG_DrawConvexHull(sites, opt);

title "Convex Hull with Opt={2,2,0.5}";
opt = {2,2,0.5};
run CG_DrawConvexHull(sites, opt);

title "Convex Hull with Opt={0,0,1}";
opt = {0,0,1};
run CG_DrawConvexHull(sites, opt);

title "Convex Hull with Opt={1,1,0.5}";
opt = {1,1,0.5};
run CG_DrawConvexHull(sites, opt);

title "Convex Hull with Opt={2,2,0}";
opt = {2,2,0};
run CG_DrawConvexHull(sites, opt);

title "Convex Hull with Opt={0,0,0.5}";
opt = {0,0,0.5};
run CG_DrawConvexHull(sites, opt);

title "Convex Hull with Opt={1,1,0}";
opt = {1,1,0};
run CG_DrawConvexHull(sites, opt);

title "Convex Hull with Opt={2,2,1}";
opt = {2,2,1};
run CG_DrawConvexHull(sites, opt);

/* test on larger set of bivariate normal data */
call randseed(123);
mu = {0 0};
Sigma = {1 -0.8, -0.8 1};
x = randnormal(150, mu, Sigma);
title "Convex Hull for Bivariate Normal Data";
run CG_DrawConvexHull(x, {2, 1, ., 1}); /* display the PROC SGPLOT statements in the log */
QUIT;
