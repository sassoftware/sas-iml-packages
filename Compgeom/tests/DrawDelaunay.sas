proc iml;
%include "CG_define.sas";
%include "CG_proc.sas";
QUIT;

/* test calls to the CG_DrawConvexHull subroutine */
PROC IML;
load module=(CG_DrawDelaunay);

/**********************************************/
/* Test the DrawDelaunay upper-level routine  */
/**********************************************/
points={0  3.1, 0.5 2, 3.3 2.5, 0.5 1.7, 0 0, 0.45 -1, -2 1};
/* the Voronoi vertices should be as follows:
v = {
4    6    3,
6    5    7,
5    4    7,
4    5    6,
2    4    3,
1    2    3,
4    2    7,
2    1    7};
*/

title "Delaunay Triangulation";
title2 "Default Options";
run CG_DrawDelaunay(points);

title2 "Options={2,0,0}";
opt = {2,0,0};
run CG_DrawDelaunay(points, opt);

title2 "Options={0,1,1,1}";
opt = {0,1,1,1};
run CG_DrawDelaunay(points, opt);

title2 "Options={1,1,0.8,0}";
opt = {1,1,0.8,0};
run CG_DrawDelaunay(points, opt);
QUIT;
