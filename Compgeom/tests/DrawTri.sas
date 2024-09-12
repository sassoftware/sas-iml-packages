proc iml;
%include "CG_define.sas";
%include "CG_proc.sas";
QUIT;

/* test calls to the CG_DrawConvexHull subroutine */
PROC IML;
load module=(CG_DrawTri);

/******************************************/
/* Test the DrawTri lower-level routine   */
/******************************************/
points={0  0, 3 0, 0 1, 1 1, 1 2};
v = {
1 2 3,
3 2 4,
3 4 5,
4 5 2};

title "DrawTri";
title2 "Default Options";
run CG_DrawTri(points, v);

opt = {2,0,0};
title2 "Options={2,0,0}";
run CG_DrawTri(points, v, opt);

opt = {0,1,1,1};
title2 "Options={0,1,1,1}";
run CG_DrawTri(points, v, opt);
QUIT;

/***********************************************************/
/* test the low-level helper function Compgeom_Make2DTri   */
/***********************************************************/

PROC IML;
load module=(Compgeom_Make2DTri CG_DrawTri);

/* Example of generating a regular grid */
x = do(0,4,0.5);
y = do(0,3,0.5);
L = Compgeom_Make2DTri(x,y);
verts = L$1;       * vertices;
tri   = L$2;       * indices ;

title "Triangulation from Unrotated 8x6 Grid (Nt=96)";
run CG_DrawTri(verts, tri);

title "Rotation by 5 Degrees";
L = Compgeom_Make2DTri(x, y, 5);
verts = L$1;      * vertices ;
tri   = L$2;      * indices ;
run CG_DrawTri(verts, tri);
QUIT;
