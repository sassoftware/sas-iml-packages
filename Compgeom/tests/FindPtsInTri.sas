proc iml;
%include "CG_define.sas";
%include "CG_proc.sas";
QUIT;


/***********************************************************/
/* test the function CG_FindPtsInTri */
/***********************************************************/
PROC IML;
load module=(Compgeom_Make2DTri CG_FindPtsInTri);

/* Example of generating a regular grid */
x = do(0,4,0.5);
y = do(0,3,0.5);
L = Compgeom_Make2DTri(x,y);
verts = L$1;       * vertices;
tri   = L$2;       * indices ;


/* now generate many points uniformly at random in the grid */
call randseed(123);
Nquery = 10000;
qmin = verts[><,];
qmax = verts[<>,];
xRange = qmin[1] || qmax[1];
yRange = qmin[2] || qmax[2];


/*
* OPTIONAL: enlarge range by 10% ;
dx = (qmax[1] - qmin[1]) / 10;
dy = (qmax[2] - qmin[2]) / 10;
xRange = (qmin[1] - dx) || (qmax[1] + dx);
yRange = (qmin[2] - dy) || (qmax[2] + dy);
*/

qx = randfun(Nquery, "uniform", xRange[1], xRange[2]);
qy = randfun(Nquery, "uniform", yRange[1], yRange[2]);
query = qx || qy;

t0 = time();
ID = CG_FindPtsInTri(verts, tri, query);
time = time()-t0;
print "Time to find 10,000 points = " time;

/* Tabulate the number of points in each triangle.
   Since the points are distributed uniformly, we expect a
   uniform distribution of triangle IDs.
   (In general, the expected counts in a triangle are proportional 
   to the area of the triangle.)
*/
title "Number of Points in Each Triangle";
title2 = "N = 10,000";
call tabulate(TriID, freq, ID, "nomiss");
mean = mean(colvec(freq));
refStmt = "refline " + char(mean,6.2) + " / axis=y;";
call scatter(TriID, Freq) other=refStmt;
QUIT;
