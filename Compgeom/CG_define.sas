*proc iml;

/**************************************************************
 * The 'private' helper functions start with the Compgeom_ prefix.
 * The 'public' functions start with the CG_ prefix are documented.
 *************************************************************/

/* Check the SYSVER macro to see if SAS 9.4 is running. 
   In SAS 9.4, the macro defines a module that calls the CVEXHULL function.
   In SAS Viya, the macro defines a module that calls the CONVEXHULL function.
*/
%macro DefineCHIndices;
%if %sysevalf(&sysver = 9.4) %then %do; /* SAS 9: call CVEXHULL */
start Compgeom_CHIndices(sites);
    indices = cvexhull(sites);
    CHIndices = indices[loc(indices>0)];
    return CHIndices;
finish;
%end;
%else %do;  /* SAS Viya: call CONVEXHULL */
start Compgeom_CHIndices(sites);
   call convexhull(indices, va, sites);
   CHIndices = indices[,1];
   return CHIndices;
finish;
%end;
%mend;
%DefineCHIndices;

/* Check the SYSVER macro to see if SAS 9.4 is running.
   In SAS Viya, the macro is empty and does nothing.
   In SAS 9.4, the macro defines a function that emulates the 
   COMPLETECASES function in SAS Viya.
*/
%macro DefineCompleteCases;
%IF %sysevalf(&sysver = 9.4) %THEN %DO;
/* Return the rows of a matrix M that are nonmissing. 
   Default: return a vector of the row numbers (indices) that are complete.
   If method="EXTRACT", return the data values for the complete rows.

   EXAMPLE:
   M = {1 2, 3 ., . 5, . ., 9 10};
   nonMissIdx = CompleteCases(M);
   nonMissData = CompleteCases(M, "Extract");
*/
start CompleteCases(M, method="ROW");
   idx = loc( countn(M, "row")=ncol(M) );
   if upcase(method)="ROW" then 
      return( idx );
   return( M[idx, ] );
finish;
store module=(CompleteCases);
%END;
%mend;
%DefineCompleteCases;


/*********************************************************/
/* The Polygon package functions */
/* Rick Wicklin (2016). The Polygon package. SAS Global Forum 2016 */
/* PolyPtInside: Determine whether a point is inside a polygon. This function
      supports non-simple (intersecting) polygons.
   Input: 
   P     an N x 2 matrix of (x,y) values for vertices. N > 2.
         Or an N x 3 matrix where the third column is an ID variable
         that identifies different polygons.
   pts   a k x 2 matrix of points to test 
   Output:
   S     a k x G matrix, where G is the number of unique ID values.
         S[i,j] = 0 if pts[i,] is outside the jth polygon
         S[i,j] = 1 if pts[i,] is strictly inside the jth polygon
         S[i,j] = . if pts[i,] is on the boundary of the jth polygon
*/

/* Find determinant at each vertex. R is a single point. 
   Let v_i = vector from R to P[i,].
   Let v_{i+1} = vector from R to P[i+1,].
   This function returns the determinant of the matrix with columns 
   [v_1  v_{i+1}].
   The i_th element is 0 when R = P[i,].
   */
start PolyWindingDet(P, R);
   lagIdx = 2:nrow(P) || 1;
   xi   = P[ ,1];       yi = P[ ,2];
   xip1 = xi[lagIdx]; yip1 = yi[lagIdx]; 
   det = (xi-R[1])#(yip1-R[2]) - (xip1-R[1])#(yi-R[2]);
   *print lagIdx, xi  yi xip1 yip1, R, det;
   return ( det );
finish;


/* Quadrant classification of vertices. R is a single point, 
   thought of as the origin. Return 0-3 when R is not a vertex.
   Return missing value for vertices that equal R. */
start PolyQuadrant(P, R);
   xi   = P[ ,1];       yi = P[ ,2];
   xPos = (xi > R[1]);  yPos = (yi > R[2]);
   xNeg = (xi < R[1]);  yNeg = (yi < R[2]);
   q = j(nrow(P), 1, .);                    /* origin is not in any quadrant */
   idx = loc( xPos & ^yNeg);  if ncol(idx)>0 then q[idx] = 0; /* angle [0, pi/2) */
   idx = loc(^xPos &  yPos);  if ncol(idx)>0 then q[idx] = 1; /* angle [pi/2, pi) */
   idx = loc( xNeg & ^yPos);  if ncol(idx)>0 then q[idx] = 2; /* angle [pi, 3pi/2) */
   idx = loc(^xNeg &  yNeg);  if ncol(idx)>0 then q[idx] = 3; /* angle [3pi/2, 2pi) */
   return ( q );
finish;

/* Which points are strictly within the bounding box around P
   [min(P[,1]), max(P[,1])] x [min(P[,2]), max(P[,2])] */
start PolyPtInBoundingBox(P, pts);
   in = (min(P[,1]) <= pts[,1]) & (pts[,1] <= max(P[,1])) &
        (min(P[,2]) <= pts[,2]) & (pts[,2] <= max(P[,2]));
   return ( in );
finish;

/* winding number w(P, R) = sum( deltaBar ) as given on p. 135-136 of
   Horman and Agathos (2001), "The point in polygon problem for arbitrary
   polygons," _Computational Geometry_, vol 20, 131-144.
   R is a single point.
*/
start Poly1WindingNumber(P, R);
   det = PolyWindingDet(P, R);
   q = PolyQuadrant(P, R);
   if any(q=.) then return (.); /* R is a vertex of P */
   N = nrow(q);
   qDif = dif(q);          /* q[i+1] - q[i] */
   qDif = qDif[2:N] // (q[1] - q[N]);  /* all we care about is sum(w), so put last diff into element that has missing value */
   /* Algorithm 2: p. 136 */
   idx = loc(((qDif =  2) | (qDif = -2)) & (det = 0));  if ncol(idx)>0 then return(.);  /* R is on edge of P */
   w = j(nrow(q), 1, 0);
   idx = loc(qDif = -3);                if ncol(idx)>0 then w[idx] =  1;
   idx = loc(qDif =  3);                if ncol(idx)>0 then w[idx] = -1;
   idx = loc((qDif = -2) & (det > 0));  if ncol(idx)>0 then w[idx] =  1;
   idx = loc((qDif =  2) & (det < 0));  if ncol(idx)>0 then w[idx] = -1;
   *print R, (P-R)[c={Cx Cy}] det q qDif w;
   return( sum(w) );
finish;

/* Given Nx2 matrix of vertices, return a vector that spcifies the 
   winding number of the polygon with respect to pts[i,] */
start PolyWindingNumber(P, pts);
   /* any pts that are outside of the bounding box have winding number 0 */
   w = PolyPtInBoundingBox(P, pts);
   idx = loc(w = 1);
   if ncol(idx)=0 then return( w );
   do i = 1 to ncol(idx);
      j = idx[i];
      w[j] = Poly1WindingNumber(P, pts[j,]);
   end;   
   return(w);
finish;

/* Given Nx2 matrix of vertices, return vector that spcifies whether pts[i,] 
   is inside the polygon (1), outside the polygon (0), or on the boundary (.) */
start _PolyPtInside(P, pts);
   wn = PolyWindingNumber(P, pts);
   /* pt in polygon if winding number is odd */
   /* if winding number is negative, mod(wn,2) is in {0 -1}. Take ABS value. */
   in = choose(wn=., ., abs(mod(wn,2)) = 1); /* winding number is odd */
   return(in);
finish;

/* PolyPtInside: Determine whether a point is inside a polygon. This function
   supports non-simple (intersecting) polygons.
*/
start PolyPtInside(P, pts);
   if ncol(P)=2 then
      return( _PolyPtInside(P, pts) );

   ID = P[,3];
   u = uniqueby(ID);         /* starting index for each group */
   result = j(nrow(pts), nrow(u)); /* allocate vector to hold results */
   u = u // (nrow(ID)+1);    /* append (N+1) to end of indices */
   do i = 1 to nrow(u)-1;    /* for each group... */
      idx = u[i]:(u[i+1]-1); /* get rows in group */
      result[,i] = _PolyPtInside(P[idx, 1:2], pts);
   end;
   return( result );
finish;

store module = (
PolyWindingDet
PolyQuadrant
PolyPtInBoundingBox
Poly1WindingNumber
PolyWindingNumber
_PolyPtInside
PolyPtInside
);
/*****************************************************************/

/* Compute the centroid of a triangle, which is the mean of its vertices */
start _TriCentroid(T);  /* T is (3 x 2) matrix */
   return( mean(T) );
finish;

/* Compute the centroid for a set of triangles.
   P is a (3*k x 3) matrix that specifies k triangles
   P[,1:2] is the set of (x,y) coordinates and 
   P[,3] is the ID variable, such as {1,1,1, 2,2,2, 3,3,3,...}

   Return (k x 3) matrix C where C[,1:2] are (x,y) coords of centroids
   and C[,3] contains the ID values.
 */
start Compgeom_TriCentroid(P);
   if ncol(P)=2 then  return( _TriCentroid(P) );
   ID = P[,3];
   u = uniqueby(ID);         /* starting index for each group */
   k = nrow(u);
   result = j(k, 2);   /* allocate vector to hold results */
   u = u // (nrow(ID)+1);    /* append (N+1) to end of indices */
   do i = 1 to k;    /* for each group... */
      idx = u[i]:(u[i+1]-1); /* get rows in group */
      result[i,] = _TriCentroid( P[idx, 1:2] );
   end;
   return( result || ID[u[1:k]] );
finish;

/*************************************************************/
/* UTILITY FUNCTIONS FOR LOCATION PROBLEM */
/* For each row of a data matrix, X, return whether row is in BBox.
   X : N x d  data matrix whose rows are the points to test
   min : 1 x d vector that represents the lower-left corner of a bounding box
   max : 1 x d vector that represents the upper-right corner of a bounding box
   RETURN : b, an N x 1 vector such that b[i]=1 iff X[i,] is in the bounding box.
   
   Ex:
   X = {1 2 3,  4 5 6};
   min = {0 0 0};
   max = {4 4 10};
   b = Compgeom_InBBox(X, min, max);
*/ 
start Compgeom_InBBox(X, min, max);
   d = ncol(X);
   bMin = ( (min <= X)[,+] = d );  /* bMin[i]=1 iff min  <=  X[i,] */ 
   bMax = ( (X <= max)[,+] = d );  /* bMax[i]=1 iff X[i, ]>= max   */ 
   return (bMin & bMax);           /* = 1 iff row is in BBox */
finish;

/* Compute barycentric coordinates for x with respect to the triangle 
   whose vertices are the three rows of P. See
   https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Vertex_approach
   X: an N x 2 data matrix. Each row is Euclidean coords of points in the plane
   P: a 3 x 2 matrix. Each rows is a vertex of a planar triangle.
   RETURN : Barycentric coords of X w/r/t P.
*/
start Compgeom_Barycentric(x, P);
   A = P` // j(1, nrow(P), 1);       /* add constraint coefficients for sum(barycentric) */
   b = x` // j(1,nrow(x),1);         /* add constraint: sum(barycentric) = 1 */
   v = solve(A, b);
   return( v` );                     /* return a row vetor */
finish;

/* For each input point, return 1/0 if the point, x, is inside/outside the 
   triangle whose vertices are the rows of P. */
start Compgeom_PtInTriangle(x, P);
   c = Compgeom_Barycentric(x, P);    /* barycentric coordinates */
   return( c[, ><] >=0 );             /* are all coords nonnegative? */
finish;

/* For testing: Create a triangulation on a regular (x,y) grid.
   The grid is created as ExpandGrid(_X, _Y),
   _X : 1 x Nx vector that contains ordered set of points for X grid.
        For example, you could use 
        X = do(0, 5, 0.5);    * evenly spaced ;
        or
        X = {0 1 1.5 2 4 5};  * ireegular spacing ;
   _Y : 1 x Ny vector that contains ordered set of points for Y grid.
        Same examples as for X, although X and Y can be different.
   rotDeg : a number (usually in [0, 90]) that specifies the degree 
        of counterclockwise rotation (about the origin) to apply to the grid.
        
   RETURN a list, L, such that 
   L$1 is a (Nx*Ny) x 2 set of vertices, and 
   L$2 is an Nt x 3 set of indices such that each row specifies the vertices of a triangle.
*/
start Compgeom_Make2DTri(_x,_y, rotDeg=0);
   x = rowvec(_x); y = rowvec(_y);
   Nx = ncol(x);
   Ny = ncol(y);  
   xy = expandgrid(y,x);
   xy = xy[, {2 1}];
   if rotDeg ^=0 then do;
      theta = rotDeg / 180 * constant('pi');
      R = ( cos(theta) || -sin(theta) ) //
          ( sin(theta) ||  cos(theta) );
      xy = xy*R`;
   end;
   *print xy[r=(1:(Nx*Ny))];
   
   dim = Ny || Nx;
   T = j((Nx-1)*(Ny-1), 3);
   S = j((Nx-1)*(Ny-1), 3);
   r = 1;
   do i = 1 to Ny-1;
      do j = 1 to Nx-1;
          /* T triangles look like |\
                                   |_\    'base' on bottom   */
         T[r,1] = sub2ndx(dim,     i||j    );
         T[r,2] = sub2ndx(dim,     i||(j+1));
         T[r,3] = sub2ndx(dim,  (i+1)||j   );
          /* S triangles look like ___
                                   \ |
                                    \|    'base' on top  */
         S[r,1] = sub2ndx(dim,     i||(j+1));
         S[r,2] = sub2ndx(dim,  (i+1)||(j+1));
         S[r,3] = sub2ndx(dim,  (i+1)||j   );
         r = r+1;
      end;
   end;
   *print S, T;
   tri = T // S;
   return( [xy, tri] );
finish;
              
/*************************************************************/
/* CG_FindPtsInTri: Solve the "location problem" for a triangulation.
   Given:
   vert = Nx x 2 matrix set of vertices, X, for a triangulation
   triIdx = Nt x 3 matrix of indices, T.
       Each row determines a triangle. The vertices of the i_th triangle are the points
       X[ T[i,1], ]
       X[ T[i,2], ]
       X[ T[i,3], ]
   query = Nq x 2 matrix of query points, Q. For the i_th row, we want to find the 
       triangle T[m,] that contains Q[i,]

   This routine is written to be efficient when there are many query points
   and a small or moderate number of triangles. The algorithm is:
   1. Initialize list of query points that are yet to be classified.
   2. For each triangle:
   3.    Compute bounding box 
   4.    From remaining query points, find those in BBox
   5.    For that subset, find if any query pts are in triangle
   6.    Update list of query points that are yet to be classified
*/
start CG_FindPtsInTri(vert, triIdx, query);
   LocOfPt = j(nrow(query), 1, .);
   /* 1. Initialize list of query points that are yet to be classified. */
   qIdx = 1:nrow(query);
   /* 2. For each triangle: */
   do i = 1 to nrow(triIdx) while(ncol(qIdx)>0);
      /* 3. Compute bounding box */
      tri = vert[ triIdx[i,], ];
      min = tri[ ><, ];  /* (minX, minY) is lower left corner of BBox */
      max = tri[ <>, ];  /* (maxX, maxY) is upper right corner of BBox */
   
      /* 4. From remaining query points, find those in BBox */
      Q = query[qIdx,];
      BBIdx = loc(Compgeom_inBBox(Q, min, max)); /* in BBox */
      if ncol(BBIdx) > 0 then do;
         testIdx = qIdx[,BBIdx];        /* indices to test */
         Q = query[testIdx,];
         /* 5. For that subset, find if any query pts are in triangle */
         TIdx = loc( Compgeom_PtInTriangle(Q, tri) );
         if ncol(TIdx)>0 then do;
            foundIdx = testIdx[,TIdx];  /* in the triangle */
            *print i tri, i foundIdx;
            LocOfPt[foundIdx] = i;      /* save the location */
            /* 6. Update list of query points that are yet to be classified */
            qIdx = setdif(qIdx, foundIdx);
         end;
      end;
   end;
   return LocOfPt;
finish;

/*************************************************************/
/* CG_FindNearestNbr: Find the nearest neighbor to each site
   Given:
   sites = N x d matrix set of points in d-dimensional Euclidean space

   Return:
   An N x 2 matrix, NN, where
   NN[,1] is a vector of indices with values 1-N. The value NN[i,1] is the 
          index for the nearest neighbor to sites[i,]. 
   NN[,2] is a vector of distances. The value NN[i,2] is the distance from 
          sites[i,] to its nearest neighbor, sites[NN[i,2], ].
*/
start CG_FindNearestNbr( sites );
   return Compgeom_NNDistance( sites );
finish;

/* A brute-force computation for the nearest neighbor:
   Compute all N**2 distances, then find smallest distance
   in each row. Return the Index (column) of minimum and the distance itself.
*/   
start Compgeom_NNDistance(X);
   D = distance(X);              /* N x N distance matrix */
   nr = nrow(D);
   diagIdx = do(1, nr*nr, nr+1); /* index diagonal elements */
   D[diagIdx] = .;               /* set diagonal elements to missing */
 
   dist = D[ ,><];               /* smallest distance in each row */
   nbrIdx = D[ ,>:<];            /* column of smallest distance in each row */
   /* coords of closest neighbor is X[nbrIdx, ] */
   return nbrIdx || dist;
finish;


/*****************************************************/
/* FUNCTIONS TO GRAPH A BOUNDED VORONOI TESSELLATION */
/*****************************************************/

/* Functions and helper function for plotting a 2-D Voronoi diagram. 
   The diagram will overlay the following objects:
   1. A collection of sites
   2. (Optional) The convex hull of the sites
   3. The vertices of the Voronoi diagram
   4. The edges of the Voronoi diagram.

   Because the Voronoi diagram contains unbounded polygons, the 
   visualization intersects the Voronoi cells with a bounding box 
   that contains all sites and Voronoi vertices.
*/

/* TODO:
   1. Do we need SIDE?
   2. Do we need to return the Angles?
   3. Handle the degenerate case where there are three sites. Then there 
      is one Voronoi vertex and the cells are {1 2, 1 2, 1 2}.
   4. Add option to Vor_InflateBBox to round the inflated BBox to a convenient
      unit such as an integer, 0.5, or 0.1
   5. Handle degenerate cases where points are not in general position:
      5A) A site can be the CH but not a vertex. This can occur when three or more
      sites are collinear. One option: detect the case 
      and manually add the point as a new "vertex" in the convex hull.
      5B) Handle symmetric configurations such as sites on a square.
      This results in Voronoi vertices that have more than 
      three edges emerging from them. One idea: detect this 
      situation and sightly perturb the sites to destroy the symmetry.
*/

/*-----------------------------------*/
/* Helper functions for bounding box */
/* Use Vor_ prefix for the helper functions for the Voronoi visualization */
/*-----------------------------------*/

/* Check whether all 2-D points are strictly within the bounding box 
   defined by BBox = {minX  minY  maxX  maxY} 

   EXAMPLE:
   BBox = {0 1, 2 3};
   Pts = {1 1, 1 2, 2 2, 1.5 2.5};
   allIn = Vor_PtsInBBox(Pts, BBox);   *returns TRUE;
*/
start Vor_PtsInBBox(_xy, BBox);
   xy = CompleteCases(_xy, "extract");
   x = xy[,1]; y = xy[,2];
   allIn = all( (BBox[1] <= x      ) & 
                (x       <= BBox[3]) &
                (BBox[2] <= y) & 
                (y       <= BBox[4]) );
   return ( allIn );
finish;


/* BBox = {minX  minY  maxX  maxY}
   The rectangle has the corners
   (bb[1], bb[2]),  (bb[3], bb[2]),  (bb[3], bb[4]),  (bb[1], bb[4])
   |__________ Side 1 ____________|
                    |__________ Side 2 ____________|
                                     |__________ Side 3 ____________|
   Given a ray, v, based at a point, p, inside a BBox,
   return the point on the rectangle at which the ray 
   intersects the BBox.
  
   The ray is the parametric form l(t) = p + t*v.
   Intersect VERT line x=X for t* s.t. 
         l(t*)_x = p_x + t* * v_x
     ==> t* = (X - p_x) / v_x
   1) t* must be > 0 if it intesects line at X.
   2) If t* > 0, check that yMin <= l(t*)_y <= yMax

   Similarly, intersect HORZ line y=Y for t* s.t. 
         l(t*)_y = p_y + t* * v_y
     ==> t* = (Y - p_y) / v_y
   1) t* must be > 0 
   2) If t* > 0, check that xMin <= l(t*)_x <= xMax

   EXAMPLE:
   BBox = {-1 -1 1 1};          *square with side length 2;
   origin = {0 0};              *base of vector;
   v = {1 0, 1 1, 1 0, -1 1, -1 0, -1 -1, 0 -1, 1 -1};  *the intersection pts;
   dir = v / sqrt(v[,##]);      *unit vectors;
   xsectPts = j(nrow(v), 2, .); *we should end up with the intersection pts;
   do i = 1 to nrow(v);
      xsectPts[i,] = Vor_RayIntersectRect(dir[i,], origin, BBox); 
   end;
   print v xsectPts;
 
*/
start Vor_RayIntersectRect(v, p, BBox);
   if ^Vor_PtsInBBox(p, BBox) then do;
       print "ERROR: Something is wrong! The point p is not in BBox";
       print p, BBox;
       return( {. .} );
   end;
   xMin = BBox[1]; yMin = BBox[2]; xMax = BBox[3]; yMax = BBox[4];
   /* does ray based at p intersect Side 1 (bottom) or Side 3 (top)? */
   if v[2] ^= 0 then do;
      Yval = yMin || yMax;
      do i = 1 to 2;
         Y = Yval[i];
         tStar = (Y - p[2]) / v[2];
         if tStar < 0 then
            continue;
         z = p + tStar * v;
         *print i Y tStar z;
         if xMin <= z[1] & z[1] <= xMax then 
            return( z );
      end;
   end;
   /* does ray based at p intersect Side 4 (left) of Side 2 (right)? */
   if v[1] ^= 0 then do;
      Xval = xMin || xMax;
      do i = 1 to 2;
         X = Xval[i];
         tStar = (X - p[1]) / v[1];
         if tStar < 0 then
            continue;
         z = p + tStar * v;
         *print i X tStar z;
         if yMin <= z[2] & z[2] <= yMax then 
            return( z );
      end;
   end;
   return( {. .} );   /* something is wrong; should never get here if p in BBox */
finish;

/* Given a point, P, on a rectangle, Rect, return a number that 
   indicates which side of the rectangle P is on:
   Side numbers: Bottom=1, Right=2, Top=3, Left=4
   Corners: SW=1, SE=2, NE=3, NW=4   

   NOT NEEDED. USE FOR DEBUGGING AND FOR MANUAL VERIFICATION.

   EXAMPLE:
   BBox = {-1 -1 1 1};          *square with side length 2;
   v = {1 0, 1 1, 1 0, -1 1, -1 0, -1 -1, 0 -1, 1 -1};  *the intersection pts;
   side = j(nrow(v), 1, .); 
   do i = 1 to nrow(v);
      side[i,] = Vor_WhichSideRect(v[i,], BBox); 
   end;
   print v side;
*/
start Vor_WhichSideRect(P, Rect, delta=1e-4);
   xMin = Rect[1]; yMin = Rect[2]; xMax = Rect[3]; yMax = Rect[4];
   vertEps = (yMax - yMin)*delta;
   horzEps = (xMax - xMin)*delta;
   onBottom = ( abs(P[2] - yMin) < vertEps );
   onTop    = ( abs(P[2] - yMax) < vertEps );
   onLeft   = ( abs(P[1] - xMin) < horzEps );
   onRight  = ( abs(P[1] - xMax) < horzEps );
   /* negative values for corners, enumerated SW, SE, NE, NW */
   if onBottom & onLeft  then return( -1 );
   if onBottom & onRight then return( -2 );
   if onTop    & onRight then return( -3 );
   if onTop    & onLeft  then return( -4 );
   /* side numbers: Bottom=1, Right=2, Top=3, Left=4 */
   if onBottom then return( 1 );
   if onRight  then return( 2 );
   if onTop    then return( 3 );
   if onLeft   then return( 4 );
   return( . ); /* error or delta too small */
finish;

/* Compute the bounding box for a set of 2-D points. Then 
   inflate the BBox by 1+Factor in all directions.
   By default, Factor = 10% = 0.1

   EXAMPLE;
   Pts = {0 0, 1 0, 0 1};
   BBox = Vor_InflateBBox(Pts);
   print BBox;
*/
start Vor_InflateBBox(Pts, Factor = 0.1); 
   if Factor <= 0 then do;
      print "Invalid Factor in Vor_InflateBBox function";
      return( {. . . .} );
   end;
   multFactor = 1 + 2*Factor;
   xMin = min(Pts[,1]);
   xMax = max(Pts[,1]);
   xC = (xMin + xMax) / 2;
   xWidth = multFactor * (xC - xMin);
   BxMin = xC - xWidth;
   BxMax = xC + xWidth;

   yMin = min(Pts[,2]);
   yMax = max(Pts[,2]);
   yC = (yMin + yMax) / 2;
   yWidth = multFactor * (yC - yMin);
   ByMin = yC - yWidth;
   ByMax = yC + yWidth;

   BBox = BxMin || ByMin || BxMax || ByMax;
   /* TODO: round BBox to a convenient unit. For now, round to integer away from 0 */
   BBox = int(BBox) + sign(BBox);
   return BBox;
finish;

/*-----------------------------------*/
/* Other helper functions            */
/*-----------------------------------*/

/* X is k x 2 matrix of 2-D points; Y is N x 2 matrix.
   Neither X nor Y should have missing values.
   Return row numbers of X that are in common with Y (up to delta)
   That is, return all i such that X[i,] is a row of Y.
   Return an empty matrix if X and Y are disjoint.
   Example:
   X = {1 1, 2 2, 3 0, 4 0};
   Y = {3 4, -1 2, -1 1, 2 2, 4 0, 5 3};
   k = Vor_PtsInCommon2D(X, Y);
   print k;  * k=3 bc X[3,] is a row of Y;
*/
start Vor_PtsInCommon2D(X, Y, delta=1E-4);
   commonIdx = {};
   do i = 1 to nrow(X);
      diff = abs(X[i,] - Y)[,<>];  /* max of abs diff for each row */
      if any( diff < delta ) then 
         commonIdx = commonIdx || i;
   end;
   return commonIdx;
finish;

/* dynamic reallocation of a matrix when you are appending rows.
   If you want to append more rows than the matrix contains,
   call M=DynamicRealloc(M), which adds rows to M then copies over the
   previous data.

   EXAMPLE:
   M = {1 2 3, 4 5 6};
   M = _DynamRealloc(M);
   print M;
*/
start _DynamRealloc(M, factor=2);
    mult = max(factor, 1);
    newM = j(int(mult*nrow(M)), ncol(M), .);
    newM[ 1:nrow(M), ] = M;
    return newM;
finish;

/*-------------------------------------------------------*/
/* Intersect unbounded Voronoi cells with a bounding box */
/*-------------------------------------------------------*/

/* NOTES: Find and mark special points on the BBox, including corners. 
These will be used to create the bounded polygons.
Use the fact that there is one Voronoi cell for each site.
The i_th Voronoi cell contains all points that are 
closer to site i than to any other site.

Construct a reference table (BdyPts). Find points(x,y) on the BBox that are 
equidistant between two sites (on the CH)
        BBox
   Cell Side x y Theta
   1    3    - -   -  
   1    3    - -   -  
   2    3    - -   -    <= repeat of prev coords
   2   -3    - -   -    <= corner
   2   -3    - -   -    <= might have two corners
   2    1    - -   -  
   3    1    - -   -    <= repeat of prev coords
   3    1    - -   -  

Generically, we expect a corner to be close to a unique vertex, 
but rarely it could be on the boundary between two cells.
We append the corners to the table. We can use the ATAN2 function 
to sort the points by angle from the centroid. See
https://blogs.sas.com/content/iml/2021/11/17/order-vertices-convex-polygon.html
*/

/* Find points on the BBox that represent the intersection of the infinite Voronoi
   cells with the BBox.  Return matrix BdyPts that has 5 columns:
   Col1: Cell = the unbounded Voronoi cell that we intersect with the BBox
   Col2: side = which side of the box? See Vor_WhichSideRect function for details.
   Col3: x = horiz coordinate of intersection point on BBox
   Col4: y = vert coordinate of intersection point on BBox
   Col5: Angle = angle in [0, 2*pi] that vector from center of BBox to the 
         intersection point makes RELATIVE to the first intersection point, 
         which is assigned Angle=0. This column is used to sort the intersection 
         points on the BBox in counterclockwise order.Probably not needed as output.

   Algorithm: Loop over the sites on the CH of the sites. Each site is a vertex of the CH.
   Find the midpoint of the line segment of the CH that comes before/after
   the site. Shoot a ray from the perpendicular bisector until it hits the
   BBox. Store the intersection point, what BBox side it is on,
   and the number of the site (=number of Voronoi vertex, too).
   Add in the corners of the BBox.

   INPUT:
   CHS : an N x 3 matrix. Each row is a point on the convex hull of a set of sites.
         Col1 = x coord; 
         Col2 = y coord
         Col3 = numeric ID var 1,2,...,N
   BBox : 1 x 4 matrix {xMin yMin xMax yMax} for bounding rectangle that contains CH

   EXAMPLE: 
   CHS = {2.0 3.4  1,
          1.2 2.3  2,
          3.5 2.5  3,
          1.3 1.5  4};
   BBox = {1 1  4 4};
   BdyPts = Vor_CHRaysIntersectRect(CHS, BBox);
   print BdyPts[c={'Cell' 'side' 'x' 'y' 'Angle'}];;
*/
start Vor_CHRaysIntersectRect(CHS, BBox);
   CHSidx = CHS[,3]; CHSxy = CHS[,1:2];
   nCHEdges = nrow(CHS);
   xMin = BBox[1]; yMin = BBox[2]; xMax = BBox[3]; yMax = BBox[4];
   BdyPts = j(2*(nCHEdges + 4), 5, .);  /* two for each CH edge and corner */
   tblRow = 1;
   C = mean(CHSxy); /* the center of a convex polygon is INSIDE the polygon */
   L = CHSxy[nCHEdges, ];
   do i = 1 to nCHEdges;
      R = CHSxy[i,];
      P = (L + R)/2;          /* midpoint of edge */
      w = R - L;              /* vector along edge */
      wPerp = -w[2] || w[1];  /* a vector perpendicular to the edge */
      /* check to make sure we have an OUTWARD pointing normal vector */
      v = P - C;              /* points outward from CH */
      if wPerp*v` < 0 then        /* oops! wPerp points into CH */
         wPerp = -wPerp;          /* wPerp is now outward pointing vector */
      /* where does wPerp intersect the sides of the rectangle? */
      s = Vor_RayIntersectRect(wPerp, P, BBox);
      side = Vor_WhichSideRect(s, BBox);
      *print i (CHSidx[i,]) P wPerp s side;
      /* add this point to the table TWICE: once for each end of segment */
      BdyPts[tblRow, 1:4] = CHSidx[i] || side || s;
      tblRow = tblRow + 1;  /* incr counter */
      if i=1 then 
         BdyPts[tblRow, 1:4] = CHSidx[nCHEdges] || side || s;
      else 
         BdyPts[tblRow, 1:4] = CHSidx[i-1] || side || s;
      tblRow = tblRow + 1;  /* incr counter */
      L = R;
   end;
   *print BdyPts[c={'Cell' 'side' 'x' 'y' 'Angle'}];

   /* Next, compute the cells for each corner of the BBox */
   Rect = (xMin || yMin) //
          (xMax || yMin) //
          (xMax || yMax) //
          (xMin || yMax);
   /* if a corner is in the boundary tale already, we do not need to 
      add it */
   dupRows = Vor_PtsInCommon2D( Rect, BdyPts[1:tblRow-1,3:4] );
   keepRows = setdif(1:nrow(Rect), dupRows);

   if ncol(keepRows) > 0 then do;
      Rect = Rect[keepRows,];
      D = distance(Rect, CHSxy);
      *print D[c=CHSidx];
      closestIdx = D[ , >:<];             /* columns of minimum */
      vertNum = CHSidx[closestIdx, 1];    /* index of minimum */
      *print closestIdx;
      side = j(nrow(Rect), 1, .);
      do i = 1 to nrow(Rect);
         side[i] = Vor_WhichSideRect(Rect[i,], BBox);
      end;
      /* add this information to the table */
      BdyPts[tblRow:tblRow+nrow(Rect)-1, 1:4] = vertNum || side || Rect;
   end;

   /* remove any unused rows */
   keepIdx = loc( countmiss(BdyPts,"row") < ncol(BdyPts) );
   BdyPts = BdyPts[keepidx,];
   *print BdyPts ;

   /* Assign angles to the points on the BBox.
      Use the ATAN2 function to order the points on the rectangle
      in counterclockwise order, which is the same order as the 
      CH of the sites. Compute angles from the center of the bounding rectangle. */
   pi = constant("pi");
   CBBox = (xMin+xMax)/2 || (yMin+yMax)/2;
   Angle = atan2(BdyPts[,4]-CBBox[2], BdyPts[,3]-CBBox[1]); /* ATAN2(y,x): Notice order! */
   /* ATAN2 outputs values in [-pi, pi). */
   /* For continuity, adjust the angles into [0, 2*pi) */
   idx = loc(Angle<0);
   if ncol(idx)>0 then 
      Angle[idx] = Angle[idx] + 2*pi;
   /* Make the first angle (T0) the reference angle. */
   Angle = Angle - Angle[1];   /* in [0, 2*pi-T0) U [-T0, 0) */
   /* For continuity, AGAIN adjust the angles into [0, 2*pi) */
   idx = loc(Angle<0);
   if ncol(idx)>0 then 
      Angle[idx] = Angle[idx] + 2*pi;
   /* The first angle has two entries. Use T0=0 for the first and T0=2*pi for the last */
   Angle[2] = 2*pi;
   /* finally, assign angles to 5th column */

   BdyPts[,5] = Angle;
   call sort(BdyPts, {1 5}); /* sort by cell, then by angle */
   *print BdyPts[c={'Cell' 'side' 'x' 'y' 'Angle'}];
   return BdyPts;
finish;

/* assume Poly is sorted by Poly[,3], which is the ID.
   For each Polygon, sort the vertices counterclockwise according to the 
   angle made by the line segment between the centroid and the vertex.
   The polygons are updated in place, so Poly is overwritten. 

   EXAMPLE:
   Poly = {0 0, 0 1, 1 1, 1 0};
   run Vor_SortPolyByAngle(Poly);            
   print Poly;
   *can also sort multiple polygons;
   Poly1 = {0 0, 0 1, 1 1, 1 0};
   Poly2 = {1 0, 1 1, 2 1, 3 0};
   Poly = (Poly1 || j(nrow(Poly1),1,1)) //
          (Poly2 || j(nrow(Poly2),1,1));
   run Vor_SortPolyByAngle(Poly);            
   print Poly;
*/
start Vor_SortPolyByAngle(Poly) global(g_debugBVP);
   if ncol(g_debugBVP)=0 then 
      g_debugBVP = 0;
   DEBUG = g_debugBVP;
   pi = constant("pi");
   /* 2. Obtain row numbers for the first observation in each level. */
   if ncol(Poly)=2 then 
      b = 1 // (nrow(Poly)+1);
   else if ncol(Poly)=3 then do;
      b = uniqueby(Poly, 3);     /* b[i] = beginning of i_th category */
      u = Poly[b, 3];            /* get unique values (if needed) */
      b = b // (nrow(Poly)+1);   /* trick: append (n+1) to end of b */
   end;
   do i = 1 to nrow(b)-1;        /* 4. For each level... */
      idx = b[i]:(b[i+1]-1);     /* 5. Find observations in level */
      xy = Poly[idx, 1:2];       /* extract vertices for this polygon */
      if DEBUG then PRINT i, xy[L='Orig Order' r=(1:nrow(xy))];
      seg = xy - mean(xy);
      Angle = atan2(seg[,2], seg[,1]); /* ATAN2(y,x) in [-pi,pi]: Notice order! */
      /* For continuity, adjust the angles into [0, 2*pi) */
      jdx = loc(Angle<0);
      if ncol(jdx)>0 then 
         Angle[jdx] = Angle[jdx] + 2*pi;
      call sortndx(ndx, Angle);
      sortedxy = xy[ndx,];
      if DEBUG then PRINT sortedxy[L='Sorted Order' r=(1:nrow(xy))];
      Poly[idx, 1:2] = sortedxy;  /* vertices are now in counterclockwise order */
   end;
finish;

/* In the call
     call voronoi(V, vertIdx, sites);
   the output V is a kx2 where each row is a Voronoi vertex:
   V = {. .,   <= point at infinity
        x2 y2,
        ...
        xk, yk};
   The output vertIdx is a set of indices (padded by missing values)
   that describe the Voronoi cells. A row of vertIdx that contains
   1 is an unbounded cell. A row that does not contain 1 is bounded.
   For example:
   {2 3 4}  is a bounded triangular cell
   {1 3 4} is an unbounded cell that has a line segment between vertices 
           3 & 4 and rays from vertices 3 and 4 that are perpendicular to 
           an edge of the CH of the sites.
   {1 2}   is a special cell that has one vertex from which one or more rays emerge
*/

/* Assume the vertIdx vector does NOT contain a 1.
   Return a bounded voronoi polygon 
   EXAMPLE:
   V = {. ., 1.5 -0.5, -1.5 2.5, 3.1666667 1.1666667, 0.5 2.5};
   idx = {5 3 2 4};
   Poly = Vor_ExtractBoundedPoly(idx, V);
   print Poly;
*/
start Vor_ExtractBoundedPoly(_vertIdx, V);
   vertIdx = _vertIdx[ ,loc(_vertIdx > 1) ];  /* get rid of missing values */
   return( V[vertIdx, ] );
finish;

/* Assume the vertIdx vector contains a 1, which means it is an 
   unbounded polygon. Return a bounded Voronoi polygon by interceting
   the cell with a bounding box (BBox) that contains all Voronoi vertices 

   EXAMPLE:
   V = {. ., 1.5 -0.5, -1.5 2.5, 3.1666667 1.1666667, 0.5 2.5};
   idx = {1 2 4 .};
   BBox = {-2 -2  4 5};
   CHS = {0  0 1,
          3 -1 2,
          2  3 3,
          1  4 4}; 
   AllBdyPts = Vor_CHRaysIntersectRect(CHS, BBox); 
   *idx corresponds to 2nd VCell, which is rows 4:6;
   BdyPtsxy = AllBdyPts[4:6,{3 4}];     *extract (x,y) coordinates;
   Poly = Vor_ExtractUnboundedPoly(idx, V, BdyPtsxy);
   call Vor_SortPolyByAngle(Poly);
   print Poly;
*/
start Vor_ExtractUnboundedPoly(_vertIdx, V, BdyPts);
   vertIdx = _vertIdx[ ,loc(_vertIdx >= 1) ];  /* get rid of missing values */
   K = ncol(vertIdx);

   /* find Voronoi vertex at infinity */
   infinityIdx = loc(vertIdx=1);
   if infinityIdx = 1 then do;                   /* first index */
      Poly = BdyPts // V[vertIdx[2:K], ];
   end;
   else if infinityIdx = K then do;              /* last index */
      Poly = V[vertIdx[1:K-1], ] // BdyPts;
   end;
   else do;                                      /* last index */
      Poly = V[vertIdx[1:infinityIdx], ] // BdyPts // V[vertIdx[infinityIdx+1:K], ];
   end;
   return( Poly );
finish;

/* Main function that returns a set of polygons for the 
   bounded Voronoi tesellation within a rectangle.

   The function takes the following args:
   V     : Voronoi vertices
   Cells : Voronoi cells (indices for the rows of the vertices)
   CHS   : Sites that generate Voronoi cells; CHS[,1:2]=(x,y); CHS[,3]=ID
   BBox  : A bounding rectangle that includes all points. If not specified,
           use inflated bounding box of union of vertices and sites.

   EXAMPLE:
   sites = {1 1, 4 2, 3 4, 2 3};
   CHSites = {2 3, 3 4, 4 1, 1 2};
   V = {. ., 2.5 1.5, 3.1666667 2.8333333};
   VCells = {1 2 .,
             3 1 2,
             1 3 .,
             1 3 2};
   BBox = {0 0  5 5};
   CHSxy = sites[ CHSites[,1], ];
   CHS = CHSxy || CHSites[,1];     * (x,y,rowNum) ;
   Poly = Vor_BoundedVPoly(V, VCells, CHS, BBox);
   print Poly;
*/
start Vor_BoundedVPoly(V, Cells, CHS, BBox=) global(g_debugBVP);
   if ncol(g_debugBVP)=0 then 
      g_debugBVP = 0;
   DEBUG = g_debugBVP;
   CHSidx = CHS[,3]; CHSxy = CHS[,1:2];
   if nrow(Cells)<4 then do;
      /* special case: The Voronoi diagram for a triangle has
         one interior vertex with THREE rays that go to infinity.
         This might get complicated? Or maybe we can handle it? */
      print "SPECIAL VORONOI DIAGRAM. NOT YET SUPPORTED.";
      return ({. . .});
   end;

   /* if BBox not specified, compute BBox of sites and Voronoi vertices,
      then increase by 10% on all sides */
   if IsSkipped(BBox) then
      BBox =  Vor_InflateBBox(V // CHSxy);
   else do;      
      /* Check that all points are inside the BBox */
      if ^Vor_PtsInBBox(V, BBox) | ^Vor_PtsInBBox(CHSxy, BBox) then do;
         print "ERROR: In Vor_BoundedVPoly, bounding box must contain all points";
         return({. . .});
      end;
   end;
   if DEBUG then PRINT BBox[c={xmin ymin xmax ymax}];

   /* Allocate room for the polygons. If there are m sites,
      we expect 6*m vertices (on average) for the Voronoi polygons.
      Round up to 10, but monitor the number of sizes and reallocate,
      if necessary. */
   numPolyRows = 10*nrow(Cells);
   Poly = j(numPolyrows, 3, .);

   /* Every Voronoi diagram in the plane has K unbounded 
      cells, where K is the number of vertices on the convex hull (CH)
      of the sites. Intersect these unbounded cells with the 
      bounding box to create bounded polygons that we can plot. 
      The following call gets the intersection points with the BBox for 
      each cell. */
   BdyPts = Vor_CHRaysIntersectRect(CHS, BBox); /* (Cell, side, x, y, Angle) */
   if DEBUG then PRINT BdyPts[c={Cell side x y Angle}];

   /* loop over Voronoi cells. 
      If the cell does not contain vertex 1 (the point at 
      infinity), then it is a bounded polygon. 
      If the cell DOES contain vertex 1, then it is an unbounded 
      polygon. Use the BdryPts and the interior points to obtain a bounded polygon.
   */
   startRow = 1;
   do i = 1 to nrow(Cells);
      if DEBUG then PRINT i (Cells[i,])[L='Cells[i,]'];
      if any(Cells[i,] = 1) then do;
         idx = loc(BdyPts[,1]=i);  /* extract the boundary points for this cell */
         B = BdyPts[idx, 3:4];     /* (x,y) coords to replace for "point at infinity" */
         P = Vor_ExtractUnboundedPoly(Cells[i,], V, B);
      end;
      else 
         P = Vor_ExtractBoundedPoly(Cells[i,], V);
      P = P || j(nrow(P), 1, i);    /* append polygon ID */
      if DEBUG then PRINT i P[c={x y ID}];
      k = nrow(P);
      endRow = startRow + k - 1;
      if endRow > numPolyRows then do;
         Poly =  _DynamRealloc(Poly);
         numPolyRows = nrow(Poly);
      end;
      Poly[startRow:endRow, ] = P;
      startRow = endRow + 1;
   end;
   Poly = CompleteCases( Poly, "extract" );     /* remove unused rows */
   run Vor_SortPolyByAngle(Poly);
   return Poly;
finish;



/*********************************************************/

/* CG_EmptyCircle returns the center and radius of the largest empty circle
   whose center is interior to the CH of a set of sites in two dimensions.
   That is, the interior of the circle does not contain any of the sites.
   
   INPUT:
   sites : (n x 2) matrix of n d-dimensional points in general position without duplicates
   OUTPUT:
   c : (2 x 2) matrix where
       c[1,] is a Voronoi vertex
       c[2,1] is the radius of the largest circle centered at c[1,] that
              does not contain any sites in its interior.
*/
start CG_EmptyCircle(sites);
   n = nrow(sites);
   d = ncol(sites);
   if n < 3 | d ^= 2 then 
      stop "The sites must have two columns and at least three rows";
      
   result = j(2, d, .);   
   call voronoi(V, VCells, sites);
   /*
   call convexhull(indices, va, sites);
   P = sites[indices[,1],];
   */
   P = sites[ Compgeom_CHIndices(sites), ];

   /* PolyPtInside returns 
      1 if pt inside polygon
      . if pt on boundary
      0 if pt outside polygon
      We could handle Empty Spheres in 3-D, etc, if we write a geneal function that 
      detects whether a point is outside a convec hull.      
   */
   isInside = PolyPtInside(P, V);
   inNdx = loc(isInside ^= 0);    /* index of interior and boundary points */
   if ncol(inNdx)=0 then 
      /* something's wrong. Shouldn't happen */
      stop "All sites are outside the convex hull"; 
   
   Vin = V[inNdx, ];             /* keep k vertices that are in closure(CH) */
   dist = distance(Vin, sites);  /* k x n matrix of distances */
   /* generically, there are three columns in each row that are closest to vertex */
   minDist = dist[,><];          /* for each row, find min distance */
   r = minDist[<>];              /* radius of largest empty circle */
   maxNdx = minDist[<:>];        /* row index of vertex that is center of empty circle */
   center = Vin[maxNdx,];
   
   result[1,] = center;
   result[2,1] = r;
   return result;
finish;

start Compgeom_EquateAxesStr(sites, ninc=5);
   z = colvec(sites);
   call gscale(scale, z, ninc);
   s = strip(char(scale[1])) + " to " + 
       strip(char(scale[2])) + " by " + 
       strip(char(scale[3]));
   return s;
finish;

store module=(
   Compgeom_CHIndices 
   Compgeom_TriCentroid 
   _TriCentroid 
   Compgeom_InBBox 
   Compgeom_Barycentric 
   Compgeom_PtInTriangle
   Compgeom_Make2DTri
   Compgeom_EquateAxesStr
   Compgeom_NNDistance

   Vor_PtsInBBox
   Vor_RayIntersectRect
   Vor_WhichSideRect
   Vor_InflateBBox
   Vor_PtsInCommon2D
   _DynamRealloc
   Vor_CHRaysIntersectRect
   Vor_SortPolyByAngle
   Vor_ExtractBoundedPoly
   Vor_ExtractUnboundedPoly
   Vor_BoundedVPoly

   CG_EmptyCircle
   CG_FindPtsInTri
   CG_FindNearestNbr
   );

*QUIT;
