*proc iml;
/* Visualization modules for PROC IML:
   - CG_DrawConvexHull
   - CG_DrawDelaunay
   - CG_DrawEmptyCircle
   - CG_DrawTri
   - CG_DrawVoronoi
*/

/* Visualize a 2-D convex hull of n sites in the plane.
   Compute the convex hull by calling the CONVEXHULL subroutine, 
   then call PROC SGPLOT to draw the CH and (optionally) vertices, 
   interior points & other points.
   SYNTAX:
   CG_DrawConvexHull(sites, opt=) 
   sites : an (n x 2) matrix of 2-D points that define the CH
   opt   : an option 3-element vector that specifies visualization options
           Default is opt = {1, 1, 0.5}
         opt[1] : Show the sites that are vertices on the CH
                  0: Do not display vertices
                  1: Display vertices (default)
                  2: Display vertices and label them by using the rows of sites
         opt[2] : Show the sites that are interior to the CH
                  0: Do not display interior sites
                  1: Display interior sites (default)
                  2: Display interior sites and label them by using the rows of sites
         opt[3] : Transparency level for the CH polygon
                  alpha: (0<=alpha<=1) Display fill with transparency alpha (default=0.5)
                  0: 0% transparent ==> 100% opaque ==> Display opaque fill
                  1: 100% transparent ==> No fill ==> Display only the outline
         opt[4] : Display PROC SGPLOT statements in log
                  0: Do not display PROC SGPLOT statements (default)
                  1: Display PROC SGPLOT statements
*/ 
start CG_DrawConvexHull(sites, opt=, debug=0);
   /* get indices for the vertices of the convex hull of the sites */
   n = nrow(sites);
   CHIndices = CompGeom_CHIndices(sites);
   CHVerts = sites[CHIndices, ];
   Poly = CHVerts || j(nrow(CHVerts),1,1) || CHIndices;
   create _Poly from Poly[c={'vx' 'vy' 'vID' 'vLabel'}];
   append from Poly;
   close;
   DSSETSTMT = "set _Poly";

   option = {1, 1, 0.5, 0}; 
   if ^IsSkipped(opt) then do;
      k = 1:min(nrow(option),nrow(opt)*ncol(opt));
      option[k] = opt[k];
   end;
   /* set properties for displaying vertices */
   ShowCHVerts = option[1];
   if ShowCHVerts=. | ShowCHVerts=0 then 
      CHVERTSSTMT = " ";
   else if ShowCHVerts>=1 then /* display vertices */
      CHVERTSSTMT = "scatter x=vx y=vy / markerattrs=GRAPHDATA1(symbol=CircleFilled)";
   if ShowCHVerts=2 then 
      CHVERTSSTMT = CHVERTSSTMT + " datalabel=vLabel labelstrip";

   /* set properties for displaying interior sites */
   showIntSites = option[2];
   if showIntSites=. | showIntSites=0 then 
      INTSITESSTMT = " ";
   else if showIntSites>=1 then /* display vertices */
      INTSITESSTMT = "scatter x=sx y=sy / markerattrs=GRAPHDATA2";
   if showIntSites=2 then
      INTSITESSTMT = INTSITESSTMT + " datalabel=sLabel labelstrip";
  
   if showIntSites then do;
      intIdx = colvec( setdif(1:n, CHIndices) );
      intSites = sites[intIdx,];
      create _IntSites from intSites intIdx[c={'sx' 'sy' 'sLabel'}];
      append from intSites intIdx;
      close;
      DSSETSTMT = DSSETSTMT + " _IntSites";
   end;

   /* set properties for displaying CH polygon (NONE, FILL, or transparency) */
   transparency = option[3];
   if transparency < 0 then transparency=0;
   else if transparency > 1 then transparency=1;
   /* if transparency = 1, then show only the outline */
   POLYOPTSTMT = "outline transparency=0 lineattrs=(thickness=2)";
   if transparency=0 then 
      POLYOPTSTMT = POLYOPTSTMT + " fill";
   else if (0 < transparency & transparency < 1) then 
      POLYOPTSTMT = POLYOPTSTMT + " fill fillattrs=(transparency=" 
                    + strip(char(transparency)) + ")";

   showStmts = option[4];
   if showStmts=. | showStmts=0 then 
      SHOWSTMT = " ";
   else 
      SHOWSTMT = "options source";
 
   if debug then DEBUGSTMT = "options source";
   else DEBUGSTMT = "options nosource";

   submit dsSetStmt CHVERTSSTMT INTSITESSTMT POLYOPTSTMT DEBUGSTMT SHOWSTMT;
      &DEBUGSTMT;
      data _DCH;
      &dsSetStmt;
      run;

      &SHOWSTMT;
      proc sgplot data=_DCH noautolegend;
         polygon x=vx y=vy ID=vID / &POLYOPTSTMT;
         &INTSITESSTMT;
         &CHVERTSSTMT;
         xaxis grid offsetmin=0.1 offsetmax=0.1 label="x";
         yaxis grid offsetmin=0.1 offsetmax=0.1 label="y";
      run;
   endsubmit;
finish;

/* Visualize a triangulation.
   SYNTAX:
   CG_DrawTri(sites, tri, opt=) 
   sites : an (n x 2) matrix of 2-D points that define the (x,y) coords of vertices
   tri   : a (k x 3} matrix of integers such that each row of tri 
         specifies a triangle. The i-th triangle has vertices
         sites[ tri[i,], ]

   opt   : an option vector that specifies visualization options
           Default is opt = {1, 1, 0.5, 0}
         opt[1] : Show vertices of the triangulation
                  0: Do not display vertices
                  1: Display vertices (default)
                  2: Display vertices and label them by using the rows of sites
         opt[2] : Label the triangles at their centroid
                  0: Do not label triangles
                  1: Label the triangles (default)
         opt[3] : Transparency level for the polygons
                  alpha: (0<=alpha<=1) Display fill with transparency alpha (default=0.5)
                  0: 0% transparent ==> 100% opaque ==> Display opaque fill
                  1: 100% transparent ==> No fill ==> Display only the outline
         opt[4] : Display PROC SGPLOT statements in log
                  0: Do not display PROC SGPLOT statements (default)
                  1: Display PROC SGPLOT statements
*/
start CG_DrawTri(sites, tri, opt=, debug=0);
   /* write vertices and ID */
   vID = T(1:nrow(sites));
   create _Vert from sites vID[c={'vx' 'vy' 'vLabel'}];
      append from sites vID;
   close;

   /* triangulation */
   P = sites[tri, ];
   ID = colvec( repeat( T(1:nrow(tri)), 1, 3) );
   create _Tri from P ID[c={'tx' 'ty' 'tID'}];
      append from P ID;
   close;

   /* centroids of each triangle */
   Centroids = Compgeom_TriCentroid(P||ID);   /* (x,y,ID) */
   create _Centroids from Centroids[c={'cx' 'cy' 'cID'}];
      append from Centroids;
   close;

   option = {1, 1, 0.5, 0}; 
   if ^IsSkipped(opt) then do;
      k = 1:min(nrow(option), nrow(opt)*ncol(opt));
      option[k] = opt[k];
   end;
   if debug then 
      print option;

   /* set properties for displaying vertices */
   VERTSSTMT = " ";
   ShowVerts = option[1];
   if ShowVerts>=1 then /* display vertices */
      VERTSSTMT = "scatter x=vx y=vy / markerattrs=GRAPHDATA1(symbol=CircleFilled)";
   if ShowVerts=2 then 
      VERTSSTMT = VERTSSTMT + " datalabel=vLabel labelstrip datalabelattrs=GRAPHDATA1(size=9)";

   /* label triangles at centroid */
   LABELTRISTMT = " ";
   LabelTri = option[2];
   if LabelTri>=1 then /* display vertices */
      LABELTRISTMT = "text x=cx y=cy text=CID / strip textattrs=GRAPHDATA2(size=12)";

   /* set properties for displaying polygon (NONE, FILL, or transparency) */
   transparency = option[3];
   if transparency < 0 then transparency=0;
   else if transparency > 1 then transparency=1;
   /* if transparency = 1, then show only the outline */
   POLYOPTSTMT = "outline transparency=0 lineattrs=(thickness=1)";
   if transparency=0 then 
      POLYOPTSTMT = POLYOPTSTMT + " fill";
   else if (0 < transparency & transparency < 1) then 
      POLYOPTSTMT = POLYOPTSTMT + " fill fillattrs=(transparency=" 
                    + strip(char(transparency)) + ")";

   SHOWSTMT = " ";
   showStmts = option[4];
   if showStmts >= 1 then 
      SHOWSTMT = "options source";
 
   if debug then DEBUGSTMT = "options source";
   else DEBUGSTMT = "options nosource";

   submit VERTSSTMT LABELTRISTMT POLYOPTSTMT DEBUGSTMT SHOWSTMT;
      &DEBUGSTMT;
      data _DTri;
      set _Vert _Tri _Centroids;
      run;

      &SHOWSTMT;
      proc sgplot data=_DTri noautolegend; 
         polygon x=tx y=ty ID=tID / &POLYOPTSTMT;
         &VERTSSTMT;     /* display vertices */
         &LABELTRISTMT;  /* label triangles */
         xaxis grid offsetmin=0.05 offsetmax=0.05 label="x";
         yaxis grid offsetmin=0.05 offsetmax=0.05 label="y";
      run;
   endsubmit;    
finish;


/* Compute a Delaunay triangulation and then call CG_DrawTri */
start CG_DrawDelaunay(sites, opt=, debug=0);
   tri = delaunay(sites);
   option = {1, 1, 0.5, 0}; 
   if ^IsSkipped(opt) then do;
      k = 1:min(nrow(option), nrow(opt)*ncol(opt));
      option[k] = opt[k];
   end;
   call CG_DrawTri(sites, tri, option, debug);
finish;

/* CG_DrawEmptyCircle: 
   Call the CG_EmptyCircle function to get the center and radius.
   Call Compgeom_EquateAxesStr to figure out how to scale the axes.
   Overlay the empty circle on the scatter plot of the sites.
*/
start CG_DrawEmptyCircle(sites);
   result = CG_EmptyCircle(sites);
   center = result[1,];
   r = result[2,1];
   *print center[c={'x' 'y'}] r;   

   /* graph the largest empty circle and the sites */
   create _EmptyCircle from sites center r [c={'x' 'y' 'cx' 'cy' 'r'}];
   append from sites center r;
   close;
      
   valueList = Compgeom_EquateAxesStr(sites);

   submit valueList;
   proc sgplot data=_EmptyCircle aspect=1 noautolegend;
      ellipseparm semimajor=r semiminor=r / xorigin=cx yorigin=cy fill fillattrs=(transparency=0.3) outline;
      scatter x=x y=y / markerattrs=(symbol=X);
      scatter x=cx y=cy / markerattrs=(symbol=CircleFilled);
      xaxis grid values=(&valuelist);
      yaxis grid values=(&valuelist);
   run;
   endsubmit;
finish;

/* Compgeom_DrawVoronoiNoCells: 
   Create a barebones diagram that has sites, CH of sites, 
   and Voronoi vertices, but not the Voronoi cells.

   EXAMPLE:
   sites = {0 0, 3 -1, 2 3, 1 4, 1 1};
   BBox = {-2 -2  4 5};
   CHSites = {1 2, 2 3, 3 4, 4 1};
   V = {. ., 1.5 -0.5, -1.5 2.5, 3.1666667 1.1666667, 0.5 2.5};
   VCells = {
         3 1 2 .,
         1 2 4 .,
         1 5 4 .,
         1 3 5 .,
         5 3 2 4};
   run Compgeom_DrawVoronoiNoCells(sites, BBox, V, VCells, CHsites); 
   run Compgeom_DrawVoronoiNoCells(sites, , V, VCells, CHsites); 
   run Compgeom_DrawVoronoiNoCells(sites, BBox, V, VCells, CHsites); 
*/
start Compgeom_DrawVoronoiNoCells(sites, BBox=, V=, VCells=, CHsites=);
%IF %sysevalf(&sysver = 9.4) %THEN %DO;
   if isSkipped(V) | isSkipped(VCells) | isSkipped(CHsites) then do;
      print "In SAS 9, you must provide the Vononoi vertices, cells, and CH of sites";
      return;
   end;
%END;
%ELSE %DO;
   if isSkipped(V) then 
      call voronoi(V, VCells, sites);
   if isSkipped(CHsites) then 
      call convexhull(CHsites, volume_area, sites) ;
%END;
   Sp = sites || T(1:nrow(sites));
   create _Sites from Sp[c={'sx' 'sy' 'sID'}]; append from Sp; close;

   Vp = V || T(1:nrow(V));
   create _VVerts from Vp[c={'vx' 'vy' 'vID'}]; append from Vp; close;

   CHSxy = sites[ CHSites[,1], ];   *(x,y) pairs;
   CHS = CHSxy || CHSites[,1];      *(x,y,rowNum);
   CHp = CHSites[,1] || CHSxy || j(nrow(CHSites),1,1);
   create _CHSites from CHp[c={'CHLabel' 'chx' 'chy' 'CHID' }]; 
      append from CHp; close;

   if isSkipped(BBox) then 
      BBox = InflateBBox(V // CHSxy);
   BBp = BBox[,{1 2}] //  BBox[,{3 2}] //
         BBox[,{3 4}] //  BBox[,{1 4}];
   BBp = BBp || {1,1,1,1};
   create _BBox from BBp[c={'bx' 'by' 'BID'}];   append from BBp; close;
   aspect = (BBox[3]-BBox[1]) / (BBox[4]-BBox[2]);
   if 0.5 <= aspect | aspect <= 2 then 
       aspectOpt = strip("aspect=" + char(1/aspect));
   else 
       aspectopt=" ";

   submit aspectOpt;
   data _Voronoi;
      format sID VID 2.;
      set _Sites _BBox _VVerts _CHSites;
   run; 

   proc sgplot data=_Voronoi &aspectOpt;
      scatter x=sx y=sy / datalabel=sID labelstrip markerattrs=GraphData2(symbol=X size=10) 
                          legendlabel="Sites" name="sites";
      /* Voronoi cells */
      scatter x=vx y=vy / datalabel=VID labelstrip markerattrs=GraphData3(symbol=CircleFilled size=12) 
                          legendlabel="Voronoi vertices" name="voronoi";
      keylegend "sites" "voronoi";
      xaxis grid label="x";
      yaxis grid label="y";
   run; 
   endsubmit;
finish;

/* CG_DrawVoronoi: Create a Voroinoi diagram that has sites, CH of sites, 
   Voronoi vertices, and (bounded) Voronoi cells.

   EXAMPLE:
   sites = {0 0, 3 -1, 2 3, 1 4, 1 1};
   BBox = {-2 -2  4 5};
   CHSites = {1 2, 2 3, 3 4, 4 1};
   V = {. ., 1.5 -0.5, -1.5 2.5, 3.1666667 1.1666667, 0.5 2.5};
   VCells = {
         3 1 2 .,
         1 2 4 .,
         1 5 4 .,
         1 3 5 .,
         5 3 2 4};
   run CG_DrawVoronoi(sites, BBox, V, VCells, CHsites); 
*/
start CG_DrawVoronoi(sites, BBox=, V=, VCells=, CHsites=, DEBUG=0);
%IF %sysevalf(&sysver = 9.4) %THEN %DO;
   if isSkipped(V) | isSkipped(VCells) | isSkipped(CHsites) then do;
      print "In SAS 9, you must provide the Vononoi vertices, cells, and CH of sites";
      return;
   end;
%END;
%ELSE %DO;
   if isSkipped(V) then 
      call voronoi(V, VCells, sites);
   if isSkipped(CHsites) then 
      call convexhull(CHsites, volume_area, sites) ;
%END;
   Sp = sites || T(1:nrow(sites));
   if DEBUG then 
      print Sp[c={'sx' 'sy' 'sID'}];
   create _Sites from Sp[c={'sx' 'sy' 'sID'}]; append from Sp; close;

   Vp = V || T(1:nrow(V));
   if DEBUG then 
      print Vp[c={'vx' 'vy' 'vID'}];
   create _VVerts from Vp[c={'vx' 'vy' 'vID'}]; append from Vp; close;

   CHSxy = sites[ CHSites[,1], ];   *(x,y) pairs;
   CHS = CHSxy || CHSites[,1];      *(x,y,rowNum);
   CHp = CHSites[,1] || CHSxy || j(nrow(CHSites),1,1);
   if DEBUG then 
      print CHp[c={'CHLabel' 'chx' 'chy' 'CHID' }];
   create _CHSites from CHp[c={'CHLabel' 'chx' 'chy' 'CHID' }]; 
      append from CHp; close;

   if isSkipped(BBox) then 
      BBox = Vor_InflateBBox(V // CHSxy);
   BB = rowvec(BBox);
   BBp = BB[,{1 2}] //  BB[,{3 2}] //
         BB[,{3 4}] //  BB[,{1 4}];
   BBp = BBp || {1,1,1,1};
   if DEBUG then 
      print BBp[c={'bx' 'by' 'BID'}];
   create _BBox from BBp[c={'bx' 'by' 'BID'}];   append from BBp; close;

   Poly = Vor_BoundedVPoly(V, VCells, CHS, BBox);
   if DEBUG then 
      print Poly[c={'vcx' 'vcy' 'VCID' }];
   create _BVoronoi from Poly[c={'vcx' 'vcy' 'VCID' }];
   append from Poly;
   close;
   aspect = (BBox[3]-BBox[1]) / (BBox[4]-BBox[2]);
   if 0.5 <= aspect | aspect <= 2 then 
       aspectOpt = strip("aspect=" + char(1/aspect));
   else 
       aspectopt=" ";

   if DEBUG then DEBUGSTMT = "options source";
   else DEBUGSTMT = "options nosource";
   
   submit aspectOpt DEBUGSTMT;
   &DEBUGSTMT;
   data _Voronoi;
      format sID VID 2.;
      set _Sites _BBox _VVerts _CHSites _BVoronoi;
   run; 

   proc sgplot data=_Voronoi &aspectOpt;
      /*CH*/
      *polygon id=CHID x=CHx y=CHy / fill outline fillattrs=(color=LightGray transparency=0.8); 
      /* Voronoi cells */
      scatter x=vx y=vy / datalabel=VID labelstrip markerattrs=GraphData1(symbol=CircleFilled size=12) legendlabel="Voronoi vertices" name="voronoi";
      polygon id=VCID x=vcx y=vcy / nofill lineattrs=GraphData1; 
      /* sites */
      scatter x=sx y=sy / datalabel=sID labelstrip markerattrs=GraphData2(symbol=X size=10) legendlabel="Sites" name="sites";
      /*add BBox*/ 
      *polygon id=BID x=bx y=by;
      *scatter x=bx y=by / markerattrs=(symbol=Square); 
      keylegend "sites" "voronoi";
      xaxis grid label="x" offsetmin=0.01 offsetmax=0.01;
      yaxis grid label="y" offsetmin=0.01 offsetmax=0.01;
   run; 
   endsubmit;
finish;

store module=(
   Compgeom_DrawVoronoiNoCells 
   CG_DrawConvexHull 
   CG_DrawDelaunay 
   CG_DrawEmptyCircle
   CG_DrawTri
   CG_DrawVoronoi
);
*QUIT;

