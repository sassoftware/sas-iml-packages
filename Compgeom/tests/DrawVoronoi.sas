proc iml;
%include "CG_define.sas";
%include "CG_proc.sas";
QUIT;

/* The DOC example */
proc iml;
load module=_all_;

/* doc example */
sites  = {1 2, 2.5 3, 3 8, 4 5.5, 5.5 3, 
          6.5 6, 6 1, 8 3};
ods graphics / width=400px height=400px;
title "Voronoi Diagram";
run CG_DrawVoronoi(sites);

/* zoom out */
BBox = {-5 -5 15 15};   
title "Voronoi Diagram (Zoom Out)";
run CG_DrawVoronoi(sites, BBox); 
QUIT;
