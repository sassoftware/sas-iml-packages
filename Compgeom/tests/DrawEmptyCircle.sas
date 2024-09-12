proc iml;
%include "CG_define.sas";
%include "CG_proc.sas";
QUIT;

/****************************************/
/* test the function CG_DrawEmptyCircle */
/****************************************/
proc iml;
load module=(CG_DrawEmptyCircle);

/* Example data from 
   https://youtu.be/dijZkOCNMo0?si=sa4G_Ggp_kDg9AY5
*/
sites = {-7  8, -9  4, -3 -4,  9  0, 9 12};

title "Largest Empty Circle Problem";
title2 "Find point in CH(sites) that is farthest from any site";
run CG_DrawEmptyCircle(sites);
QUIT;

/* Example from 
   https://blogs.sas.com/content/iml/2016/10/12/empty-space-distance-plot.html
*/
data BigCities;
set maps.uscity;                          /* part of SAS/GRAPH */
where statecode not in ("AK" "HI" "PR") & pop > 200000;
if pop>300000 then Name=City;             /* label really big cities */
else Name = "        ";
run;

title  "Large US Cities";
title2 "Population > 200,000";
footnote J=L "Center is in southern Montana";
proc sgplot data=BigCities;
  scatter x=x y=y / markerattrs=(symbol=CircleFilled) datalabel=Name;
  xaxis grid;
  yaxis grid;
run;

   
proc iml;
load module=(CG_DrawEmptyCircle);
use BigCities;
read all var {x y} into sites;
close;

title "Largest Empty Circle for Large US Cities";
title2 "Find point in CH farthest from any large city";
footnote J=L "Center is in southern Montana";
run CG_DrawEmptyCircle(sites);
footnote;
QUIT;
