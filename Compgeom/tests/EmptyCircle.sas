proc iml;
%include "CG_define.sas";
%include "CG_proc.sas";
QUIT;

/***********************************************************/
/* test the helper function CG_EmptyCircle */
/***********************************************************/

proc iml;
load module=(CG_EmptyCircle);

/* Example data from 
   https://youtu.be/dijZkOCNMo0?si=sa4G_Ggp_kDg9AY5
*/
sites = {-7  8, -9  4, -3 -4,  9  0, 9 12};
result = CG_EmptyCircle(sites);
center = result[1,];
radius = result[2,1];
print center[c={'x' 'y'} L=""] radius; 
QUIT;
