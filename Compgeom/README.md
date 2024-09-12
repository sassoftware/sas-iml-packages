# The Compgeom repo

## Description

This project is a set of helper routines for computational geometry in two dimensions. The functions are written in the SAS IML language. They are intended to show users how to use the built-in computational geometry subroutines: CONVEXHULL, DELAUNAY, and VORONOI.  Of particular interest are the subroutines that visualize these geometric objects: The **CG_DrawConvexHull** routine, the **CG_DrawDelaunay** routine, and the **CG_DrawVoronoi** routine.

The SAS IML functions are described in the documentation for the Compgeom package.

## Main functions

The following functions help to solve or visualize of computational geometry problems in 2-D. In the following list, the _sites_ are a set of unique 2-D points that generate the convex hull (CH), the Delaunay tringulation, and the Voronoi diagram.

- **CG_DrawConvexHull**: Create a graph that visualizes the convex hull of the _sites_ matrix. Optionally, display and label the sites that are on or inside the convex hull.
- **CG_DrawDelaunay**: Call the DELAUNAY function to obtain the Delaunay triangulation of the _sites_ matrix. Visualize the Delaunay triangulation. Optionally, display and label the sites and the triangles.
- **CG_DrawEmptyCircle**: Display a scatter plot of the sites and overlay the largest circle that contains no sites and that is centered inside the convex hull of the sites. 
- **CG_DrawTri**: Create a graph that visualizes a triangulation. Optionally, display and label the sites and label the triangles.
- **CG_DrawVoronoi** Call the VORONOI function to obtain the Voronoi tessellation of the _sites_ matrix. Visualize the Voronoi cells. 
- **CG_EmptyCircle**: Find the largest circle that contains no sites and that is centered inside the convex hull of the sites. Geometrically, the center is one of the Voronoi vertices determined by the set of sites.
- **CG_FindNearestNbr**: Let P be an  Nxd matrix and let i be any row of P. Return j &ne; i such that dist(P[j,], P[i,]) is as small as possible. 
- **CG_FindPtsInTri**: For each query point, _q_, find the triangle in a triangulation that contains _q_. If _q_ is not in any triangle, return a missing value.

## Documentation

The file CompGeom_Doc.docx is a Word file that describes the syntax of each public function. The documentation shows how to call the functions and provides examples of output.
