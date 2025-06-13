// -----------------------------------------------------------------------------

// Mesh with fine sampling at the y boundaries

// -----------------------------------------------------------------------------
lc = 0.01;

Point(1) = {0.5, -1, 0, lc};
Point(2) = {0.5, 1, 0, lc};
Point(3) = {1.5, 1, 0, lc};
Point(4) = {1.5, -1, 0, lc};

Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};

Curve Loop(5) = {1,2,3,4}; Plane Surface(6) = {5};
// Physical Surface("mesh") = {6};

Field[1] = Distance;
Field[1].CurvesList = {2,4};
Field[1].Sampling = 100;


// We then define a `Threshold' field, which uses the return value of the
// `Distance' field 1 in order to define a simple change in element size
// depending on the computed distances
//
// SizeMax -                     /------------------
//                              /
//                             /
//                            /
// SizeMin -o----------------/
//          |                |    |
//        Point         DistMin  DistMax
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc / 10;
Field[2].SizeMax = lc;
Field[2].DistMin = 0.02;
Field[2].DistMax = 0.3;

Background Field = 2;

// This will prevent over-refinement due to small mesh sizes on the boundary.

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm = 5;