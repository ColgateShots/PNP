// Gmsh project created on Tue Jan 19 17:45:28 2021
SetFactory("OpenCASCADE");
//+
Point(1) = {500, 500, 0, 1};
//+
Point(2) = {-500, 500, 0, 1};
//+
Point(3) = {-500, -500, 0, 1};
//+
Point(4) = {500, -500, 0, 1};
//+
Point(5) = {0, 0, 0, 1};
//+
Point(6) = {111.6, 500, 0, 1};
//+
Point(7) = {-111.6, 500, 0, 1};
//+
Line(1) = {6, 1};
//+
Line(2) = {1, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 2};
//+
Line(5) = {2, 7};
//+
Spline(6) = {7, 5, 6};
//+
Curve Loop(1) = {1, 2, 3, 4, 5, 6};
//+
Plane Surface(1) = {1};
//+
Physical Curve(1) = {1, 2, 3, 4, 5, 6};
//+
Physical Surface(2) = {1};
