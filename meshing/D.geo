// Gmsh project created on Thu Mar 18 14:46:02 2021
SetFactory("OpenCASCADE");
//+
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {500, 0, 0, 1.0};
//+
Point(3) = {500, 500, 0, 1.0};
//+
Point(4) = {0, 500, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Physical Curve(1) = {1};
//+
Physical Curve(2) = {2};
//+
Physical Curve(3) = {3};
//+
Physical Curve(4) = {4};
//+
Physical Surface(5) = {1};
