// Gmsh project created on Thu Mar  4 11:17:11 2021
SetFactory("OpenCASCADE");

//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {500, 0, 0, 1.0};
//+
Point(3) = {500, 1000, 0, 1.0};
//+
Point(4) = {223.9, 1000, 0, 1.0};
//+
Point(5) = {0, 500, 0, 1.0};
//+
Point(6) = {50, 525, 0, 1.0};
//+
Point(7) = {90, 581, 0, 1.0};
//+
Point(8) = {110, 621, 0, 1.0};
//+
Point(9) = {220, 984, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Spline(4) = {4,9,8,7,6,5};
//+
Line(5) = {5, 1};
//+
Curve Loop(1) = {3, 4, 5, 1, 2};
//+
Plane Surface(1) = {1};
//+
Physical Curve(1) = {3};
//+
Physical Curve(2) = {2};
//+
Physical Curve(3) = {1};
//+
Physical Curve(4) = {5};
//+
Physical Curve(5) = {4};
//+
Physical Surface(6) = {1};

