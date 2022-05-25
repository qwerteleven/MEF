

h = 1.e-3; w = 4.72e-3; t = 0.035e-3;
xBox = w/2. * 6.; yBox = h * 12.;


s = 1. ;

p0 = h /10.*s;
pLine0 = w/2. / 10. * s; pLine1 = w/2. / 50. * s;
pxBox = xBox / 10. * s ; pyBox = yBox / 8. * s;


Point(1) = { 0,    0, 0, p0};
Point(2) = { xBox, 0, 0, pxBox};
Point(3) = { xBox, h, 0, pxBox};
Point(4) = { 0,    h, 0, pLine0};
Point(5) = { w/2., h, 0, pLine1};
Point(6) = { 0,    h+t, 0, pLine0};
Point(7) = { w/2., h+t, 0, pLine1};
Point(8) = { 0,    yBox, 0, pyBox};
Point(9) = { xBox, yBox, 0, pyBox};



Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 9};
Line(4) = {9, 8};
Line(5) = {8, 6};

Line(7) = {4, 1};
Line(8) = {5, 3};
Line(9) = {4, 5};
Line(10) = {6, 7};
Line(11) = {5, 7};



Line Loop(12) = {8, -2, -1, -7, 9}; Plane Surface(13) = {12};
Line Loop(14) = {10, -11, 8, 3, 4, 5}; Plane Surface(15) = {14};



Physical Line(16) = {1, 2, 4};
Physical Line(17) = {9, 10};
Physical Line(18) = {5};

Physical Surface (19) = {13};
Physical Surface (20) = {15};

Physical Line(21) = {7};