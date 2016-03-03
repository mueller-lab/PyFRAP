//Parameters
volSize_px=40;
volSize_fine=2.5;
center_x=256.5;
center_y=256.5;
radius=200.863218314;
height=90.3332804037;
sl_refined=20;
sl_height=36.1333121615;
sl_width=5.;

// Upper circle
Point(1) = {center_x, center_y, 0, volSize_px}; 
Point(2) = {-radius+center_x, center_y, 0, volSize_px};
Point(3) = {center_x, center_y+radius, 0, volSize_px};
Point(4) = {center_x+radius, center_y, 0, volSize_px};
Point(5) = {center_x, -radius+center_y, 0, volSize_px};

//Lower Circle
Point(6) = {center_x, center_y, -height, volSize_px}; 
Point(7) = {-radius+center_x, center_y, -height, volSize_px};
Point(8) = {center_x, center_y+radius, -height, volSize_px};
Point(9) = {center_x+radius, center_y, -height, volSize_px};
Point(10) = {center_x, -radius+center_y, -height, volSize_px};

//Circles describing upper circle
Circle(41) = {2, 1, 3}; 
Circle(42) = {3, 1, 4}; 
Circle(43) = {4, 1, 5}; 
Circle(44) = {5, 1, 2};

//Circles describing lower circle
Circle(45) = {7, 6, 8}; 
Circle(46) = {8, 6, 9}; 
Circle(47) = {9, 6, 10}; 
Circle(48) = {10, 6, 7};

//Lines connecting circles 
Line(49) = {2,7};
Line(50) = {3,8};
Line(51) = {4,9};
Line(52) = {5,10};

//Surfaces for disc
Line Loop(53) = {44, 49, -48, -52};
Ruled Surface(54) = {53};
Line Loop(55) = {52, -47, -51, 43};
Ruled Surface(56) = {55};
Line Loop(57) = {42, 51, -46, -50};
Ruled Surface(58) = {57};
Line Loop(59) = {50, -45, -49, 41};
Ruled Surface(60) = {59};
Line Loop(61) = {42, 43, 44, 41};
Ruled Surface(62) = {61};
Line Loop(63) = {45, 46, 47, 48};
Ruled Surface(64) = {63};
Surface Loop(65) = {64, 60, 58, 62, 56, 54};

//Volumes
Volume(80) = {65};




Field[5] = Box;
Field[5].VIn = volSize_fine*3;
Field[5].VOut = volSize_px/3;
Field[5].XMin = center_x-2*sl_refined; 
Field[5].XMax = center_x+2*sl_refined;
Field[5].YMin = center_y-2*sl_refined;
Field[5].YMax = center_y+2*sl_refined;
Field[5].ZMin = -sl_height-2*sl_width;
Field[5].ZMax = -sl_height+2*sl_width;

Field[6] = Box;
Field[6].VIn = volSize_fine;
Field[6].VOut = volSize_fine*3;
Field[6].XMin = center_x-sl_refined; 
Field[6].XMax = center_x+sl_refined;
Field[6].YMin = center_y-sl_refined;
Field[6].YMax = center_y+sl_refined;
Field[6].ZMin = -sl_height-sl_width;
Field[6].ZMax = -sl_height+sl_width;


/*
Field[3] = MathEval;
Field[3].F = (Sqrt((x-center_x)^2+(y-center_y)^2)/radius)*volSize_px+volSize_fine;*/


// Field[1] = Attractor;
// Field[1].NodesList = {1,2,3,4,5};
// Field[1].NNodesByEdge = 50;
// Field[1].EdgesList = {41,42,43,44};
// 
// Field[2] = Threshold;
// Field[2].IField = 1;
// Field[2].LcMin = 10;
// Field[2].LcMax = 20;
// Field[2].DistMin = 20;
// Field[2].DistMax = 50;

/*
Field[3] = Box;
Field[3].VIn = 2*volSize_fine;
Field[3].VOut = volSize_px;
Field[3].XMin = center_x-sl_refined; 
Field[3].XMax = center_x+sl_refined;
Field[3].YMin = center_y-sl_refined;
Field[3].YMax = center_y+sl_refined;
Field[3].ZMin = -2*sl_height-sl_width;
Field[3].ZMax = -2*sl_height+sl_width;*/

Mesh.CharacteristicLengthExtendFromBoundary = 1;
Field[8] = Max;
Field[8].FieldsList = {5};

Background Field =8;







