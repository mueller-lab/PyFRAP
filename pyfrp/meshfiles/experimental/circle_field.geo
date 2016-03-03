volSize_px=24;
center_x=256.5;
center_y=256.5;
radius=289.211219521;

// Upper circle
Point(1) = {center_x, center_y, 0, volSize_px}; 
Point(2) = {-radius+center_x, center_y, 0, volSize_px};
Point(3) = {center_x, center_y+radius, 0, volSize_px};
Point(4) = {center_x+radius, center_y, 0, volSize_px};
Point(5) = {center_x, -radius+center_y, 0, volSize_px};

//Circles describing upper circle
Circle(41) = {2, 1, 3}; 
Circle(42) = {3, 1, 4}; 
Circle(43) = {4, 1, 5}; 
Circle(44) = {5, 1, 2};

Line Loop(61) = {42, 43, 44, 41};
Plane Surface(62) = {61};

Field[1] = Attractor;
Field[1].NodesList = {1};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 5;
Field[2].LcMax = 24;
Field[2].DistMin = 50;
Field[2].DistMax = 100;

Field[3] = Min;
Field[3].FieldsList = {2};
Background Field = 3;
