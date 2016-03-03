from fipy import *


volSize_px=50
radius=300
height=1000

center_embr_px_x=256
center_embr_px_y=256

mesh = Gmsh3D('''
	//Parameters
	volSize_px = %(volSize_px)g; 
	radius = %(radius)g;
	height = %(height)g;
	center_embr_px_x = %(center_embr_px_x)g;
	center_embr_px_y = %(center_embr_px_y)g;
	//Points of upper circle
	Point(1) = {center_embr_px_x, center_embr_px_y, 0, volSize_px}; 
	Point(2) = {-radius+center_embr_px_x, center_embr_px_y, 0, volSize_px};
	Point(3) = {center_embr_px_x, center_embr_px_y+radius, 0, volSize_px};
	Point(4) = {center_embr_px_x+radius, center_embr_px_y, 0, volSize_px};
	Point(5) = {center_embr_px_x, -radius+center_embr_px_y, 0, volSize_px};
	//Points of lower circle
	Point(6) = {center_embr_px_x, center_embr_px_y, -height, volSize_px}; 
	Point(7) = {-radius+center_embr_px_x, center_embr_px_y, -height, volSize_px};
	Point(8) = {center_embr_px_x, center_embr_px_y+radius, -height, volSize_px};
	Point(9) = {center_embr_px_x+radius, center_embr_px_y, -height, volSize_px};
	Point(10) = {center_embr_px_x, -radius+center_embr_px_y, -height, volSize_px};
	//Circles
	Circle(41) = {2, 1, 3}; 
	Circle(42) = {3, 1, 4}; 
	Circle(43) = {4, 1, 5}; 
	Circle(44) = {5, 1, 2};
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
	''' % locals()) 

x,y,z=mesh.cellCenters

print len(x)