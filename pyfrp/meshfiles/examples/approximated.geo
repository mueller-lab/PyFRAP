//vertices
Point(1)= {6.12323399574e-17,-1.0,0.0,1};
Point(2)= {0.15643446504,-0.987688340595,0.0,1};
Point(3)= {0.309016994375,-0.951056516295,0.0,1};
Point(4)= {0.45399049974,-0.891006524188,0.0,1};
Point(5)= {0.587785252292,-0.809016994375,0.0,1};
Point(6)= {0.707106781187,-0.707106781187,0.0,1};
Point(7)= {0.809016994375,-0.587785252292,0.0,1};
Point(8)= {0.891006524188,-0.45399049974,0.0,1};
Point(9)= {0.951056516295,-0.309016994375,0.0,1};
Point(10)= {0.987688340595,-0.15643446504,0.0,1};
Point(11)= {1.0,0.0,0.0,1};
Point(12)= {0.987688340595,0.15643446504,0.0,1};
Point(13)= {0.951056516295,0.309016994375,0.0,1};
Point(14)= {0.891006524188,0.45399049974,0.0,1};
Point(15)= {0.809016994375,0.587785252292,0.0,1};
Point(16)= {0.707106781187,0.707106781187,0.0,1};
Point(17)= {0.587785252292,0.809016994375,0.0,1};
Point(18)= {0.45399049974,0.891006524188,0.0,1};
Point(19)= {0.309016994375,0.951056516295,0.0,1};
Point(20)= {0.15643446504,0.987688340595,0.0,1};
Point(21)= {6.12323399574e-17,1.0,0.0,1};
Point(22)= {-1,0,0,0.01};

//lines
Line(43)= {21,22};
Line(44)= {22,1};

//arcs

//bSplines
BSpline(45)= {1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,21};

//lineLoops
Line Loop(44)= {45,43,44};

//ruledSurfaces

//surfaceLoops

//volumes

//fields

