#pragma once

#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>

namespace Another
{
	#define WIREFRAME 5
	#define FLAT 4
	#define SOLID 6
	#define SMOOTH 7
	#define EDGE -1
	#define FACE -2
	#define CORNER -3
	#define UNDEFINED -5
	#define PI 3.14159265358979323846
	#define POLYCUBEDRAW 10
	#define ORIGINALDRAW 11
	#define BOTH            12
	#define DEBUGMODE       13
	#define NONDEBUGMODE    14
	#define BOTHDRAW 15
	#define ZEROANGLE -1 
	#define PIANGLE 1 
	#define NORMALANGLE 0

	const int SELECTION_MODE = 15;

	typedef CGAL::Simple_cartesian<double>  Kernel;
	typedef Kernel::Point_3         Point;
	typedef Kernel::Vector_3        Vector_3;
	typedef Kernel::Vector_3        Normal;
	typedef Kernel::Plane_3         Plane;
	typedef Kernel::Point_2         PlanePoint;
	typedef Kernel::Segment_3         Segment;


	struct  Another_Traits
	{
		typedef Kernel::Point_3 Point_3; 
		typedef Kernel::FT FT;
		typedef Kernel::Vector_3 Vector_3; 
		typedef Kernel::Plane_3 Plane_3; 
		typedef Kernel::Point_2		Point_2;
		typedef Kernel::Triangle_3 Triangle_3;
		typedef Kernel::Triangle_2 Triangle_2;
		typedef Kernel::Segment_3 Segment_3;
		typedef Kernel::Segment_2 Segment_2;
		typedef Kernel::Ray_2 Ray_2;
		typedef Kernel::Ray_3 Ray_3;
		typedef Kernel::Vector_2 Vector_2;
		typedef Kernel::Tetrahedron_3 Tetrahedron_3;
		typedef Kernel::Sphere_3 Sphere_3;
		typedef Kernel::Object_2 Object_2;
		typedef Kernel::Line_2 Line_2;
		typedef Kernel::Line_3 Line_3;
		typedef Kernel::Iso_rectangle_2 Iso_rectangle_2;

		typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;
		typedef Kernel::Direction_2 Direction_2;
		typedef Kernel::Direction_3 Direction_3;
		typedef Kernel::Circle_2 Circle_2;
		typedef Kernel::Circle_3 Circle_3;
	};

}