#pragma once
#include <CGAL/Triangulation_vertex_base_2.h>

namespace Another
{
	
	//class to hold the vertex information of the triangulation
	//This vertex can be a surface point on the mesh
	//This vertex is a model of concept of Triangulation_vertex_base_2
	class Another_triangulation_vertex
	{
	public:
		Another_triangulation_vertex(void);
		~Another_triangulation_vertex(void);
	};

}