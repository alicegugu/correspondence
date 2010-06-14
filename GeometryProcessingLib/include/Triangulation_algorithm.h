#pragma once

#include <vector>;
#include <string>
#include "Another_triangulation_data_structure.h"

namespace Another
{
	struct Geodesic_point{	};

	struct Geodesic_vertex_point:Geodesic_point
	{
		int m_vertex_index;
	};

	struct Geodesic_edge_point:Geodesic_point
	{
		int m_vertex_index_1;
		int m_vertex_index_2;
		float m_ratio;
	};
	struct Geodesic_triangle
	{
		int m_vertex_index1;
		int m_vertex_index2;
		int m_vertex_index3;
	};
	struct Geodesic_path
	{
		int m_vertex_index_from;
		int m_vertex_index_to;
		std::vector<Geodesic_point> m_path;
	};

	/************************************************************************/
	/* /brief Utility class for triangulation algorithm                     */
	/************************************************************************/
	class Triangulation_algorithm
	{
		std::vector<Geodesic_path> m_geodesic_paths;
		std::vector<Geodesic_point> m_vertex;
		std::vector<Geodesic_triangle> m_triangles;
		Another_triangulation_data_structure m_tds;
		typedef Another_triangulation_data_structure::Vertex_handle Vertex_handle;
		typedef Another_triangulation_data_structure::Face_handle Face_handle;
	public:
		Triangulation_algorithm(void);
		~Triangulation_algorithm(void);

		//Read the edges of triangulations from shiqing program's result
		bool read_triangulation_file(std::string filename);
		//Triangulate the model using feature points
		void triangulate(const Another_HDS_model& model, const vector<int>& features, Another_triangulation_data_structure& tds);
		//Segment the mesh according to the triangulation data structure
		int segment_mesh(const Another_HDS_model& model, Another_triangulation_data_structure& tds, vector<Another_HDS_model>& patches);
	};
}