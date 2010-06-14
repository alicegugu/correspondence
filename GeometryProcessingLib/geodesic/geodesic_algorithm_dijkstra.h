#pragma once

#include "geodesic_algorithm_graph_base.h"
#include "geodesic_mesh_elements.h"
#include "Another_HDS_model.h"
#include <vector>
#include <set>
#include <assert.h>

namespace Another
{
	class DijkstraNode
	{
		typedef DijkstraNode* node_pointer;
	public: 
		DijkstraNode(){};
		~DijkstraNode(){};

		double& distance_from_source(){return m_distance;};
		node_pointer& previous(){return m_previous;};
		unsigned& source_index(){return m_source_index;};
		Vertex_handle vertex(){return m_vertex;};
		void set_vertex(Vertex_handle v) {m_vertex = v;}

		void clear()
		{
			m_distance = GEODESIC_INF;
			m_previous = NULL;
		}

		bool operator()(node_pointer const s1, node_pointer const s2) const
		{
			return s1->distance_from_source() != s2->distance_from_source() ?
				s1->distance_from_source() < s2->distance_from_source() :
			s1->vertex()->GetIndex() < s2->vertex()->GetIndex();
		};

		double distance(Surface_point* p)
		{
			return m_vertex->distance(p->get_point());
		}

		Surface_point surface_point()
		{
			return Surface_point(m_vertex);
		}

	private: 
		double m_distance;					//distance to the closest source
		unsigned m_source_index;			//closest source index
		node_pointer m_previous;			//previous node in the geodesic path
		Vertex_handle m_vertex;			//correspoding vertex
	};

	class GeodesicAlgorithmDijkstra: public GeodesicAlgorithmGraphBase<DijkstraNode>
	{
	public:
		typedef DijkstraNode Node;
		typedef Node* node_pointer;

		GeodesicAlgorithmDijkstra(Another_HDS_model* mesh):
		GeodesicAlgorithmGraphBase<Node>(mesh)
		{
			m_type = DIJKSTRA;
			m_nodes.resize(mesh->size_of_vertices());
			int i = 0;
			for (Vertex_iterator v = mesh->vertices_begin(); v!= mesh->vertices_end(); ++v)
			{
				m_nodes[i].set_vertex(v);
				v->m_node_index = i;
				i++;
			}
		};

		~GeodesicAlgorithmDijkstra(){};
		void assgin_distance();

	protected:

		void list_nodes_visible_from_source(Mesh_element_base* p, 
			std::vector<node_pointer>& storage);		//list all nodes that belong to this mesh element

		void list_nodes_visible_from_node(node_pointer node,			//list all nodes that belong to this mesh element
			std::vector<node_pointer>& storage,
			std::vector<double>& distances, 
			double threshold_distance);	//list only the nodes whose current distance is larger than the threshold
	};

}