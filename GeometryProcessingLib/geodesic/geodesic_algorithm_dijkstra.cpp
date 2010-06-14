#include "geodesic_algorithm_dijkstra.h"

namespace Another
{

	void GeodesicAlgorithmDijkstra::list_nodes_visible_from_node(node_pointer node, //list all nodes that belong to this mesh element
		std::vector<node_pointer>& storage,
		std::vector<double>& distances,
		double threshold_distance)
	{
		Vertex_handle v = node->vertex();
		assert(storage.size() == distances.size());

		bool first = true;
		Halfedge_handle firste = v->halfedge();
		Halfedge_handle tmpe = firste;
		int i = 0;
		while( tmpe != firste || first )
		{
			first = false;
			Vertex_handle new_v = tmpe->opposite_vertex(v);
			node_pointer new_node = &m_nodes[new_v->m_node_index];
			if(new_node->distance_from_source() > threshold_distance + tmpe->get_length())
			{
				storage.push_back(new_node);
				distances.push_back(tmpe->get_length());
			}

			tmpe = tmpe->next();
			tmpe = tmpe->opposite();
			i++;
			if (i > 100)
			{
				cerr<<"infinite loop for 1-ring neighor"<<endl;
			}
		}
	}


	void GeodesicAlgorithmDijkstra::list_nodes_visible_from_source(Base_handle p, std::vector<node_pointer>& storage)
	{
		assert(p->type() != UNDEFINED_POINT_TYPE);

		if(p->type() == FACE_TYPE)
		{
			Face_handle f = Another_CGAL_cast::face_cast(p);
			for(unsigned i=0; i<3; ++i)
			{
				vector<Vertex_handle> vs; 
				f->get_all_vertex(vs);
				for ( unsigned int i = 0; i< vs.size(); i++ )
				{
					Vertex_handle tmpv = vs[i];
					storage.push_back(&m_nodes[tmpv->m_node_index]);
				}
			}
		}
		else if(p->type() == EDGE_TYPE)
		{
			Halfedge_handle e = Another_CGAL_cast::halfedge_cast(p);
			Vertex_handle tmpv = e->vertex();
			storage.push_back(&m_nodes[tmpv->m_node_index]);
			tmpv = e->opposite_vertex();
			storage.push_back(&m_nodes[tmpv->m_node_index]);
		}
		else			//VERTEX
		{
			Vertex_handle v = Another_CGAL_cast::vertex_cast(p);
			storage.push_back(&m_nodes[v->m_node_index]);
		}
	}

	void GeodesicAlgorithmDijkstra::assgin_distance()
	{
		for (size_t i = 0; i<m_nodes.size(); i++)
		{
			Another::DijkstraNode node = m_nodes[i];
			Vertex_handle v = node.vertex();
			v->m_distance = node.distance_from_source();
		}
	}


}