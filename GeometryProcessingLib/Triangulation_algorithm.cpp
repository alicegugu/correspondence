#include "Triangulation_algorithm.h"
#include <iostream>
#include <sstream>

namespace Another
{
// 
// 	/************************************************************************/
// 	/* /brief A abstract class for scanning a line of file                                             */
// 	/************************************************************************/
// 	class My_scanner
// 	{
// 	private:
// 		static vector<>
// 		public:
// 		virtual void scan(stringstream& line) = 0;
// 	}
// 
// 	/************************************************************************/
// 	/* /brief A scanner class for scanning geodesic path info               */
// 	/************************************************************************/
// 	class Geodesic_path_scanner:public My_scanner 
// 	{
// 		int m_vertex_from;
// 		int m_vertex_to;
// 	public:
// 		virtual void scan(stringstream& line) 
// 		{
// 			line >> m_vertex_from >> m_vertex_to;
// 		}
// 		int Get_vertex_from() { return m_vertex_from; }
// 		int Get_vertex_to() { return m_vertex_to; }
// 	}
// 
// 	/************************************************************************/
// 	/* /brief A scanner class for scanning vertex point info                */
// 	/************************************************************************/
// 	class Vertex_point_scanner : public My_scanner 
// 	{
// 
// 	public:
// 		virtual void scan(stringstream& line) 
// 		{
// 			string vert
// 		}
// 	}
// 
// 	/************************************************************************/
// 	/* /brief A scanner class for scanning edge point info                  */
// 	/************************************************************************/
// 	class Edge_point_scanner : public My_scanner 
// 	{
// 	public:
// 		virtual void scan(stringstream& line) 
// 		{
// 			string vert
// 		}
// 	}

	/************************************************************************/
	/* /brief read the text file                                             */
	/************************************************************************/
	bool Triangulation_algorithm::read_triangulation_file(std::string filename)
	{
		//Read the file content in memory
		if ( filename.c_str() != NULL )
		{
			std::ifstream stream;
			stream.open(filename.c_str());
			if (stream.is_open())
			{
				while(stream.good())
				{
					std::string s;  
					getline(stream,s); 
					std::stringstream oss;
					oss << s;
					std::string prefix;
					oss>>prefix;

					//A new path found
					if ( !prefix.compare("GeodesicPath") )
					{
						Geodesic_path path;
						oss >> path.m_vertex_index_from>> path.m_vertex_index_to;

						while( getline(stream,s))
						{
							std::stringstream ss2(s);
							ss2>>prefix;

							if ( !prefix.compare("VertexPoint"))
							{
								Geodesic_vertex_point vertexPoint;
								ss2>>vertexPoint.m_vertex_index;	
								path.m_path.push_back(vertexPoint);
							}
							if ( !prefix.compare("EdgePoint"))
							{
								Geodesic_edge_point edgePoint;
								ss2>>edgePoint.m_vertex_index_1>>edgePoint.m_vertex_index_2>>edgePoint.m_ratio;
								path.m_path.push_back(edgePoint);
							}
							if (!prefix.compare("GeodesicPath"))
							{

								int offset = s.size()+2;
								stream.seekg(-offset,std::ios_base::cur);
								break;
							}
						}
						m_geodesic_paths.push_back(path);
					}
				}
			}
		}
		
		//Init the triangulation data structure

		//Find the root face
		Geodesic_path pathTmp;
		Vertex_handle v1;
		Vertex_handle v2;
		Vertex_handle v3;
		
		vector<vector<int>> triangle_table;


		while(m_geodesic_paths.size()>0)
		{
			pathTmp = m_geodesic_paths[0];
			int v1index = pathTmp.m_vertex_index_from;
			int v2index = pathTmp.m_vertex_index_to;
			for (vector<vector<int>>::iterator it = triangle_table.begin(); it != triangle_table.end(); ++it)
			{
				if (it->size() == 4) //Skip integral triangle
				{
					continue;
				}
				else if (it->size() == 3 ) //Two edges has been detected for this triangle
				{
					int temp1 = it->at(2);
					int temp2 = it->at(0);
					if ( v1index == temp1 && v2index == temp2)
					{
						it->push_back(v2index);
					}
				}
				else if (it->size() == 2) //One edge has been detected for this triangle
				{
					int tmp = it->at(1);
					if( v1index == tmp ) //Find the second one
					{
						 it->push_back(v2index);
					}
				}
			}
		}
		return true;
	}

	/************************************************************************/
	/* /brief Triangulation based on features                               */
	/* /param[in] const Another_HDS_model& the model to be triangulation	*/
	/* /param[in] const vector<int>& features which are vertex index		*/
	/* /param[out] Another_triangulation_data_structure& the triangulation	*/
	/*				data structure											*/
	/************************************************************************/
	void Triangulation_algorithm::triangulate(const Another_HDS_model& model, const vector<int>& features, Another_triangulation_data_structure& tds)
	{
		string filename = "D:\\result.txt";
		read_triangulation_file(filename);
		//tds = m_tds;
	}

	/************************************************************************/
	/* /brief Segment the mesh according to the triangulation data structure*/
	/* /param[in] const Another_HDS_model& the model to be segmented		*/
	/* /param[in] Another_triangulation_data_structure& the triangulation	*/
	/* /param[out] vector<Another_HDS_model>& patches after segmentation	*/
	/************************************************************************/
	int Triangulation_algorithm::segment_mesh(const Another_HDS_model& model, Another_triangulation_data_structure& tds, vector<Another_HDS_model>& patches)
	{
		return -1;
	}

	Triangulation_algorithm::Triangulation_algorithm(void)
	{
	}

	Triangulation_algorithm::~Triangulation_algorithm(void)
	{
	}

}