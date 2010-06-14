#include ".\include\Sampling_algorithm.h"

namespace Another
{
	bool compare_function (Another::Surface_point i,Another::Surface_point j) 
	{
		Vertex_handle vi = Another_CGAL_cast::vertex_cast(i.base_element());
		Vertex_handle vj = Another_CGAL_cast::vertex_cast(j.base_element());
		return (vi->m_distance < vj->m_distance ); 
	}

	Sampling_algorithm::Sampling_algorithm(void)
	{
	}

	Sampling_algorithm::~Sampling_algorithm(void)
	{
	}

	/************************************************************************/
	/* /brief Farthest Point Sampling                                       */
	/* /param[in and out] vector<Vertex_handle>& features: Store source		*/
	/* points and output samples											*/
	/* /param[in] Another_HDS_model& model									*/
	/* /param[in] int desired sample number									*/
	/************************************************************************/
	void Sampling_algorithm::FPS(vector<Vertex_handle>& features, Another_HDS_model& model, int sampleNumber )
	{
		//Assert the number of desired samples is larger than source points
		if (features.size()>sampleNumber)
		{
			cout<<"sources point number is larger than sample number"<<endl;
			return;
		}

		//Mark source and targets points on model, using a property map
		//Source marked 1 and target marked -1
		vector<Surface_point> sources;
		Vertex2int_map_type vertex2props;
		Vertex_PM_type vpm(vertex2props);

		Vertex_iterator vitb = model.vertices_begin(); 
		Vertex_iterator vite = model.vertices_end();
		CGAL_For_all(vitb, vite) put(vpm, &(*vitb) , -1);

		for (size_t i = 0; i<features.size(); i++)
		{
			Vertex_handle vertexP = features[i];
			Surface_point surfaceP(vertexP);
			sources.push_back(surfaceP);
			put(vpm, &(*vertexP) , 1);
		}

		vector<Another::Surface_point> targets;
		for (Another::Vertex_iterator i = model.vertices_begin(); i!= model.vertices_end(); ++i )
		{
			Another::Vertex_handle tmpv = i;
			int tmpi = get(vpm, &(*i));
			if ( tmpi == -1)
			{
				targets.push_back( Surface_point(tmpv) );
			}
		}

		// recursive voronoi diagram construction
		int currentSampleNumber = sources.size();
		
		model.calculate_edge_length();
		while( currentSampleNumber < sampleNumber )
		{
			GeodesicAlgorithmDijkstra geodesicAlgorithm(&model);
			geodesicAlgorithm.propagate(sources);
			//geodesicAlgorithm.print_statistics();
			// 				std::vector<Another::Surface_point> path;
			// 				for(int i=0; i<targets.size(); ++i)
			// 				{
			// 					algorithm->trace_back(targets[i], path);
			// 				}

			//Get the largest distance point and add it into source, remove it from target
			geodesicAlgorithm.assgin_distance();
			sort (targets.begin(), targets.end(), compare_function); 

			sources.push_back( targets[targets.size()-1] );
			currentSampleNumber++;
			targets.pop_back();
		}

		//send back samples
		features.clear();
		for (size_t i = 0; i<sources.size(); i++)
		{
			Surface_point tmpSurfacePoint = sources[i];
			Vertex_handle tmpVertexHandle = Another_CGAL_cast::vertex_cast(tmpSurfacePoint.base_element());
			features.push_back(tmpVertexHandle);
		}
	}

}
