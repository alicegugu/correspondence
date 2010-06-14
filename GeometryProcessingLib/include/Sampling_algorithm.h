#pragma once

#include "Another_CGAL_cast.h"
#include "Geometry_processing_algorithm.h"
#include "geodesic_algorithm_dijkstra.h"

/************************************************************************/
/* \brief             Sampling Algorithms for triangular meshes         */
/************************************************************************/

namespace Another
{
	class Sampling_algorithm
	{
	public:
		Sampling_algorithm(void);
		~Sampling_algorithm(void); 
		void FPS(vector<Vertex_handle>& features, Another_HDS_model& model, int sampleNumber );
	//	vector<Samples> sampling_gaussian_FPS(Another_HDS_model model, Vertex_PM_type);
	};
}