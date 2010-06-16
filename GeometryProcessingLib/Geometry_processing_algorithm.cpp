#include ".\include\Geometry_processing_algorithm.h"
#include ".\include\Polyhedral_surface_rings.h"
#include ".\include\Polyhedral_surface_operations.h"
#include ".\include\Polyhedral_surface.h"

#include <vector>
#include <math.h>

const int MUTUALK = 50;
const double eps = 0.00000001;
namespace Another
{
	GeometryProcessingAlgorithm::GeometryProcessingAlgorithm(void)
	{
	}

	GeometryProcessingAlgorithm::~GeometryProcessingAlgorithm(void)
	{
	}

	double GeometryProcessingAlgorithm::fix_sine(double sine)
	{
		if(sine >= 1)
			return 1;
		else if(sine <= -1)
			return -1;
		else
			return sine;
	}
	double GeometryProcessingAlgorithm::compute_angle_rad(Another::Point P, Another::Point Q, Another::Point R)
	{
		double PI_s = 3.14159265359;

		Another::Vector_3 u = P - Q;
		Another::Vector_3 v = R - Q;

		// check
		double product = std::sqrt(u*u) * std::sqrt(v*v);
		if(product == 0)
			return 0.0;

		// cosine
		double dot = (u*v);
		double cosine = dot / product;

		// sine
		Another::Vector_3 w = CGAL::cross_product(u,v);
		double AbsSine = std::sqrt(w*w) / product;

		if(cosine >= 0)
			return std::asin(fix_sine(AbsSine));
		else
			return PI_s-std::asin(fix_sine(AbsSine));
	}

	double GeometryProcessingAlgorithm::computeCot(Another::Point position_v_i, Another::Point position_v_j, Another::Point position_v_k, Another::Point position_v_l )
	{
		// Compute the norm of v_j -> v_i vector
		Another::Vector_3 edge = position_v_i - position_v_j;
		double len = std::sqrt(edge*edge);

		// Compute angle of (v_j,v_i,v_k) corner (i.e. angle of v_i corner)
		// if v_k is the vertex before v_j when circulating around v_i

		double gamma_ij  = compute_angle_rad(position_v_i, position_v_k, position_v_j);


		// Compute angle of (v_l,v_i,v_j) corner (i.e. angle of v_i corner)
		// if v_l is the vertex after v_j when circulating around v_i

		double delta_ij = compute_angle_rad(position_v_i, position_v_l, position_v_j);


		double weight = 0.0;
		assert(len != 0.0);    // two points are identical!
		double cot_gamma = 0.0;
		double cot_delta = 0.0;
		double PI_s = 3.14159265359;
		if (gamma_ij<0.01)
		{
			cot_gamma = 10;
		}
		else if (gamma_ij>PI_s -0.01)
		{
			cot_gamma = -10;
		}
		else
		{
			cot_gamma = (double)1.0/std::tan(gamma_ij);
		}
		if (delta_ij<0.01)
		{
			cot_delta = 10;
		}
		else if(delta_ij> (PI_s-0.01))
		{
			cot_delta = -10;
		}
		else
		{
			cot_delta = (double)1.0/std::tan(delta_ij);
		}

		if(len != 0.0)
			weight = cot_gamma +cot_delta;
		else
			weight = 0;


		return weight;
	}



	/************************************************************************/
	/* \brief Calculate the Laplacian matrix of the model                   */
	/* \param																*/
	/************************************************************************/
	double* GeometryProcessingAlgorithm::LaplacianMatrix(Another_HDS_model&model, LapalacianOperatorType type)
	{
		double* lpMatrix = NULL;
		int n = model.size_of_vertices();
		lpMatrix = (double*)malloc(n*n*sizeof(double));
		
		
		for (int i=0; i<n; ++i)
		{
			for (int j=0; j<n; ++j)
			{
				lpMatrix[i*n+j] =0;
			}
		}

		//For each vertex
		Another::Vertex_handle tmpv = model.vertices_begin();
		for (;tmpv!= model.vertices_end(); ++tmpv)
		{
			int tmpIndex = tmpv->GetIndex()-1;
			Another::Point  Vi_position = tmpv->point();

			//1-ring neighbor
			Another::Halfedge_handle tmpEdge;
			Another::Halfedge_handle firstEdge = tmpv->halfedge();
			tmpEdge = firstEdge;
			bool first = true;
			double sum = 0.0;
			while(tmpEdge != firstEdge||first) //1-ring neighbor
			{
				first = false;

				Another::Vertex_handle vtmpj; //Vertex_handle vtmp;

				Another::Halfedge_handle lastEdge = tmpEdge;
				Another::Halfedge_handle lastOppositeEdge = tmpEdge->opposite();
				vtmpj =  lastOppositeEdge->vertex(); //neighbor
				Another::Point Vj_position = vtmpj->point();

				Another::Halfedge_handle Eik = tmpEdge->next();
				Another::Vertex_handle Vk = Eik->vertex();
				Another::Point  Vk_position = Vk->point();

				Another::Halfedge_handle Ejl = lastOppositeEdge->next();
				Another::Vertex_handle Vl = Ejl->vertex();
				Another::Point  Vl_position = Vl->point();

				double wij = 0.0;
				if ( type == Another::LapalacianOperatorType::COTANGENT )
				{
					wij = computeCot(Vi_position, Vj_position, Vk_position, Vl_position);
				}
				else if ( type ==  Another::LapalacianOperatorType::UNIFORM )
				{
					wij = 1.0;
				}
				int nIndex = vtmpj->GetIndex()-1;
				lpMatrix[tmpIndex*n+nIndex] = -wij;
				sum = sum + wij;
				tmpEdge = Ejl->next();
			
			}
			lpMatrix[tmpIndex*n+tmpIndex] = sum;
		}
		return lpMatrix;
	}

	/************************************************************************/
	/* \brief Gather points around the vertex using one-neighhor-ring.		*/
	/*  copied from CGAL/examples/Jet_fitting_3/Mesh_estimation.cpp			*/
	/*	@param[in] v the vertex												*/
	/*	@param[out] neighbor the neighbors vertexes of v					*/
	/*	@param[in] The property map containing the							*/
	/************************************************************************/
	void GeometryProcessingAlgorithm::GatheringNeighborRings(Vertex* v, Vertex_PM_type& vpm, std::vector<Point_3>& neighhor, int neighborPointsToUse, int min_nb_points, int nb_rings)
	{
		std::vector<Vertex*> gathered;
		neighhor.clear();
		//OPTION -p nb_points_to_use, with nb_points_to_use != 0. Collect
		//enough rings and discard some points of the last collected ring to
		//get the exact "nb_points_to_use"
		if ( neighborPointsToUse != 0 ) 
		{
			Poly_rings::collect_enough_rings(v, neighborPointsToUse, gathered, vpm);
			if ( gathered.size() > neighborPointsToUse ) gathered.resize(neighborPointsToUse);
		}
		else 
		{
			// nb_points_to_use=0, this is the default and the option -p is not considered;
			// then option -a nb_rings is checked. If nb_rings=0, collect
			// enough rings to get the min_nb_points required for the fitting
			// else collect the nb_rings required
			if ( nb_rings == 0 )
				Poly_rings::collect_enough_rings(v, min_nb_points, gathered, vpm);
			else 
				Poly_rings::collect_i_rings(v, nb_rings, gathered, vpm);
		}

		//store the gathered points
		std::vector<Vertex*>::iterator itb = gathered.begin();
		std::vector<Vertex*>::iterator ite = gathered.end();
		CGAL_For_all(itb,ite) neighhor.push_back((*itb)->point());

	}

	Vertex2GuassianCurvature_PM_type& GeometryProcessingAlgorithm::GaussianCurvature(Vertex2Monge_PM_type& monges, Another_HDS_model& model, Vertex2GuassianCurvature_PM_type& GCpm)
	{
		Vertex_iterator vitb;
		Vertex_iterator vite;
		vitb = model.vertices_begin(); vite = model.vertices_end();
		for (; vitb != vite; vitb++)
		{
			double guassianCurvature;
			double k1, k2;
			Another_Monge_form monge = get(monges, vitb);
			k1 = monge.coefficients()[0];
			k2 = monge.coefficients()[1];
			guassianCurvature = k1*k2;
			put(GCpm,vitb,guassianCurvature);
			vitb->SetGaussianCurvature(guassianCurvature);
		}
		return GCpm;
	}

	Vertex2Monge_PM_type& GeometryProcessingAlgorithm::Monge_jet_fitting(Another_HDS_model& model, int min_nb_points, int d_fitting, int d_monge, int nb_rings, int neighborPointsToUse, Vertex2Monge_PM_type&monges)
	{
		//perform calculation for each vertex
		Vertex2int_map_type vertex2props;
		Vertex_PM_type vpm(vertex2props);


		//Hedge, with enriched hedge
		//HEdgePM_type hepm = get_hepm(boost::edge_weight_t(), P);
		//Hedge, using a std::map
		Hedge2double_map_type hedge2props;
		Hedge_PM_type hepm(hedge2props);

		//Facet PM, with enriched Facet
		//FacetPM_type fpm = get_fpm(boost::vertex_attribute_t(), P);
		//Facet PM, with std::map
		Face2normal_map_type facet2props;
		Face_PM_type fpm(facet2props);

		//initialize Polyhedral data : length of edges, normal of facets
		Hedge_ops::compute_edges_length(model, hepm);
		Face_ops::compute_facets_normals(model, fpm);

		std::vector<Point_3> neighbors;
		Vertex_iterator vitb;
		Vertex_iterator vite;
		vitb = model.vertices_begin(); 
		vite = model.vertices_end();
		CGAL_For_all(vitb, vite) put(vpm, &(*vitb) , -1);

		vitb = model.vertices_begin(); vite = model.vertices_end();
		for (; vitb != vite; vitb++) 
		{
			//initialize
			Vertex* v = &(*vitb);
			neighbors.clear();
			Another_Monge_form monge_form;

			//gather points around the vertex using rings
			GatheringNeighborRings(v, vpm,neighbors, neighborPointsToUse, min_nb_points, nb_rings);

			//skip if the nb of points is to small
			if ( neighbors.size() < min_nb_points )
			{std::cerr << "not enough pts for fitting this vertex" << neighbors.size() << std::endl;
			continue;}

			// perform the fitting
			Another_Monge_via_jet_fitting monge_fit;
			monge_form = monge_fit(neighbors.begin(), neighbors.end(),
				d_fitting, d_monge);
			//switch min-max ppal curv/dir wrt the mesh orientation
			const Normal normal_mesh = Face_ops::compute_vertex_average_unit_normal(v, fpm);
			monge_form.comply_wrt_given_normal(normal_mesh);

			put(monges, vitb->halfedge()->vertex(), monge_form);
		}
		return monges;
	}

	/************************************************************************/
	/* /brief Find the local maxim and minim of Gaussian curvature field    */
	/* /param[in] Another_HDS_model& model									*/
	/* /param[out] Local maxim and minim feature point list					*/
	/************************************************************************/
	void GeometryProcessingAlgorithm::LocalMaxandMinGuassian(Another_HDS_model& model, vector<Vertex_handle>& featurePoints)
	{
		//Init a property map and set it with default value -1
		Vertex2int_map_type vertex2props;
		Vertex_PM_type vpm(vertex2props);
		bool clear = ClearPropertyMap(model,vpm,-1);

		//For each vertex, check the neighbor's values to determine whether it is a maxim or minim
		Vertex_iterator vit;
		Vertex_iterator vitb;
		Vertex_iterator vite;
		vitb = model.vertices_begin(); 
		vite = model.vertices_end();

		for ( vit = vitb; vit != vite; ++vit)
		{
			//Gather 1-ring neighbor
			std::vector<Vertex*> gathered;
			Poly_rings::collect_i_rings( &(*vit), 1, gathered, vpm); 
			assert(gathered.size());

			double gc = vit->GetGaussianCurvature();
			bool isFeature = true;
			
			Vertex* lastv = gathered[1];
			double lastGc = lastv->GetGaussianCurvature();			
			double lastDeltaGc = gc - lastGc;

			for (size_t i = 2; i<gathered.size(); i++)
			{
				Vertex* tmpv = gathered[i];
				double tmpGc = tmpv->GetGaussianCurvature();		
				double deltaGc = gc - tmpGc;
				if (deltaGc*lastDeltaGc>0.1 ) //lastgc and gc has same sign
				{
					lastDeltaGc = deltaGc;
					continue;
				}
				else
				{
					isFeature = false;
					break;
				}
			}
			if (isFeature)
			{
				Vertex_handle tmpv = vit;
				featurePoints.push_back(tmpv);
			}
		}
	}


	/************************************************************************/
	/* /brief Find the local maxim and minim of Gaussian curvature field    */
	/* /param[in] Another_HDS_model& model									*/
	/* /param[in] The minim distance between feature points					*/
	/* /param[out] Local maxim and minim feature point list					*/
	/************************************************************************/
	void GeometryProcessingAlgorithm::LocalMaxandMinGuassian(Another_HDS_model& model, float radius, vector<Vertex_handle>& featurePoints)
	{
		//Init a property map for gathering neighbor rings and set it with default value -1
		Vertex2int_map_type vertex2props;
		Vertex_PM_type vpm(vertex2props);
		ClearPropertyMap(model,vpm,-1);

		//Init a property map for record feature point allocation and set it with default value 0
		//0 signifies the vertex is available as feature point
		//1 signifies the vertex is not available as feature point
		Vertex2int_map_type vertex2propsFeaturePoint;
		Vertex_PM_type vpmFeaturePoint(vertex2propsFeaturePoint);
		ClearPropertyMap(model,vpmFeaturePoint,0);

		//For each vertex, check the neighbor's values to determine whether it is a maxim or minim
		Vertex_iterator vit;
		Vertex_iterator vitb;
		Vertex_iterator vite;
		vitb = model.vertices_begin(); 
		vite = model.vertices_end();

		for ( vit = vitb; vit != vite; ++vit)
		{
			//Check whether this vertex is available as feature point
			int signify = get(vpmFeaturePoint,vit);
			if (signify)
			{
				continue;
			}

			//Gather 1-ring neighbor
			std::vector<Vertex*> gathered;
			Poly_rings::collect_i_rings( &(*vit), 1, gathered, vpm); 
			assert(gathered.size());

			double gc = vit->GetGaussianCurvature();
			bool isFeature = true;

			Vertex* lastv = gathered[1];
			double lastGc = lastv->GetGaussianCurvature();			
			double lastDeltaGc = gc - lastGc;

			for (size_t i = 2; i<gathered.size(); i++)
			{
				Vertex* tmpv = gathered[i];
				double tmpGc = tmpv->GetGaussianCurvature();		
				double deltaGc = gc - tmpGc;
				if (deltaGc*lastDeltaGc>0.1 ) //lastgc and gc has same sign
				{
					lastDeltaGc = deltaGc;
					continue;
				}
				else
				{
					isFeature = false;
					break;
				}
			}
			if (isFeature)
			{
				Vertex_handle tmpv = vit;
				featurePoints.push_back(tmpv);

				//All the points near the new feature point within the radius cannot be set as new feature point
				//So we set the signify of availability for 1 rings neighbor of the new feature point as 1
				std::vector<Vertex*> affectedNeighbors;
				int rad = radius;
				Poly_rings::collect_i_rings( &(*tmpv), rad, affectedNeighbors, vpm);
				for (size_t i = 0; i< affectedNeighbors.size(); i++)
				{
					Vertex* neighborv = affectedNeighbors[i];
					put(vpmFeaturePoint, neighborv , 1);
				}
			}
		}
	}

	/************************************************************************/
	/* /brief Set a property type to default value                          */
	/* /param[in] Another_HDS_model& model									*/
	/* /param[in] Vertex_PM_type& property map								*/
	/* /param[in] int default value											*/
	/* /param[out] bool success or not										*/
	/************************************************************************/
	bool GeometryProcessingAlgorithm::ClearPropertyMap(Another_HDS_model& model, Vertex_PM_type& vpm, float defaultValue)
	{
		Vertex_iterator vitb;
		Vertex_iterator vite;
		vitb = model.vertices_begin(); 
		vite = model.vertices_end();
		CGAL_For_all(vitb, vite) put(vpm, &(*vitb) , defaultValue);
		return true;
	}

	/************************************************************************/
	/* /brief embed source and target model in planar domain                */
	/* /param source[in] source model index									*/
	/* /param target[in] target model index									*/
	/************************************************************************/
	void GeometryProcessingAlgorithm::Embed(Another_HDS_model& source, Another_HDS_model& target)
	{
		Another::GeometryProcessingAlgorithm geoAl;
		vector<Vertex_handle>& sourcePoints = source.get_samples();
		vector<Vertex_handle>& targetPoints = target.get_samples();
		if (sourcePoints.size() == 0 && targetPoints.size()==0)
		{
			cout <<"Error GeometryProcessingAlgorithm::Embed: No source and target feature points"<<endl;
			return;
		}
		vector<complex<double>> sourceMidEdgeSample;
		vector<complex<double>> targetMidEdgeSample;
		geoAl.ProjectSamples(sourcePoints,sourceMidEdgeSample);
		geoAl.ProjectSamples(targetPoints,targetMidEdgeSample);
		source.set_mid_edge_features(sourceMidEdgeSample);
		target.set_mid_edge_features(targetMidEdgeSample);
	}

	/************************************************************************/
	/* /brief Mobius Voting on two model                                    */
	/* /prama [in] source source model										*/
	/* /prama [in] target target model										*/
	/* /return correspondence matrix										*/
	/************************************************************************/
	double* GeometryProcessingAlgorithm::MobiusVoting(Another_HDS_model& source, Another_HDS_model& target)
	{
		source.calculate_edge_length();
		target.calculate_edge_length();

		vector<Vertex_handle>& sourcePoints = source.get_samples();
		vector<Vertex_handle>& targetPoints = target.get_samples();

		vector<complex<double>> sourceMidEdgeSample;
		vector<complex<double>> targetMidEdgeSample;

		ProjectSamples(sourcePoints,sourceMidEdgeSample);
		ProjectSamples(targetPoints,targetMidEdgeSample);
		
		if (sourcePoints.size() < 10 || targetPoints.size() <10)
		{
			cerr<<"MobiusVoting Error: sample size is too small"<<endl;
			return NULL;
		}
		if (sourcePoints.size() != targetPoints.size())
		{
			cerr <<"MobiusVoting Error: sample size is not equal"<<endl;
			return NULL;
		}

		int iteration = 0;
		int samplesize = sourcePoints.size();
		int interationNum = samplesize*(samplesize-1);
		double* correspondenceMatrix = (double*)malloc(sizeof(double)*samplesize*samplesize);
		for (size_t i = 0; i< samplesize*samplesize;i++)
		{
			correspondenceMatrix[i] = 0;
		}
		while( iteration < interationNum )
		{
			//Random z1, z2, z3 in source samples
			int z[3] = { -1, -1, -1};
			srand ( time(NULL) );
			for (int i = 0; i< 3;)
			{
				int ramdom = rand()%samplesize;	
				if (ramdom == z[0] || ramdom == z[1] || ramdom == z[2])
				{
					continue;
				}
				z[i] = ramdom;
				i++;
			}
			
			vector<complex<double>> sourceTriplet;
			for (int i = 0; i<3; i++)
			{
				  sourceTriplet.push_back(sourceMidEdgeSample[z[i]]);
			}

			//Random w1, w2, w3 in target samples
			int w[3] = { -1, -1, -1};
			srand ( time(NULL)+ 1 );
			for (int i = 0; i< 3;)
			{
				int ramdom = rand()%samplesize;	
				if (ramdom == w[0] || ramdom == w[1] || ramdom == w[2])
				{
					continue;
				}
				w[i] = ramdom;
				i++;
			}

			vector<complex<double>> targetTriplet;
			for (int i = 0; i<3; i++)
			{
				targetTriplet.push_back(targetMidEdgeSample[w[i]]);
			}

			//Find the mobius transfromation
			vector<complex<double>> m1, m2;
			GetMobiusTransformation(sourceTriplet,m1);
			GetMobiusTransformation(targetTriplet,m2);

			vector<pair<int, int>> pairs;
			pairs.clear();
			//Todo FindMutualNeightbor
			FindMutualNeightbor(sourceMidEdgeSample,m1, targetMidEdgeSample,m2, pairs);
			int n = sourceMidEdgeSample.size();
			if ( pairs.size()>= n/2 )
			{
				//To do Measuring deformation error
				double energy = CalculateDeformationEnergy(sourceMidEdgeSample,m1, targetMidEdgeSample,m2, pairs);
				for (size_t i = 0; i < pairs.size(); i++)
				{
					int k = pairs[i].first;
					int l = pairs[i].second;
					double vvv = 1.0/(eps+energy/((double)pairs.size()));
					double oldvalue= correspondenceMatrix[k*samplesize+l] ;
					double aaaa = correspondenceMatrix[k*samplesize+l] + vvv;
					correspondenceMatrix[k*samplesize+l] = aaaa;
					int asdf = 0;
				}
			}

			iteration++;

		}
		FILE* refile;	
		refile = fopen("D:\\matrixCorrespoendence.txt","w+");
		for (int i = 0; i<samplesize;i++)
		{
			for (int k = 0; k<samplesize; k++)
			{
				fprintf(refile,"%f ",correspondenceMatrix[i*samplesize+k]);
			}
			fprintf(refile,"\n");
			
		}
		fclose(refile);
		return correspondenceMatrix;
	}

	/************************************************************************/
	/* /brief Find the mobius transformation which transform the point z1,	*/
	/*	      z2, z3 to y_j = e{i 2*pi*j/3} j = 1,2,3                       */
	/*	/param input point													*/
	/************************************************************************/
	void GeometryProcessingAlgorithm::GetMobiusTransformation(vector<complex<double>>& z, vector<complex<double>>& transformation)
	{
		assert(z.size() == 3);
		complex<double> z11 = z[1] -z[2];
		complex<double> z12 = z[0]*z[2]-z[0]*z[1];
		complex<double> z21 = z[1] - z[0];
		complex<double> z22 = z[0]*z[2] - z[2]*z[1];
		transformation.clear();
		
		complex<double> y11(-1.0/3.0,0.0);
		complex<double> y12(1.0/6.0,sqrt(3.0)/6.0);
		complex<double> y21(-1.0/3.0,0.0);
		complex<double> y22(1.0/6.0,-sqrt(3.0)/6.0);

		complex<double> a;
		a = y11*z11 + y12*z21;
		complex<double> b;
		b = y11*z12 + y12*z22;
		complex<double> c;
		c = y21*z11 + y22*z21;
		complex<double> d;
		d = y21*z12 + y22*z22;

		transformation.clear();
		transformation.push_back(a);
		transformation.push_back(b);
		transformation.push_back(c);
		transformation.push_back(d);
	}

	/************************************************************************/
	/* /brief Project point in Points to the nearest mid-edge vertex        */
	/* /param [in]Points input point										*/
	/* /param [out]MidEdgeSample output mid-edge vertex						*/
	/************************************************************************/
	void GeometryProcessingAlgorithm::ProjectSamples(vector<Vertex_handle>& Points,vector<complex<double>>& MidEdgeSample)
	{
		if ( Points.size()==0 )
		{
			cout<<"Error GeometryProcessingAlgorithm::ProjectSamples: No feature points for project."<<endl;
			return;
		}
		for (size_t i = 0; i<Points.size(); i++)
		{
			Vertex_handle v = Points[i];
			Halfedge_handle firstEdge = v->halfedge();
			Halfedge_handle tmpEdge = firstEdge;
			bool first = true;
			double len = tmpEdge->get_length();
			double tmplen;
			complex<double> closetMidEdgeP(tmpEdge->GetU(),tmpEdge->GetConjU() );
			while(tmpEdge != firstEdge||first)
			{
				tmplen = tmpEdge->get_length();
				if (tmplen > len&&tmpEdge->visited == false )
				{
					complex<double> midEdgeP(tmpEdge->GetU(),tmpEdge->GetConjU() );
					closetMidEdgeP = midEdgeP;
					tmpEdge->visited = true;
				}	

				first = false;
				tmpEdge = tmpEdge->opposite();
				tmpEdge = tmpEdge->next();
			}
			MidEdgeSample.push_back(closetMidEdgeP);
		}
	}


	/************************************************************************/
	/* /brief Find mutually nearst neighors to formulate candidate			*/
	/*		correspodence                                                   */
	/* /param [sourceMidEdgeSample]
	/************************************************************************/
	void GeometryProcessingAlgorithm::FindMutualNeightbor(vector<complex<double>>& sourceMidEdgeSample,const vector<complex<double>>& m1, vector<complex<double>>& targetMidEdgeSample, const vector<complex<double>>&m2,vector<pair<int, int>>& pairs)
	{
		vector<complex<double>> transformedSource;
		vector<complex<double>> transformedTarget;
		for (size_t i = 0 ; i< sourceMidEdgeSample.size(); i++)
		{
			complex<double> sampleSource = sourceMidEdgeSample[i];
			complex<double> sampleTarget = targetMidEdgeSample[i];
			complex<double> transformedSampleSource;
			complex<double> transformedSampleTarget;
			MobiusTransform(sampleSource, m1, transformedSampleSource);
			MobiusTransform(sampleTarget, m2, transformedSampleTarget);
			transformedSource.push_back(transformedSampleSource);
			transformedTarget.push_back(transformedSampleTarget);
		}

		complex<double> infiniteComplex(INFINITE,INFINITE);
// 		for (size_t i = 0; i<transformedSource.size(); i++)
// 		{
// 			complex<double> sampleSource = transformedSource[i];
// 			int indexTarget = MatchClosest(transformedTarget, sampleSource);
// 			if (indexTarget >= 0)
// 			{
// 				pair<int, int> pa(i,indexTarget);
// 				pairs.push_back(pa);
// 				transformedTarget[indexTarget] = infiniteComplex;
// 			}
// 		}

		//caculate the distance between points
		vector<vector<double>> mutralMatrix;
		mutralMatrix.resize(transformedSource.size());
		for (size_t i = 0; i<transformedSource.size(); i++)
		{
			mutralMatrix[i].resize(transformedSource.size() -i );
			for (size_t k = i; k<transformedTarget.size(); k++)
			{
				if (i == k)
				{
					mutralMatrix[i][k-i]= 1000000.0;
				}
				else
				{
					complex<double> sourceTmp = transformedSource[i];
					complex<double> targetTmp = transformedTarget[k];
					double realPart = abs(sourceTmp.real() - targetTmp.real());
					double imagePart = abs(sourceTmp.imag() - targetTmp.imag());
					double dis = sqrt( realPart*realPart + imagePart*imagePart );
					if (dis < 0.000000001)
					{
						dis = 0.0;
					}
					mutralMatrix[i][k-i] = dis;
				}
			}
		}

		double smallest = 100.0;
		double tmpEntry;
		pair<int,int> tmpPair;
		do
		{
			smallest = 1000000.0;
			for ( int k = 0; k<mutralMatrix.size(); k++ )
			{
				for ( int l = 0; l<mutralMatrix[k].size(); l++ )
				{
					tmpEntry = mutralMatrix[k][l];
					if (tmpEntry < smallest)
					{
						tmpPair.first = k;
						tmpPair.second = l+k;
						smallest = tmpEntry;
					}
				}

			}
			if (smallest < 1000000.0) //if find a entry is larger than zero and is the largest
			{
				pairs.push_back(tmpPair);

				for (size_t i = 1; i<mutralMatrix[tmpPair.first].size(); i++)
				{
					mutralMatrix[tmpPair.first][i] = 1000000.0;
				}
				for (size_t i = 0; i<tmpPair.second; i++)
				{
					mutralMatrix[i][tmpPair.second-i] = 1000000.0;
				}
			}

		}while(smallest<1000000.0);

	}

	/************************************************************************/
	/* /brief Mobius transform the input point                                                                     */
	/************************************************************************/
	void GeometryProcessingAlgorithm::MobiusTransform(const complex<double>&sample, const vector<complex<double>>&m1, complex<double>& transformedSample)
	{
		assert( m1.size() == 4 );
		complex<double> numerator = m1[0]*sample + m1[1];
		complex<double> denominator =m1[2]*sample + m1[3];
		transformedSample = numerator/denominator;
	}

	/************************************************************************/
	/* /brief Find the closet point to sample in samplelist					*/
	/* /param [in]samplelist input sample list								*/
	/* /param [in]sample the sample point									*/
	/* /return the index of closest point in sample list, if failed to find	*/
	/*			return -1													*/
	/************************************************************************/
	int GeometryProcessingAlgorithm::MatchClosest(vector<complex<double>>& samplelist, complex<double>& sample)
	{
		assert(samplelist.size()>0);

		complex<double> tmpSample = samplelist[0];
		double realPart = abs(tmpSample.real() - sample.real());
		double imagePart = abs(tmpSample.imag() - sample.imag());
		double dis = sqrt( realPart*realPart + imagePart*imagePart );
		double tmpdis = dis;
		int matched = -1;
		for (size_t i = 0; i<samplelist.size(); i++)
		{
			complex<double> tmpSample = samplelist[i];
			double realPart = abs(tmpSample.real() - sample.real());
			double imagePart = abs(tmpSample.imag() - sample.imag());
			dis = sqrt( realPart*realPart + imagePart*imagePart );
			if ( dis <= tmpdis && dis < 1.0 )
			{
				tmpdis = dis;
				matched = i;
			}
		}
		return matched;
	}


	/************************************************************************/
	/* /brief Calculate the deformation energy with the transformation m1,m2*/
	/*																		*/
	/************************************************************************/
	double GeometryProcessingAlgorithm::CalculateDeformationEnergy(vector<complex<double>>&sourceMidEdgeSample,vector<complex<double>>&m1, vector<complex<double>>&targetMidEdgeSample,vector<complex<double>>&m2, vector<pair<int, int>>pairs)
	{
		double engery = 0;
		for (size_t i = 0; i<pairs.size(); i++)
		{
			pair<int, int> tmpPair = pairs[i];
			int indexSource = tmpPair.first;
			int indexTarget = tmpPair.second;
			double realPart = sourceMidEdgeSample[indexSource].real() - targetMidEdgeSample[indexTarget].real();
			double imagePart = sourceMidEdgeSample[indexSource].imag() - targetMidEdgeSample[indexTarget].imag();
			double dist = sqrt(realPart*realPart + imagePart*imagePart);
			engery = engery + dist;
		}
		return engery;
	}
}