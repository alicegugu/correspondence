#pragma once

#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Lapack/Linear_algebra_lapack.h>

#include <boost/property_map/property_map.hpp>

#include <stdexcept>
#include <map>
#include <complex>

#include "Traits.h"
#include "Another_HDS_model.h"
#include "Polyhedral_surface_rings.h"
#include "Polyhedral_surface_operations.h"


namespace Another
{
	// \struct used to comapre the address of Half edge
	struct Hedge_cmp 
	{
		bool operator()(Halfedge_handle a, Halfedge_handle b) const
		{
			return &*a <&*b;
		}
	};

	struct Face_cmp
	{
		bool operator()(Face_handle a, Face_handle b) const
		{
			return &*a < &*b;
		}
	};

	struct Vertex_cmp
	{
		bool operator()(Vertex_handle a, Vertex_handle b) const
		{
			return &*a < &*b;
		}
	};


	typedef std::map<Vertex_handle, double> Vertex2double_map_type;
	typedef boost::associative_property_map<Vertex2double_map_type> Vertex_PM_double;

	typedef std::map<Vertex_handle, int, Vertex_cmp> Vertex2int_map_type;
	/// \typedef Vertex property map, with std::map
	typedef boost::associative_property_map<Vertex2int_map_type> Vertex_PM_type;
	typedef T_Polyhedral_surface_rings<Another_HDS_model, Vertex_PM_type> Poly_rings;


	typedef std::map<Halfedge_handle, double, Hedge_cmp> Hedge2double_map_type;
	/// \typedef Hedge property map, with enriched halfedge with its length
	typedef boost::associative_property_map<Hedge2double_map_type> Hedge_PM_type;
	typedef T_Polyhedral_surface_hedge_operations<Another_HDS_model, Hedge_PM_type> Hedge_ops;


	typedef std::map<Face_handle, Normal, Face_cmp> Face2normal_map_type;
	/// \typedef Face property map with enriched face with its normal
	typedef boost::associative_property_map<Face2normal_map_type> Face_PM_type;
	typedef T_Polyhedral_surface_face_operations<Another_HDS_model, Face_PM_type> Face_ops;

	typedef CGAL::Monge_via_jet_fitting<Kernel> Another_Monge_via_jet_fitting;
	typedef Another_Monge_via_jet_fitting::Monge_form Another_Monge_form;

	typedef std::map<Vertex_handle,Another_Monge_form, Vertex_cmp> Vertex2Monge_map_type;
	typedef boost::associative_property_map<Vertex2Monge_map_type> Vertex2Monge_PM_type;

	typedef std::map<Vertex_handle,double, Vertex_cmp> Vertex2GuassianCurvature_map_type;
	typedef boost::associative_property_map<Vertex2GuassianCurvature_map_type> Vertex2GuassianCurvature_PM_type;
	/************************************************************************/
	/* \brief Geometry Processing Algorithms for Paint Model                */
	/************************************************************************/
	class GeometryProcessingAlgorithm
	{
	public:
		GeometryProcessingAlgorithm(void);
		~GeometryProcessingAlgorithm(void);

		//Gather points around the vertex using one-neighhor-ring. copied from CGAL/examples/...
		void GatheringNeighborRings(Vertex* v, Vertex_PM_type& vpm ,std::vector<Point_3>&neighhor, int neighborPointsToUse, int min_nb_points, int nb_rings);

		//Caculate Gaussian curvature for each point using Monge jet fitting
		Vertex2GuassianCurvature_PM_type& GaussianCurvature(Vertex2Monge_PM_type& monges, Another_HDS_model& model, Vertex2GuassianCurvature_PM_type& GCpm);
		//Caculate Monge (k1, k2, d1, d2,..etc)
		Vertex2Monge_PM_type& Monge_jet_fitting(Another_HDS_model& model, int min_nb_points, int d_fitting, int d_monge, int nb_rings, int neighborPointsToUse, Vertex2Monge_PM_type&monges);
		void LocalMaxandMinGuassian(Another_HDS_model& model, vector<Vertex_handle>& featurePoints);
		void LocalMaxandMinGuassian(Another_HDS_model& model, float radius, vector<Vertex_handle>& featurePoints);
		double* MobiusVoting(Another_HDS_model& source, Another_HDS_model& target);
		void Embed(Another_HDS_model& source, Another_HDS_model& target);
	private:
		void GetMobiusTransformation(vector<complex<double>>& z, vector<complex<double>>& transformation);
		void ProjectSamples(vector<Vertex_handle>& Points,vector<complex<double>>& MidEdgeSample);
		void MobiusTransform(const complex<double>&sample, const vector<complex<double>>&m1, complex<double>& transformedSample);
		int MatchClosest(vector<complex<double>>& samplelist, complex<double>& sample);
		void FindMutualNeightbor(vector<complex<double>>& sourceMidEdgeSample,const vector<complex<double>>& m1, vector<complex<double>>& targetMidEdgeSample, const vector<complex<double>>&m2,vector<pair<int, int>>& pairs);
		double CalculateDeformationEnergy(vector<complex<double>>&sourceMidEdgeSample,vector<complex<double>>&m1, vector<complex<double>>&targetMidEdgeSample,vector<complex<double>>&m2, vector<pair<int, int>>pairs);
		bool ClearPropertyMap(Another_HDS_model& model, Vertex_PM_type& vpm, float defaultValue);
		//HarmonicMap
	};

	//template function, todo
	//vector<Vertex_handle> LocalMaxandMin(Another_HDS_model& model, Vertex_PM_double& vpm);
}