#pragma once

#include "Traits.h"
using namespace std;

namespace Another
{


	//----------------------------------------------------------------
	// some functors
	//----------------------------------------------------------------

	//computing the edge length
	struct Edge_length 
	{
		template <class HalfEdge_handle>
		double operator() (HalfEdge_handle h) {
			double l = CGAL::sqrt(CGAL::squared_distance(h->prev()->vertex()->point(),
				h->vertex()->point()));
			return l;
		}
	};

	//computing a the facet normal
	struct Facet_unit_normal {
		template < class Facet >
		typename Facet::Vector_3 operator() (Facet & f) {
			typename Facet::Halfedge_handle h = f.halfedge();
			typename Facet::Vector_3 normal =
				CGAL::cross_product(h->vertex()->point() -
				h->opposite()->vertex()->point(),
				h->next()->vertex()->point() -
				h->opposite()->vertex()->point());
			return normal / CGAL::sqrt(normal * normal);
		}
	};


	//----------------------------------------------------------------
	// operations on hedges, facets etc, handled using proeprty maps on the whole surface
	//----------------------------------------------------------------

	template <class TPoly , class HEdgePropertyMap> class T_Polyhedral_surface_hedge_operations
	{
	protected:
		typedef typename TPoly::Vertex Vertex;
		typedef typename TPoly::Halfedge_handle Halfedge_handle;
		typedef typename TPoly::Halfedge Halfedge;
		typedef typename TPoly::Facet_handle Facet_handle;
		typedef typename TPoly::Facet Facet;

		typedef typename TPoly::Vertex_iterator Vertex_iterator;
		typedef typename TPoly::Halfedge_iterator Halfedge_iterator;
		typedef typename TPoly::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;

	public:

		static void compute_edges_length(TPoly& P, HEdgePropertyMap& hepm);
		static double compute_mean_edges_length_around_vertex(Vertex * v, HEdgePropertyMap& hepm);
	};


	/************************************************************************/
	/* /brief Calculate the Ecurid length of each halfedge in P and save the*/
	/* length in the property map hepm,which contains (Halfedge_handle, len)*/
	/*	/param[in] TPoly& the input HDS										*/
	/*	/param[out] HEdgePropertyMap the property map used to store length	*/
	/************************************************************************/
	template <class TPoly , class HEdgePropertyMap>
	void T_Polyhedral_surface_hedge_operations<TPoly, HEdgePropertyMap>::
		compute_edges_length(TPoly& P, HEdgePropertyMap& hepm)
	{
		Halfedge_iterator itb = P.halfedges_begin(), ite = P.halfedges_end();
		for(; itb!=ite; itb++) {
			Halfedge_handle h = itb;
			double l = Edge_length()(h);
			put(hepm, h, l);
			//cerr << "[" << l << "," << get(hepm, h) << "] ";
		}
	}

	template <class TPoly , class HEdgePropertyMap>
	double T_Polyhedral_surface_hedge_operations<TPoly, HEdgePropertyMap>::
		compute_mean_edges_length_around_vertex(Vertex * v, HEdgePropertyMap& hepm)
	{
		Halfedge_around_vertex_circulator  hedgeb = v->vertex_begin(), hedgee = hedgeb;
		int  count_he = 0;
		double sum = 0.;
		CGAL_For_all(hedgeb, hedgee){
			sum += get(hepm, hedgeb);
			count_he++;
		}
		return sum/count_he;
	}


	//-------------------------------------------------------------------------------
	// operations on facets
	//-------------------------------------------------------------------------------
	template <class TPoly , class FacetPropertyMap> class T_Polyhedral_surface_face_operations
	{
	protected:
		typedef typename TPoly::Vertex Vertex;
		typedef typename TPoly::Halfedge_handle Halfedge_handle;
		typedef typename TPoly::Halfedge Halfedge;
		typedef typename TPoly::Facet_handle Facet_handle;
		typedef typename TPoly::Facet Facet;

		typedef typename TPoly::Vertex_iterator Vertex_iterator;
		typedef typename TPoly::Halfedge_iterator Halfedge_iterator;
		typedef typename TPoly::Facet_iterator Facet_iterator;
		typedef typename TPoly::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;

	public:
		static void compute_facets_normals(TPoly& P, FacetPropertyMap& fpm);
		static typename TPoly::Traits::Vector_3 compute_vertex_average_unit_normal(Vertex * v, const FacetPropertyMap& fpm);
	};




	template <class TPoly , class FacetPropertyMap>
	void T_Polyhedral_surface_face_operations<TPoly, FacetPropertyMap>::
		compute_facets_normals(TPoly& P, FacetPropertyMap& fpm)
	{
		//iterator returning facet handles
		Facet_iterator itb = P.facets_begin(), ite = P.facets_end();
		for(; itb!=ite; itb++)
			put(fpm, itb, Facet_unit_normal()(*itb));
	}


	template <class TPoly , class FacetPropertyMap>
	typename TPoly::Traits::Vector_3 T_Polyhedral_surface_face_operations<TPoly, FacetPropertyMap>::
		compute_vertex_average_unit_normal(Vertex * v, const FacetPropertyMap& fpm)
	{
		Facet_handle f;
		Vector_3 sum(0., 0., 0.), n;

		Halfedge_around_vertex_circulator hedgeb = v->vertex_begin(), hedgee = hedgeb;
		do{
			if (hedgeb->is_border_edge()){
				hedgeb++;
				continue;
			}

			f = hedgeb->facet();
			n = get(fpm, f);
			sum = (sum + n);
			hedgeb++;
		}
		while (hedgeb != hedgee);
		sum = sum / std::sqrt(sum * sum);
		return sum;
	}

}
