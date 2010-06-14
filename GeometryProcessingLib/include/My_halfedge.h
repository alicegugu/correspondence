#pragma once

#include "Traits.h"
#include "MeshElementBase.h"

namespace Another
{

    /***********************************************************************************************/
    /*********************A halfedge with the following additional information:***************/
    /**********************1. edge_type *****************************************************************/
    /**********************2. index ****************************************************************/
    /**********************3. edgeID **************************************************************/
    /************************************************************************************************/
    
    
    template<class Refs,  class Tprev, class Tvertex, class Tface, class Traits>
    struct My_halfedge:public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>, public Mesh_element_base
    {
    
        typedef CGAL::HalfedgeDS_halfedge_base<Refs,CGAL::Tag_true,CGAL::Tag_true,CGAL::Tag_true> BASECLASS;
		typedef typename Traits::Point_3 Point_3;
		typedef typename Refs::Vertex_handle         Vertex_handle;
		typedef typename Refs::Face_handle         Face_handle;
        private:       
            int edge_type;
            int index;
            int edgeID;
			double m_conjU;
			double m_U;
			double m_corner_angle; //the corner angle of the vertex that the halfedge point to
       
        public:
            bool detected;
			double m_len;
			bool visited;

         public: 
			 My_halfedge() : m_len(-1) { detected =  false; edge_type =UNDEFINED; index = UNDEFINED;edgeID = UNDEFINED; m_type = EDGE_TYPE;}
			double& get_length()  
			{
				return m_len; 
			}
            
            ~My_halfedge() {}
            
            My_halfedge( int t):BASECLASS() { edge_type = t; }
            
            bool SetIndex( int in ) { index = in; return true; }
            
            bool SetType( int t ) { edge_type = t; return true; }
            
            bool SetEdgeID( int id ) { edgeID = id; return true; }
            
            int GetIndex() { return index; }
            
            int GetEdgeID() { return edgeID; }
            
            int GetType() { return edge_type; }

			bool SetConjU(double conjU)
			{
				m_conjU = conjU;
				return true;
			}

			double GetConjU()
			{
				return m_conjU;
			}

			bool SetU(double u)
			{
				m_U = u;
				return true;
			}

			double GetU()
			{
				return m_U;
			}

			double get_corner_angle()
			{
				return m_corner_angle;
			}

			void set_corner_angle(double angle)
			{
				m_corner_angle = angle;
			}

			Face_handle opposite_face(Face_handle f)
			{
				if( is_border() )
				{
					assert( face() == f );
					return NULL;
				}

				assert(face() == f || opposite_face() == f );
				return  face() == f ? opposite_face() : face();
			};

			Face_handle opposite_face()
			{
				assert( !is_border() );
				Halfedge_handle e = opposite();
				return e->face();
			}

			Vertex_handle opposite_vertex(Vertex_handle v)
			{
				assert(belongs(v));
				return vertex() == v ? opposite_vertex(): vertex();
			};

			Vertex_handle opposite_vertex()
			{
				Halfedge_handle e = opposite();
				return e->vertex();
			};

			bool belongs(Vertex_handle v)
			{
				return vertex() == v || opposite_vertex() == v;
			}

			/************************************************************************/
			/*  /brief calculate the local coordinates for the input point          */
			/************************************************************************/
			void local_coordinates(Point_3 point, PlanePoint& localP)
			{
				Point_3 p1 = vertex()->point();
				Point_3 p2 = opposite_vertex()->point();
				double d0squared =  squared_length( p1 - point );
				double d0 = sqrt( d0squared );
				if(d0 < 1e-50)
				{
					localP[0] = 0.0;
					localP[1] = 0.0;
					return;
				}

				double d1squared = squared_length( p2 - point );
				double d1 = sqrt ( d1squared ) ;
				if(d1 < 1e-50)
				{
					localP[0] = m_len;
					localP[1] = 0.0;
					return;
				}

				localP[0] = m_len/2.0 + (d0squared - d1squared)/(2.0*m_len);
				localP[1] = sqrt(std::max(0.0, d0squared - localP[0]*localP[0]));
				return;
			}

    };
    
    
}