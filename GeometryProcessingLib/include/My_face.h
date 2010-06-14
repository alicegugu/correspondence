#pragma once

#include "MeshElementBase.h"
#include <vector>

namespace Another
{
	class Mesh_element_base;
    template<class Refs, class GeoTraits>
    struct My_face:public CGAL::HalfedgeDS_face_base< Refs, CGAL::Tag_true, Plane>, public Mesh_element_base
    {
	public:
		typedef typename GeoTraits::Vector_3 Vector_3;
        typedef CGAL::HalfedgeDS_face_base< Refs, CGAL::Tag_true, Plane> BASECLASS;
		typedef typename Refs::Vertex_handle Vertex_handle;
		typedef typename Refs::Halfedge_handle Halfedge_handle;

		Vector_3 normal;

        int face_type;
        int index;
        int faceID;
        bool visited;
        Normal n;
        bool select;

	private:
		double m_corner_angles[3];		//triangle angles in radians; angles correspond to vertices in m_adjacent_vertices

        public:
            int ZeroOrPI; //ZEROANGLE = -1 PIANGLE = 1 NORMALANGLE = 0
            ~My_face() {}
            //  My_face(const My_face& f):BASECLASS() {}
            My_face():BASECLASS() { face_type = UNDEFINED;index = UNDEFINED; faceID = UNDEFINED;visited = false; ZeroOrPI =NORMALANGLE; m_type = FACE_TYPE; }
            //  My_face():BASECLASS() { }
            bool SetType( int t ) { face_type = t; return true; }
            bool SetIndex( int in ) { index = in; return true; }
            bool SetFaceID( int id ) { faceID = id; return true; }
            bool SetNormal (Normal& in) { n = in; return true;}
            int GetIndex() { return index; }
            int GetFaceID() { return faceID; }
            int GetType() { return face_type; }
            Normal GetNormal() { return n; }
            bool SetSelect( bool se) { select = se; return select; }
            bool GetSelect() {return select;}
			const Vector_3& get_unit_normal() const { return normal; }
			Vector_3& get_unit_normal() { return normal; }
			double* corner_angles(){return m_corner_angles;};
			Halfedge_handle opposite_edge(Vertex_handle v)
			{
				Halfedge_handle e = halfedge();
				for(unsigned i=0; i<3; ++i)
				{
					if(!e->belongs(v))
					{
						return e;
					}
					else
					{
						e = e->next();
					}
				}
				assert(0);
				return NULL;
			}
			Vertex_handle opposite_vertex(Halfedge_handle e)
			{
				Halfedge_handle tmpe = halfedge();

				for(unsigned i=0; i<3; ++i)
				{
					Vertex_handle v = tmpe->vertex();
					if(!e->belongs(v))
					{
						return v;
					}
					else
					{
						tmpe = tmpe->next();
					}
				}
				assert(0);
				return NULL;
			}

			/************************************************************************/
			/*  /brief Return the edge shared the vertex with the input edge        */
			/*  /param [input] input edge											*/
			/*	/param [input] shared vertex										*/
			/*	/param [output] output edge											*/
			/************************************************************************/
			Halfedge_handle next_edge(Halfedge_handle e, Vertex_handle v)
			{
				assert(e->belongs(v));
				Halfedge_handle next = halfedge();
				for(unsigned i=0; i<3; ++i)
				{
					if(e != next && next->belongs(v))
					{
						return next;
					}
					else
					{
						next = next->next();
					}
				}
				assert(0);
				return NULL;
			}

			/************************************************************************/
			/* /brief Calculate the angle of the corner with input vertex           */
			/************************************************************************/
			double vertex_angle(Vertex_handle v)
			{
				Halfedge_handle e = halfedge();
				for(unsigned i=0; i<3; ++i)
				{
					Vertex_handle tmpv = e->vertex();
					if( tmpv == v )
					{
						return e->get_corner_angle();
					}
					else
					{
						e = e->next();
					}
				}
				assert(0);
				return 0;
			}
			void get_all_vertex(std::vector<Vertex_handle>& vertexs )
			{
				assert(vertexs.size());
				Halfedge_handle e = halfedge();
				for (int i = 0; i< 3; i++)
				{
					Vertex_handle v = e->vertex();
					vertexs.push_back(v);
					e = e->next();
				}
				return;
			}
        };
}