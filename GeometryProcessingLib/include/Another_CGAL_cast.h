#pragma once

#include "MeshElementBase.h"
#include "Another_HDS_model.h"

namespace Another
{

	/************************************************************************/
	/* Used for cast base class to inherited class's handle                 */
	/************************************************************************/
	class Another_CGAL_cast
	{
	public:
		Another_CGAL_cast(void) {}
		~Another_CGAL_cast(void) {}

		/************************************************************************/
		/*  /brief cast a base point to Vertex_handle                           */
		/************************************************************************/
		static Vertex_handle vertex_cast(Base_handle p)
		{
			assert( p->type() == VERTEX_TYPE );
			Vertex* vp = static_cast<Vertex*>(p);
			Halfedge_handle e = vp->halfedge();
			return e->vertex();
		}

		/************************************************************************/
		/*  /brief cast a Vertex_handle to Base class point                     */
		/************************************************************************/
		static Base_handle vertex_cast(Vertex_handle p)
		{
			assert( p->type() == VERTEX_TYPE );
			return static_cast<Base_handle>(&(*p));
		}

		/************************************************************************/
		/*  /brief cast a vertex point to Vertex_handle                         */
		/************************************************************************/
		static Vertex_handle vertex_cast(Vertex* vp)
		{
			assert( vp->type() == VERTEX_TYPE );
			Halfedge_handle e = vp->halfedge();
			return e->vertex();
		}

		/************************************************************************/
		/*  /brief cast a base point to Face_handle                             */
		/************************************************************************/
		static Face_handle face_cast(Base_handle p)
		{
			assert( p->type() == FACE_TYPE );
			Face* fp = static_cast<Face*>(p);
			Halfedge_handle e = fp->halfedge();
			return e->face();

		}

		/************************************************************************/
		/*  /brief cast a Face_handle to Base_handle                            */
		/************************************************************************/
		static Base_handle face_cast(Face_handle fp)
		{
			assert( fp->type() == FACE_TYPE );
			return static_cast<Base_handle>(&(*fp));

		}

		/************************************************************************/
		/*  /brief cast a Face point to Face_handle                             */
		/************************************************************************/
		static Face_handle face_cast(Face* fp)
		{
			assert( fp->type() == FACE_TYPE );
			Halfedge_handle e = fp->halfedge();
			return e->face();

		}

		/************************************************************************/
		/*  /brief cast a base point to Halfedge_handle                         */
		/************************************************************************/
		static Halfedge_handle halfedge_cast(Base_handle p)
		{
			assert( p->type() == EDGE_TYPE );
			Edge* ep = static_cast<Edge*>(p);
			Face_handle f = ep->face();
			Halfedge_handle e = f->halfedge();
			for (unsigned i = 0; i<3; i++)
			{
				if ( &(*e) == ep )
				{
					return e;
				}
				else
				{
					e = e->next();
				}
			}
			assert(0);
		}

		/************************************************************************/
		/*  /brief cast a Halfedge point to Halfedge_handle                     */
		/************************************************************************/
		static Halfedge_handle halfedge_cast(Edge* ep)
		{
			assert( ep->type() == EDGE_TYPE );
			Face_handle f = ep->face();
			Halfedge_handle e = f->halfedge();
			for (unsigned i = 0; i<3; i++)
			{
				if ( &(*e) == ep )
				{
					return e;
				}
				else
				{
					e = e->next();
				}
			}
			assert(0);
		}

		/************************************************************************/
		/*  /brief cast a Face_handle to Base_handle                            */
		/************************************************************************/
		static Base_handle halfedge_cast(Halfedge_handle ep)
		{
			assert( ep->type() == EDGE_TYPE );
			return static_cast<Base_handle>(&(*ep));
		}
	};
}
