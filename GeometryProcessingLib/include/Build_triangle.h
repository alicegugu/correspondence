#pragma once
#include <GL/glew.h>
#include <GL/glut.h>

#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <iostream>
#include <fstream>
#include "My_one2one_map.h"

#include "td.h"

using namespace std;
extern TDobject *_tdReadOBJObject(char *pathname);
namespace Another
{
	enum FileType
	{
		MFILE = 0,
		OBJFILE = 1
	};
	// A modifier creating a polyhedron with the incremental builder.
	template <class HDS>    
	class T_Build_triangle : public CGAL::Modifier_base<HDS> 
	{
	private:
		FILE* mf;
		FileType m_type;
		char* m_filename;
	public:
		T_Build_triangle() { mf = NULL; m_type = MFILE; }
		bool SetFile(FILE* m) 
		{  
			if ( m != NULL  )
			{
				mf = m; 
				return true;
			}
			else
			{
				return false;
			}
		}

		void SetFileType( FileType type)
		{
			m_type = type;
		}

		FileType GetFileType()
		{
			return m_type;
		}

		void SetFileName(char* name)
		{
			m_filename = name;
		}

		char* GetFileName()
		{
			return m_filename;
		}

		void operator()( HDS& hds) 
		{
			switch(m_type)
			{
			case MFILE:
				read_m_file(hds);
				break;
			case OBJFILE:
				read_obj_file(hds);
				break;
			default:
				break;
			}
		}

		
		void read_obj_file(HDS& hds)
		{
			/*
			//read obj file to TDobject object
			TDobject* object = tdReadObject(m_filename);

			CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);

			B.begin_surface( object->num_texvertices, object->num_faces);
			for (int i = 0; i<object->num_texvertices; i++)
			{
				TDvertex tmpv =  object->vertices[i];
				Vertex_handle vertexHandleO = B.add_vertex( Point( tmpv.x, tmpv.y, tmpv.z));
			}

			for ( size_t i = 0; i < object->num_groups; i++ )
			{
				TDgroup g = object->groups[i];
				for (size_t i=0; i<g.num_triangles; i++)
				{
					TDface f = g.triangles[i];
					B.begin_facet();
					B.add_vertex_to_facet( f.vertices[0]-1 );
					B.add_vertex_to_facet( f.vertices[1]-1 );
					B.add_vertex_to_facet( f.vertices[2]-1 );
					B.end_facet();
				}
			}
			*/
			cout<<"Reading obj files has linkage error"<<endl;

		}

		void read_m_file(HDS& hds)
		{
			CGAL_assertion_msg(mf,"file is not null");
			// Postcondition: `hds' is a valid polyhedral surface.
			int tmp_vertexNum, tmp_faceNum;
			fscanf(mf, "Vertex Num:%d Face Num:%d", &tmp_vertexNum, &tmp_faceNum);
			CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
			B.begin_surface( tmp_vertexNum, tmp_faceNum);
			typedef typename HDS::Vertex   Vertex;
			typedef typename Vertex::Point Point;
			typedef typename HDS::Vertex_handle Vertex_handle;
			float x, y, z, ox,oy,oz,normalx, normaly, normalz,onormalx, onormaly, onormalz;
			int vertexIndex, faceIndex, father;
			int point1, point2,point3;
			char type;
			int length;
			while( ( type = fgetc(mf) ) != EOF )
			{
				switch(type)
				{
				case 'V':
					// if(fscanf(mf, "ertex %d  %f %f %f {normal=(%f %f %f) Opos=(%f %f %f) Onormal=(%f %f %f) father=(%d)}\n", &vertexIndex, &x , &y, &z,  &normalx, &normaly, &normalz, &ox, &oy, &oz,  &onormalx, &onormaly, &onormalz, &father))
					//decocube
					//if(fscanf(mf, "ertex %d  %f %f %f {normal=(%f %f %f) Opos=(%f %f %f) Onormal=(%f %f %f) father=(%d)}\n", &vertexIndex, &x , &y, &z,  &normalx, &normaly, &normalz, &ox, &oy, &oz,  &onormalx, &onormaly, &onormalz, &father))

					//bunny: Vertex 1 0.333333 0.8125 0.0208333 {normal=(0.4572689048 0.7800744386 0.4270702741) Opos=(0.155348619 0.9310593523 0.09017550973) Onormal=(0.4572689048 0.7800744386 0.4270702741)}
					if( fscanf(mf, "ertex %d %f %f %f {normal=(%f %f %f) Opos=(%f %f %f) Onormal=(%f %f %f)}", 
						&vertexIndex, &x , &y, &z,  &normalx, &normaly, &normalz, &ox, &oy, &oz,  &onormalx, &onormaly, &onormalz) )
					{
						//Add Vertex to model and polyModel and set up their map between vertexes
						Vertex_handle vertexHandleO = B.add_vertex( Point( x, y, z));
						Point p(x,y,z);
						Normal n(normalx,normaly,normalz);
						vertexHandleO->SetNormal(n); //normal
						vertexHandleO->SetIndex(vertexIndex); //index 1-based
					}
					break;

				case 'F':
					if(fscanf(mf, "ace %d  %d %d %d", &faceIndex, &point1, &point2, &point3) )
					{
						B.begin_facet();
						B.add_vertex_to_facet(point1-1);
						B.add_vertex_to_facet(point2-1);
						B.add_vertex_to_facet(point3-1);
						B.end_facet();
					}
					break;

				case '#':
					do
					{
						type = fgetc(mf);
					}
					while (type != EOF && type != '\n' );
					break;
				}
			}

			B.end_surface();
		}
	};

}
