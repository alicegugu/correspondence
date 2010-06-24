#include <math.h>

#include ".\include\Another_HDS_model.h"
#include ".\include\My_one2one_map.h"
#include "geodesic_mesh_elements.h"
#include "Polyhedral_surface_operations.h"

using namespace std;

extern double compute_angle_rad(Another::Point P, Another::Point Q, Another::Point R);

namespace Another
{
	const double color_epsilon = 10e-6;
	GLint    viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	GLfloat g_color_map[300][3];
	const int FEATURE_CORRESPONDENCE = 1;

	//My_One2OneMap<unsigned int,POLYCUBEMODEL_HDS::Vertex_handle,SORT_POLYCUBEMODEL_VERTEX> my_poly_corner_map;
	//	My_One2OneMap<unsigned int,Another_HDS::Vertex_iterator,SORT_PAINTMODEL_VERTEX> my_corner_map;
	//My_One2OneMap<unsigned int,Poly_Edge_handle,SORT_POLYCUBEMODEL_EDGE> my_poly_edge_map;
	//My_One2OneMap<unsigned int, Poly_Face_handle, SORT_POLYCUBEMODEL_FACE> my_poly_face_map;


	//     GLfloat deltaDistanceO[3] = { 0.0, 0.0, 0.0 };

	//     GLfloat yRotO = 0;
	//     GLfloat scaleO = 1;
	//     GLfloat deltaDistanceP[3] = { 0.0, 0.0, 0.0 };
	//     GLfloat scaleP = 1;
	//     //tablet position
	//     extern int mouse_x;
	//     extern int mouse_y;
	//     extern int windowWdith;
	//     extern int windowHeight;
	//     extern Point my_cursorO;
	//     extern Point my_cursorP;
	//     extern  float fovy;
	//     extern  float xy_aspect;
	//     extern  float zNear;
	//     extern  float zFar;

	//     GLfloat poly_upX;
	//     GLfloat poly_upY;
	//     GLfloat poly_upZ;
	//     GLfloat poly_eyeX;
	//     GLfloat poly_eyeY;
	//     GLfloat poly_eyeZ;
	//     GLfloat poly_centerX;
	//     GLfloat poly_centerY;
	//     GLfloat poly_centerZ;


	Another_HDS_model::Another_HDS_model( void )
	{

		//set default draw type as smooth
		m_draw_type = SMOOTH;

		//set default color as white
		m_color[0] = 0.8;
		m_color[1] = 0.8;
		m_color[2] = 0.8;
		m_color[3] = 1.0;

		//set default material
		m_ambient_matrix[0] = 0.7;
		m_ambient_matrix[1] = 0.7;
		m_ambient_matrix[2] = 0.7;
		m_ambient_matrix[3] = 1.0;

		m_diffuse_matrix[0] = 0.8;
		m_diffuse_matrix[1] = 0.8;
		m_diffuse_matrix[2] = 0.8;
		m_diffuse_matrix[3] = 1.0;

		m_shininess = 5.0;

		m_specular_matrix[0] = 1.0;
		m_specular_matrix[1] = 1.0;
		m_specular_matrix[2] = 1.0;
		m_specular_matrix[3] = 1.0;


		m_mode = DEBUGMODE;
		m_draw_axis = false;
		m_has_UV = false;
		m_has_gaussian_curvature = false;
		m_has_1D_u = false;
	}

	Another_HDS_model::~Another_HDS_model( void )
	{
	}

	int Another_HDS_model::closest_vertices(Surface_point& p, std::vector<Vertex_handle>& storage)
	{
		assert(p.type() != UNDEFINED_POINT_TYPE);
		storage.clear();

		if(p.type() == VERTEX_TYPE)
		{

			storage.push_back(Another_CGAL_cast::vertex_cast(p.base_element()));

			return 1;
		}
		else if(p.type() == FACE_TYPE)
		{

			Face_handle fp = Another_CGAL_cast::face_cast(p.base_element());
			vector<Vertex_handle> vs;
			fp->get_all_vertex(vs);
			for (unsigned int i = 0; i<3; i++)
			{
				storage.push_back(vs[i]);
			}

			return 2;
		}
		else if(p.type() == EDGE_TYPE)		//for edge include all 4 adjacent vertices
		{
			Halfedge_handle edge = Another_CGAL_cast::halfedge_cast(p.base_element());


			storage.push_back(edge->vertex());
			storage.push_back(edge->opposite_vertex());
			if (edge->is_border())
			{
				Face_handle face = edge->face();
				storage.push_back( &(*face->opposite_vertex(edge)) );
				return 3;
			}
			else
			{
				Face_handle face = edge->face();
				storage.push_back(&(*face->opposite_vertex(edge)));
				face = edge->opposite_face();
				edge = edge->opposite();
				storage.push_back(&(*face->opposite_vertex(edge)));
				return 4;
			}
		}

		assert(0);
		return 0;
	}

	/************************************************************************/
	/* /brief Draw mid edge features                                        */
	/* /param mode[in] GL_SELCET selection mode GL_SMOOTH render mode		*/
	/************************************************************************/
	void Another_HDS_model::draw_mid_edge_features(int mode)
	{
		if (m_middle_edge_features.size() == 0)
		{
			cout<<"Error:Another_HDS_model::draw_mid_edge_features: no feature point to render."<<endl;
			return;
		}

		for (int i=0; i<m_middle_edge_features.size(); i++)
		{
			complex<double> feature = m_middle_edge_features[i];
			if (mode == GL_SELECT)
			{
				glShadeModel(GL_SELECT);
			}
			else
			{
				glShadeModel(GL_SMOOTH);
			}

			glPushMatrix();
			glTranslatef(feature.real()*5.0f ,feature.imag()*5.0f,1.0f);
			if (mode == GL_SELECT)
			{
				glPushName(i);
			}
			if(m_middle_edge_features_select[i])
			{	
				glColor3f(1.0f,0.0f,0.0f);
			}
			else
			{
				glColor3f(0.0f,1.0f,0.0f);
			}
			glutSolidSphere(0.2,16,16);
			if (mode == GL_SELECT)
			{
				glPopName();
			}
			glPopMatrix();
		}
	}

	bool Another_HDS_model::draw_midedge( )
	{
		glShadeModel(GL_SMOOTH);
		glBegin(GL_TRIANGLES);

		Face_iterator faceIt = facets_begin();
		Face_iterator faceEnd = facets_end();
		for( ; faceIt!= faceEnd; faceIt++ )
		{
			//skip the face has been cut
			if ( faceIt->GetFaceID() == m_cutface )
			{
				continue;
			}

			Halfedge_handle tmpedge;
			GLfloat vertexArray1[3];
			GLfloat vertexArray2[3];
			GLfloat vertexArray3[3];
			tmpedge = faceIt->halfedge();
			vertexArray1[0]= tmpedge->GetU()*5.0;
			vertexArray1[1]= tmpedge->GetConjU()*5.0;
			vertexArray1[2] = 1.0;

			tmpedge = tmpedge->next();
			vertexArray2[0]= tmpedge->GetU()*5.0;
			vertexArray2[1]= tmpedge->GetConjU()*5.0;
			vertexArray2[2] = 1.0;

			tmpedge = tmpedge->next();
			vertexArray3[0]= tmpedge->GetU()*5.0;
			vertexArray3[1]= tmpedge->GetConjU()*5.0;
			vertexArray3[2] = 1.0;
			//glColor3i(255,0,0);
			glVertex3fv(vertexArray1);
			glVertex3fv(vertexArray2);
			glVertex3fv(vertexArray3);
		}

		glEnd();

		return true;

	}


	void Another_HDS_model::draw_HamonicMap()
	{
		glShadeModel(GL_SMOOTH);

		enable_material();

		glBegin(GL_TRIANGLES);

#pragma region Draw harmonic texture mapping

		GLfloat vertexArray[3];
		GLfloat normalArray[3];
		float tmp_u;
		Vertex_handle v;
		Point p;
		Normal n;
		Face_iterator faceIt = facets_begin();
		Halfedge_handle tmpedge,firstedge;

		//Get the  location in memory of the texture coordinates
		//	GLint texCoordinates = glGetAttribLocation( program,char *name);

		for( ; faceIt!= facets_end(); faceIt++ )
		{
			//Skip the face to cut
			if ( faceIt->GetFaceID() == m_cutface )
			{
				continue;
			}
			//Set normal and position for vertex 2
			tmpedge = faceIt->halfedge();
			int faceID = faceIt->GetFaceID();
			v = tmpedge->vertex();
			p = v->point();
			n = v->GetNormal();
			vertexArray[0] = p.x();
			vertexArray[1] = p.y();
			vertexArray[2] = p.z();
			normalArray[0] = n.x();
			normalArray[1] = n.y();
			normalArray[2] = n.z();
			tmp_u = v->GetU();

			glNormal3fv(normalArray);
			glTexCoord2f(tmp_u, 0.0);
			glVertex3fv(vertexArray);
			//		glVertexAttrib2f(texCoordinates, GLfloat v0, GLfloat v1);

			//Set normal and position for vertex 2
			tmpedge=tmpedge->next();
			v = tmpedge->vertex();
			p = v->point();
			n = v->GetNormal();
			vertexArray[0] = p.x();
			vertexArray[1] = p.y();
			vertexArray[2] = p.z();
			normalArray[0] = n.x();
			normalArray[1] = n.y();
			normalArray[2] = n.z();
			tmp_u = v->GetU();
			glNormal3fv(normalArray);
			glTexCoord2f(tmp_u, 0.0);
			glVertex3fv(vertexArray);

			tmpedge=tmpedge->next();
			v = tmpedge->vertex();

			//Set normal and position for vertex 3
			p = v->point();
			n = v->GetNormal();
			vertexArray[0] = p.x();
			vertexArray[1] = p.y();
			vertexArray[2] = p.z();
			normalArray[0] = n.x();
			normalArray[1] = n.y();
			normalArray[2] = n.z();
			tmp_u = v->GetU();
			glNormal3fv(normalArray);
			glTexCoord2f(tmp_u, 0.0);
			glVertex3fv(vertexArray);
		}
#pragma endregion 

		glEnd();
	}

	bool Another_HDS_model::draw(RenderOption renderOption)
	{
		switch(renderOption)
		{
		case RenderOption::FEATURE_POINT:
			draw_feature_points(GL_SMOOTH);
			set_curvature();
			draw(SMOOTH,1,1);
			break;
		case RenderOption::CURVATURE:
			set_curvature();
			draw(SMOOTH,1,1);
			break;
		case RenderOption::DISCRETE_HARMONIC_MAP:
			draw_HamonicMap();
			break;
		case RenderOption::MESH:
			un_set_curvature();
			draw(SMOOTH,1,1);
			break;
		case RenderOption::PLANAR_EMBEDDING:
			draw_mid_edge_features(GL_SELECT);
			draw_midedge();
			break;
		case RenderOption::TRIANGULAR:
			draw(WIREFRAME,1,1);
		case RenderOption::LAPLACE:
			draw_laplacian();
			break;
		case RenderOption::ALIGNMENT_SELECT:
			draw(SMOOTH,1,1);
			draw_feature_points(ALIGNMENT_SELECT);
			break;
		case RenderOption::DEFAULT:
			draw(SMOOTH,1,1);
		}
		return true;
	}

	void Another_HDS_model::draw_laplacian()
	{
		glShadeModel(GL_SMOOTH);
		enable_material();
		glBegin(GL_TRIANGLES);
		GLfloat vertexArray[3];
		GLfloat normalArray[3];
		float tmp_u;
		float tmp_v;
		Vertex_handle v;
		Point p;
		Normal n;
		Face_iterator faceIt = facets_begin();
		Halfedge_handle tmpedge,firstedge;

		for( ; faceIt!= facets_end(); faceIt++ )
		{
			tmpedge = faceIt->halfedge();
			int faceID = faceIt->GetFaceID();

			v = tmpedge->vertex();

			glColor4fv( m_color );
			p = v->m_laplace_position;
			n = v->GetNormal();
			vertexArray[0] = p.x();
			vertexArray[1] = p.y();
			vertexArray[2] = p.z();
// 			normalArray[0] = n.x();
// 			normalArray[1] = n.y();
// 			normalArray[2] = n.z();




			glNormal3fv(normalArray);

			if (m_has_gaussian_curvature)
			{
				double curvature = v->Get1DU();
				glTexCoord2f(curvature, curvature);
			}
			if (m_has_UV)
			{
				tmp_u = v->GetU();
				tmp_v = v->GetV();
				glTexCoord2f(tmp_v,tmp_u);
			}
			glVertex3fv(vertexArray);

			tmpedge=tmpedge->next();
			v = tmpedge->vertex();
			// 						if(v->GetType() == CORNER || v->GetType() == EDGE)
			// 						{
			// 							glColor4fv(m_corner_color);
			// 						}
			// 						else
			// 						{
			// 							glColor4fv( m_color );
			// 						}
			glColor4fv( m_color );
			p = v->m_laplace_position;
			n = v->GetNormal();
			vertexArray[0] = p.x();
			vertexArray[1] = p.y();
			vertexArray[2] = p.z();
			normalArray[0] = n.x();
			normalArray[1] = n.y();
			normalArray[2] = n.z();
			glNormal3fv(normalArray);
			if (m_has_gaussian_curvature)
			{
				double curvature = v->Get1DU();
				glTexCoord2f(curvature,curvature);
			}
			if (m_has_UV)
			{
				tmp_u = v->GetU();
				tmp_v = v->GetV();
				glTexCoord2f(tmp_v,tmp_u);
			}
			glVertex3fv(vertexArray);

			tmpedge=tmpedge->next();
			v = tmpedge->vertex();
			// 						if(v->GetType() == CORNER || v->GetType() == EDGE)
			// 						{
			// 							glColor4fv(m_corner_color);
			// 						}
			// 						else
			// 						{
			// 							glColor4fv( m_color );
			// 						}
			glColor4fv( m_color );
			p = v->m_laplace_position;
			n = v->GetNormal();
			vertexArray[0] = p.x();
			vertexArray[1] = p.y();
			vertexArray[2] = p.z();
			normalArray[0] = n.x();
			normalArray[1] = n.y();
			normalArray[2] = n.z();
			glNormal3fv(normalArray);
			if (m_has_gaussian_curvature)
			{
				double curvature = v->Get1DU();
				glTexCoord2f(curvature, curvature);
			}
			if (m_has_UV)
			{
				tmp_u = v->GetU();
				tmp_v = v->GetV();
				glTexCoord2f(tmp_v,tmp_u);
			}
			glVertex3fv(vertexArray);
		}
		glEnd();

	}
	void Another_HDS_model::get_guassian_curvature(vector<double>& gc)
	{
		for (Another::Vertex_iterator it = vertices_begin(); it != vertices_end(); ++it)
		{
			gc.push_back(it->GetGaussianCurvature());
		}
	}


	bool Another_HDS_model::set_1D_texture_coordinates(vector<double>& texture_coordinates)
	{
		if (texture_coordinates.size()!= size_of_vertices())
		{
			return false;
		}
		int i = 0;
		for (Another::Vertex_iterator it = vertices_begin(); it != vertices_end(); ++it)
		{
			it->Set1DU(texture_coordinates[i]);
			i++;
		}
		m_has_1D_u = true;
		return true;
	}


	bool Another_HDS_model::draw(int type,int PolyOrOriginal,int debug)
	{
		//draw_XYZ_Normals();

		switch(type)
		{
		case WIREFRAME:
			{
				if(m_mode == DEBUGMODE)
				{
					glBegin(GL_LINES);
					Face_iterator faceIt = facets_begin();
					Halfedge_handle tmpedge,firstedge;
					Vertex_handle v1,v2,v3;
					Point p;
					GLfloat vertexArray1[3];
					GLfloat vertexArray2[3];
					GLfloat vertexArray3[3];
					for( ; faceIt!= facets_end(); faceIt++ )
					{
						tmpedge = faceIt->halfedge();
						v1 = tmpedge->vertex();
						p = v1->point();
						vertexArray1[0] = p.x();
						vertexArray1[1] = p.y();
						vertexArray1[2] = p.z();

						tmpedge=tmpedge->next();
						v2 = tmpedge->vertex();
						p = v2->point();
						vertexArray2[0] = p.x();
						vertexArray2[1] = p.y();
						vertexArray2[2] = p.z();
						glColor4fv( m_color );
						glVertex3fv(vertexArray1);
						glVertex3fv(vertexArray2);



						tmpedge=tmpedge->next();
						v3 = tmpedge->vertex();
						p = v3->point();
						vertexArray3[0] = p.x();
						vertexArray3[1] = p.y();
						vertexArray3[2] = p.z();
						glColor4fv( m_color );

						glVertex3fv(vertexArray2);
						glVertex3fv(vertexArray3);

						//               if((v1->GetType() == CORNER || v1->GetType() == EDGE) &&(v3->GetType() == CORNER || v3->GetType() == EDGE))
						//                   {
						//                   glColor4fv(m_corner_Color);
						//                   }
						//               else
						//                   {
						//                   glColor4fv( m_color );
						//                   }
						glVertex3fv(vertexArray1);
						glVertex3fv(vertexArray3);
					}
					glEnd();
				}

				break;
			}
		case SMOOTH:
			{
				if(m_mode ==DEBUGMODE)
				{
					glShadeModel(GL_SMOOTH);
					enable_material();
					glBegin(GL_TRIANGLES);
					GLfloat vertexArray[3];
					GLfloat normalArray[3];
					float tmp_u;
					float tmp_v;
					Vertex_handle v;
					Point p;
					Normal n;
					Face_iterator faceIt = facets_begin();
					Halfedge_handle tmpedge,firstedge;

					for( ; faceIt!= facets_end(); faceIt++ )
					{
						tmpedge = faceIt->halfedge();
						int faceID = faceIt->GetFaceID();

						v = tmpedge->vertex();

						glColor4fv( m_color );
						p = v->point();
						n = v->GetNormal();
						vertexArray[0] = p.x();
						vertexArray[1] = p.y();
						vertexArray[2] = p.z();
						normalArray[0] = n.x();
						normalArray[1] = n.y();
						normalArray[2] = n.z();




						glNormal3fv(normalArray);

						if (m_has_gaussian_curvature)
						{
							double curvature = v->Get1DU();
							glTexCoord2f(curvature, curvature);
						}
						if (m_has_UV)
						{
							tmp_u = v->GetU();
							tmp_v = v->GetV();
							glTexCoord2f(tmp_v,tmp_u);
						}
						glVertex3fv(vertexArray);

						tmpedge=tmpedge->next();
						v = tmpedge->vertex();
						// 						if(v->GetType() == CORNER || v->GetType() == EDGE)
						// 						{
						// 							glColor4fv(m_corner_color);
						// 						}
						// 						else
						// 						{
						// 							glColor4fv( m_color );
						// 						}
						glColor4fv( m_color );
						p = v->point();
						n = v->GetNormal();
						vertexArray[0] = p.x();
						vertexArray[1] = p.y();
						vertexArray[2] = p.z();
						normalArray[0] = n.x();
						normalArray[1] = n.y();
						normalArray[2] = n.z();
						glNormal3fv(normalArray);
						if (m_has_gaussian_curvature)
						{
							double curvature = v->Get1DU();
							glTexCoord2f(curvature,curvature);
						}
						if (m_has_UV)
						{
							tmp_u = v->GetU();
							tmp_v = v->GetV();
							glTexCoord2f(tmp_v,tmp_u);
						}
						glVertex3fv(vertexArray);

						tmpedge=tmpedge->next();
						v = tmpedge->vertex();
						// 						if(v->GetType() == CORNER || v->GetType() == EDGE)
						// 						{
						// 							glColor4fv(m_corner_color);
						// 						}
						// 						else
						// 						{
						// 							glColor4fv( m_color );
						// 						}
						glColor4fv( m_color );
						p = v->point();
						n = v->GetNormal();
						vertexArray[0] = p.x();
						vertexArray[1] = p.y();
						vertexArray[2] = p.z();
						normalArray[0] = n.x();
						normalArray[1] = n.y();
						normalArray[2] = n.z();
						glNormal3fv(normalArray);
						if (m_has_gaussian_curvature)
						{
							double curvature = v->Get1DU();
							glTexCoord2f(curvature, curvature);
						}
						if (m_has_UV)
						{
							tmp_u = v->GetU();
							tmp_v = v->GetV();
							glTexCoord2f(tmp_v,tmp_u);
						}
						glVertex3fv(vertexArray);
					}
					glEnd();
					// 					if (my_corner_map.begin()!=my_corner_map.end())
					// 					{
					// 						DrawCorners();
					// 					}
				}

				//Selection mode, which the user can select facet to draw
				else if ( m_mode == SELECTION_MODE )
				{

					//Draw the face
					//	this->draw_all_facets(1);

				}
				if (m_draw_axis)
				{
					draw_XYZ_Normals();
				}
			}
		}
		return true;
	}

	void Another_HDS_model::enable_material()
	{
		glEnable(GL_COLOR_MATERIAL);
	}

	/************************************************************************/
	/* /brief Load model from OFF file                                      */
	/*																		*/
	/************************************************************************/
	bool Another_HDS_model::load_off_file(const char * filename)
	{
		//  		 		ifstream stream;
		//  				stream.open(filename);
		//  				stream >> this;
		// 		std::cout << "loadMesh...  "<< "Polysurf with " << m_original_model.size_of_vertices()
		// 			<< " vertices and " << m_original_model.size_of_facets()
		// 			<< " facets. " << std::endl;
		return true;
	}

	/************************************************************************/
	/*  /brief Load model from M file                                       */
	/************************************************************************/
	bool Another_HDS_model::load_m_file(FILE* mf)
	{	
		if (mf != NULL)
		{
			m_triangle.SetFile(mf);
			m_triangle.SetFileType(MFILE);
			//	m_original_model.delegate( m_triangle );
			delegate( m_triangle );
			//CGAL_assertion( m_original_model.is_triangle( m_original_model.halfedges_begin()));
			return true;
		}
		else
			return false;
	}

	/************************************************************************/
	/*  /brief Load model from wavefront obj file                                       */
	/************************************************************************/
	bool Another_HDS_model::load_obj_file(FILE* mf)
	{	
		if (mf != NULL)
		{
			m_triangle.SetFile(mf);
			m_triangle.SetFileType(OBJFILE);
			//	m_original_model.delegate( m_triangle );
			delegate( m_triangle );

			//CGAL_assertion( m_original_model.is_triangle( m_original_model.halfedges_begin()));
			return true;
		}
		else
			return false;
	}

	//Draw x,y,z axis
	void Another_HDS_model::draw_XYZ_Normals()
	{
		// Z axis
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		GLUquadricObj* quadObj1;
		quadObj1 = gluNewQuadric();
		gluQuadricDrawStyle(quadObj1,GLU_FILL);
		gluQuadricNormals(quadObj1,GL_FLAT);
		gluQuadricOrientation(quadObj1,GLU_OUTSIDE);
		gluQuadricTexture(quadObj1,GL_TRUE);
		glColor3f(0.0,0.0,1.0);
		gluCylinder(quadObj1,0.01,0.01,1,20,8);

		glPushMatrix();
		glTranslatef(0,0,1);
		GLUquadricObj* quadObjcone1;
		quadObjcone1 = gluNewQuadric();
		gluQuadricDrawStyle(quadObjcone1,GLU_FILL);
		gluQuadricNormals(quadObjcone1,GL_FLAT);
		gluQuadricOrientation(quadObjcone1,GLU_OUTSIDE);
		gluQuadricTexture(quadObjcone1,GL_TRUE);
		glColor3f(0.0,0.0,1.0);
		gluCylinder(quadObjcone1,0.03,0,0.2,20,8);
		glPopMatrix();
		gluDeleteQuadric(quadObjcone1);


		// X axis

		glRotatef(90,0,1,0);
		GLUquadricObj* quadObj2;
		quadObj2 = gluNewQuadric();
		gluQuadricDrawStyle(quadObj2,GLU_FILL);
		gluQuadricNormals(quadObj2,GL_FLAT);
		gluQuadricOrientation(quadObj2,GLU_OUTSIDE);
		gluQuadricTexture(quadObj2,GL_TRUE);
		glColor3f(1.0,0.0,0.0);
		gluCylinder(quadObj2,0.01,0.01,1,20,8);
		gluDeleteQuadric(quadObj2);

		glPushMatrix();
		glTranslatef(0.0,0.0,1.0);
		GLUquadricObj* quadObjcone2;
		quadObjcone2 = gluNewQuadric();
		gluQuadricDrawStyle(quadObjcone2,GLU_FILL);
		gluQuadricNormals(quadObjcone2,GL_FLAT);
		gluQuadricOrientation(quadObjcone2,GLU_OUTSIDE);
		gluQuadricTexture(quadObjcone2,GL_TRUE);
		glColor3f(1.0,0.0,0.0);
		gluCylinder(quadObjcone2,0.03,0,0.2,20,8);
		glPopMatrix();
		gluDeleteQuadric(quadObjcone2);

		glRotatef(-90,1,0,0);
		// Y axis
		GLUquadricObj* quadObj3;
		quadObj3 = gluNewQuadric();
		gluQuadricDrawStyle(quadObj3,GLU_FILL);
		gluQuadricNormals(quadObj3,GL_FLAT);
		gluQuadricOrientation(quadObj3,GLU_OUTSIDE);
		gluQuadricTexture(quadObj3,GL_TRUE);
		glColor3f(0.0,1.0,0.0);
		gluCylinder(quadObj3,0.01,0.01,1,20,8);
		gluDeleteQuadric(quadObj3);

		glPushMatrix();
		glTranslatef(0,0,1);
		GLUquadricObj* quadObjcone3;
		quadObjcone3 = gluNewQuadric();
		gluQuadricDrawStyle(quadObjcone3,GLU_FILL);
		gluQuadricNormals(quadObjcone3,GL_FLAT);
		gluQuadricOrientation(quadObjcone3,GLU_OUTSIDE);
		gluQuadricTexture(quadObjcone3,GL_TRUE);
		glColor3f(0.0,1.0,0.0);
		gluCylinder(quadObjcone3,0.03,0,0.2,20,8);
		glPopMatrix();
		glPopMatrix();
		gluDeleteQuadric(quadObjcone3);
	}


	// Read UV for each vertex from a file
	bool Another_HDS_model::load_UV_file(FILE* uvfile)
	{
		if (uvfile!=NULL)
		{	
			// Read uv coordinates number from file
			int UVnum ;
			if( fscanf(uvfile,"uvnum %d",&UVnum) )
			{
				//Read uv coordinates from file
				int tmp_uvindex;
				float tmp_u;
				float tmp_v;
				float* uArray = (float*)malloc( sizeof(float)*UVnum );
				float* vArray = (float*)malloc( sizeof(float)*UVnum );
				if (uArray!=NULL &&vArray!=NULL)
				{
					for (int i = 0; i< UVnum; i++)
					{
						if( fscanf(uvfile,"%d	(%f, %f)",&tmp_uvindex, &tmp_v,&tmp_u) )
						{
							uArray[tmp_uvindex] = tmp_u;
							vArray[tmp_uvindex] = tmp_v;
						}
					}

					//Assign uv to points
					int faceNum;
					int tmp_faceindex;
					int tmp_facevertexindex, tmp_vertexindex, tmp_normalindex, tmp_colorindex, tmp_uvmapindex,tmp_uvmap1;
					Vertex_iterator tmp_v;
					fscanf(uvfile,"%d",&faceNum);
					for (int j = 0; j<faceNum; j++)
					{
						//Format:  Face|faceVertexIndex|vertexIndex|normalIndex|colorIndex|| UV_map1
						for (int k = 0; k<3; k++)
						{
							//pig
							fscanf(uvfile,"%d	%d	%d	%d	%d		%d",&tmp_faceindex, &tmp_facevertexindex,&tmp_vertexindex,&tmp_normalindex,&tmp_colorindex,&tmp_uvmapindex);

							//bunny
							//	fscanf(uvfile,"%d	%d	%d	%d	%d		%d	%d",&tmp_faceindex, &tmp_facevertexindex,&tmp_vertexindex,&tmp_normalindex,&tmp_colorindex,&tmp_uvmap1,&tmp_uvmapindex);
							//m file index from 1 , uv file index from 0
							tmp_v = m_vertexMap.FindVertexhandleO(tmp_vertexindex+1);
							tmp_v->SetUV( uArray[tmp_uvmapindex], vArray[tmp_uvmapindex] );
						}
						fscanf(uvfile,"\n\r");

					}
					//	fscanf(uvfile,"")
					//free uv arrays

					free(uArray);
					free(vArray);
				}
			}
			else
			{
				return false;
			}

			m_has_UV = true;
			return true;
		}
		else
		{
			return false;
		}

	}

	//Set a one to one map between index in m file and vertex handle
	bool Another_HDS_model::setup_OneToOneVertexMap()
	{
		int tmp_index = 1;
		for( Vertex_iterator v = vertices_begin(); v != vertices_end();++v)
		{
			v->SetIndex(tmp_index);
			m_vertexMap.InsertvertexHandleO(tmp_index, v);
			tmp_index ++;
		}
		return true;
	}

	void Another_HDS_model::set_correspondence(const vector<int>& indexs)
	{
		m_index.clear();
		m_enlarge.clear();
		for (size_t i =0; i<indexs.size(); i++)
		{
			m_index.push_back(indexs[i]);
			m_enlarge.push_back(false);
		}
	}






	//     void Another_HDS_model::DrawCursorO()
	//         {
	//         glPushMatrix();
	//         glTranslatef(my_cursorO.x(), my_cursorO.y(), my_cursorO.z() );
	//         glColor3f(1,0,1);
	//         glutSolidSphere(0.01,16,16);
	//         glPopMatrix();
	//         }

	//     void Another_HDS_model::DrawCursorP()
	//         {
	//         glPushMatrix();
	//         glTranslatef(my_cursorP.x(), my_cursorP.y(), my_cursorP.z() );
	//         glColor3f(1,0,1);
	//         glutSolidSphere(0.01,16,16);
	//         glPopMatrix();
	//         }


	//     void Another_HDS_model::DrawGrid()
	//         {
	//         // draw grid
	//         glMatrixMode(GL_MODELVIEW);
	//         glPushMatrix();
	//         glBegin(GL_LINES);
	//         glColor3f(1,1,1);
	//         for(int i=-10;i<=10;++i) 
	//             {
	//             glVertex3f(i,0,-10);
	//             glVertex3f(i,0,10);
	// 
	//             glVertex3f(10,0,i);
	//             glVertex3f(-10,0,i);
	//             }
	//         glEnd();
	//         glPopMatrix();
	//         }
	bool Another_HDS_model::set_U(double* u, int m, int n)
	{
		for( Vertex_iterator v = vertices_begin(); v != vertices_end();++v)
		{
			int tmpIndex  = v->GetIndex()-1;
			v->SetU(u[tmpIndex]);
		}
		return true;
	}

	bool Another_HDS_model::PlanarEmbedding()
	{
		Halfedge_iterator root = edges_begin();

		root->SetConjU(0.0);
		double u1 = root->vertex()->GetU();
		double u2 = root->opposite()->vertex()->GetU();
		root->SetU( (u1+u2)/2.0 );
		

		Init_Faces_Unvisited();
		BFS_faces(root,true);
		return true;
	}

	bool Another_HDS_model::Init_Faces_Unvisited()
	{
		for (Face_iterator it =  facets_begin(); it != facets_end(); it++ )
		{
			it->visited = false;
		}
		return true;
	}
	void Another_HDS_model::BFS_faces(Halfedge_handle root,bool visited)
	{
		queue<Halfedge_handle> tmp_q;
		tmp_q.push(root);  //push root
		root->face()->visited = visited;


		//Processing the queue
		while(tmp_q.size()>0)
		{
			Halfedge_handle tmpEdge = (Halfedge_handle)tmp_q.front();
			tmp_q.pop();
			int size = tmp_q.size();
			Face_handle tmpFace = tmpEdge->face();
			caculate_conjuct_facet(tmpEdge);

			//push all unvisited the childs
			bool first = true;
			Halfedge_handle tmp_child_Edge = tmpEdge->opposite();
			Face_handle  tmp_child_Face  = tmp_child_Edge->face();
			Face_handle firstFace = tmp_child_Face;
			while(tmp_child_Face!=firstFace||first)
			{
				first = false;
				if(tmp_child_Face->visited == visited)
				{
					//search for next sibling
					tmpEdge = tmpEdge->next();
					tmp_child_Edge = tmpEdge->opposite();
					tmp_child_Face = tmp_child_Edge->face();
					continue;
				}
				else
				{
					//calculate child
					tmp_q.push(tmp_child_Edge);
					tmp_child_Face->visited = visited;

					//search for next sibling
					tmpEdge = tmpEdge->next();
					tmp_child_Edge = tmpEdge->opposite();
					tmp_child_Face = tmp_child_Edge->face();

					continue;
				}
			}
		}
	}

	void Another_HDS_model::caculate_conjuct_facet(Halfedge_handle seed)
	{
		Halfedge_handle Ejk = seed;
		double Us = Ejk->GetConjU(); 
		double Ur, Ut;
		double Uk, Ui, Uj;
		Face_handle tmp_face = Ejk->face();

		//Caculate the conjuct for this face
		Vertex_handle Vk = Ejk->vertex();
		Uk = Vk->GetU();
		Point Vk_position;
		Vk_position = Vk->point();

		Halfedge_handle Eki = Ejk->next();
		Vertex_handle Vi = Eki->vertex();
		Ui = Vi->GetU();
		Point Vi_position;
		Vi_position = Vi->point();

		Halfedge_handle Eij = Eki->next();
		Vertex_handle Vj = Eij->vertex();
		Uj = Vj->GetU();
		Point Vj_position;
		Vj_position = Vj->point();

		double cotuij = (Ui- Uj)/std::tan( compute_angle_rad(Vj_position,Vk_position, Vi_position) );
		double cotukj = (Uk- Uj)/std::tan( compute_angle_rad(Vj_position, Vi_position, Vk_position) );
		Ur = Us + ( cotuij + cotukj)/(double)2;

		Ut = Ur + ( (Uk- Ui)/std::tan(compute_angle_rad(Vi_position,Vj_position, Vk_position)) + (Uj-Ui)/std::tan(compute_angle_rad(Vj_position,Vk_position, Vi_position)))/(double)2;

		Eij->SetU((Ui+Uj)/(double)2);
		Eij->SetConjU(Ur);
		Eij->opposite()->SetConjU(Ur);

		Eki->SetU((Ui+Uk)/(double)2);
		Eki->SetConjU(Ut);
		Eki->opposite()->SetConjU(Ut);

		Ejk->SetU((Uj+Uk)/(double)2);
		Ejk->opposite()->SetConjU(Us);
	}

	/************************************************************************/
	/* Set feature points on mesh                                           */
	/************************************************************************/
	void Another_HDS_model::set_featurePoint(const std::vector<Vertex_handle>& features)
	{
		assert(features.size());
		m_features.clear();
		for (size_t i = 0; i < features.size(); i++)
		{
			m_features.push_back(features[i]);
		}
	}

	/************************************************************************/
	/* Set mid-edge feature points on mesh                                  */
	/************************************************************************/
	void Another_HDS_model::set_mid_edge_features(const std::vector<complex<double>>& features)
	{
		assert(features.size());
		m_middle_edge_features.clear();
		m_middle_edge_features_select.clear();
		for (size_t i = 0; i < features.size(); i++)
		{
			m_middle_edge_features.push_back(features[i]);
			m_middle_edge_features_select.push_back(false);
		}
	}

	/************************************************************************/
	/* /brief calculate the halfedge length for each edge                   */
	/************************************************************************/
	void Another_HDS_model::calculate_edge_length()
	{
		Halfedge_iterator itb = halfedges_begin(), ite = halfedges_end();
		for(; itb!=ite; itb++) 
		{
			Halfedge_handle h = itb;
			double l = Edge_length()(h);
			itb->m_len = l;
			itb->visited = false;
		}
	}


	void Another_HDS_model::set_enlarge(GLint name)
	{
		if (m_enlarge.size()!=0 &&  name >-2)
		{
			if (name == -1)
			{
				for (int i = 0; i<m_enlarge.size(); i++)
				{
					m_enlarge[i] = false;
				}
			}
			else
			{
				m_enlarge[name] = true;
			}
		}
		else if (m_enlarge.size() == 0 )
		{
			cerr <<"Error: Enlarge size is zero"<<endl;
		}
		else
		{
			cerr <<"Error: name is out range of Enlarge"<<endl;
		}
		
	}


	/************************************************************************/
	/* /brief Draw feature points                                           */
	/************************************************************************/
	void Another_HDS_model::draw_feature_points(int mode)
	{
		switch(mode)
		{
		case ALIGNMENT_SELECT:
			if (m_features.size() == 0)
			{
				cout<<"No feature point on mesh"<<endl;
				break;
			}
			for (int i = 0; i<m_features.size(); i++)
			{
				glShadeModel(GL_SMOOTH);
				glPushMatrix();
				Vertex_handle v = m_features[i];
				Point_3 v_position = v->point();
				glTranslatef(v_position.x(),v_position.y(),v_position.z());
				
			
					float color_tmp[3];
					color_tmp[0] = 0.0;
					color_tmp[1] = 0.0;
					color_tmp[2] = 0.0;
					bool selected= false;
					assert(m_three_points.size()<4);
					for(int j =0; j <m_three_points.size(); j++ )
					{
						if (i == m_three_points[j])
						{
							selected = true;
							color_tmp[j] = 1.0;
						}
					}
					
					if (selected ==false)
					{
						color_tmp[0] = 1.0;
						color_tmp[1] = 1.0;
						color_tmp[2] = 0.0;
					}
					glColor3f(color_tmp[0], color_tmp[1], color_tmp[2]);
				
				glutSolidSphere(0.5,16,16);
				glPopMatrix();
			}
			break;
		case GL_SMOOTH:
			if (m_features.size() == 0)
			{
				cout<<"No feature point on mesh"<<endl;
				break;
			}
			for (int i = 0; i<m_features.size(); i++)
			{
				Vertex_handle v = m_features[i];
				Point_3 v_position = v->point();
				glShadeModel(GL_SMOOTH);
				glPushMatrix();
				glTranslatef(v_position.x(),v_position.y(),v_position.z());
				glColor3f(1,0,1);
				glutSolidSphere(0.5,16,16);
				glPopMatrix();
			}
			break;

		case  GL_SELECT:
			if (m_features.size() == 0)
			{
				cout<<"No feature point on mesh"<<endl;
				break;
			}
			for (int i = 0; i<m_features.size(); i++)
			{
				Vertex_handle v = m_features[i];
				Point_3 v_position = v->point();
				glShadeModel(GL_SELECT);
				glPushMatrix();
				glTranslatef(v_position.x(),v_position.y(),v_position.z());
				glPushName(i);
				glutSolidSphere(0.5,16,16);
				glPopName();
				glPopMatrix();
			}
			break;

		case FEATURE_CORRESPONDENCE:
			if ( m_index.size() ==0 )
			{
				cout<<"No matching"<<endl;
				break;
			}
			else
			{
				m_color_correspoendence = 50;
				for (int i = 0; i<m_index.size(); i++)
				{
					Vertex_handle v = m_features[m_index[i]];
					Point_3 v_position = v->point();

					if (mode == GL_SELECT)
					{
						glShadeModel(GL_SELECT);
					}
					else
					{
						glShadeModel(GL_SMOOTH);
					}

					glPushMatrix();
					glTranslatef(v_position.x(),v_position.y(),v_position.z());


					if (mode == GL_SELECT)
					{
						glPushName(i);
					}

					GLfloat red = get_color_map(i%6,0);  //pt color
					GLfloat green = get_color_map(i%6,1);  //pt color
					GLfloat blue = get_color_map(i%6,2);  //pt color
					glColor3f(red,green,blue);

					if (m_enlarge[i])
					{
						glutSolidSphere(0.8,16,16);
					}
					else
					{
						glutSolidSphere(0.5,16,16);
					}
					if (mode == GL_SELECT)
					{
						glPopName();
					}
					glPopMatrix();
				}
			}
			break;
		}
	}

	void Another_HDS_model::init_color_map(int mapSize, int min_bright, int max_dark)
	{
		int step = 45;
		int m_color_correspoendence = 50;
		for (int k =0; k<mapSize; k++)
		{
// 			int t1,t2,t3,tmp;
// 			for(;;)
// 			{
// 				if (m_color_correspoendence+ step>215)
// 					m_color_correspoendence=m_color_correspoendence - 216 + step;
// 				else
// 					m_color_correspoendence+=step;
// 				tmp = m_color_correspoendence;
// 				t1 = (int)tmp/36;
// 				t2 = (int)(tmp%36) /6;
// 				t3 = (tmp%36) %6;
// 
// 				if ((t1+t2+t3 >min_bright) && (t1+t2+t3 < max_dark))
// 				{                    
// 					break;
// 				}
// 			}
// 
// 			const int k_max = 255;
// 
// 			float m1 = (float)((k_max -t1*51) / 256.0);
// 			tmp = tmp % 36;
// 			float m2 = (float)((k_max -t2*51) / 256.0);
// 			float m3 = (float)((k_max -t3*51) / 256.0);

			float m1 = random(0.01f, 1.0f);
			float m2 = random(0.01f, 1.0f);
			float m3 = random(0.01f, 1.0f);
			
			float sum = m1*m1 + m2*m2 + m3*m3;
			sum = sqrt(sum);
			if (sum > color_epsilon )
			{
				m1 = m1/sum;
				m2 = m2/sum;
				m3 = m3/sum;
			}

			g_color_map[k][0] = m1;
			g_color_map[k][1] = m2;
			g_color_map[k][2] = m3;
		}
	}
	////////////////////////////////////////////////////////////////////////////////
	//
	//              Computes a random float value in a specified range
	//
	////////////////////////////////////////////////////////////////////////////////
	float Another_HDS_model::random(float min, float max)
	{
		assert(min < max);  // This isn't really necessary

		return(min + (max - min) * (float(rand() % 1024) / 1024.0f));
	}

	GLfloat Another_HDS_model::get_color_map(int index, int component)
	{
		return g_color_map[index][component];
	}
	

	void Another_HDS_model::init_features_unselected()
	{
		for (int i = 0; i<m_features.size(); i++)
		{
			m_features_selected.push_back(false);
		}
	}
	void Another_HDS_model::set_feature_visited(int name)
	{
		if (m_features_selected.size()>name)
		{
			m_features_selected[name] = true;
		}
		else
		{
			cerr<<"Selected feature point index is out of rang"<<endl;
		}
	}

	void Another_HDS_model::set_feature_unvisited(int name)
	{
		if (m_features_selected.size()>name)
		{
			m_features_selected[name] = false;
		}
		else
		{
			cerr<<"Selected feature point index is out of rang"<<endl;
		}
	}

	/************************************************************************/
	/* /brief Get the feature point vertex index by name                    */
	/* /param[in] int the index of feature point in m_features				*/
	/* /param[out] int the vertex index of feature point in HDS				*/
	/************************************************************************/
	int Another_HDS_model::get_feature_index(int name)
	{
		Vertex_handle v = m_features[name];
		return v->GetIndex();
	}


	bool Another_HDS_model::load_laplacian_mesh(ifstream& file)
	{
		vector<Point> coordinates;
		float v1,v2,v3;
		while(file.good())
		{
			file>>v1>>v2>>v3;
			Point tmpp(v1,v2,v3);
			coordinates.push_back(tmpp);
		}

		Vertex_handle tmpv = vertices_begin();
		for ( ; tmpv!= vertices_end(); ++tmpv )
		{
			int index = tmpv->GetIndex()-1;
			tmpv->m_laplace_position = coordinates[index];
		}

		return true;
	}


	/************************************************************************/
	/* /brife                                                                     */
	/************************************************************************/
	bool Another_HDS_model::save_m_file(string filename, int type)
	{
		//Open a file and check the statu of file
		ofstream stream;
		stream.open(filename.c_str());
		if (stream.fail()||!stream.good())
		{
			cerr<<"Error Another_HDS_model::save_m_file: Can not open file: "<<filename<<endl;
			return false;
		}

		//Write vertices into file. the index is 1-based
		Vertex_iterator vib = vertices_begin();
		Vertex_iterator vie = vertices_end();
		Vertex_iterator vitmp = vib;
		for (; vitmp!=vie; vitmp++ )
		{
			//Vertex 1  -0.998566 -0.998968 0.882 
			int vertexIndexTmp = vitmp->GetIndex();
			Point_3 coordinatesTmp;
			if (type == 1)
			{
				coordinatesTmp = vitmp->m_laplace_position;
			}
			else
				 coordinatesTmp = vitmp->point();

			stream<<"Vertex "<<vertexIndexTmp<<" "<<coordinatesTmp[0]<<" "<<coordinatesTmp[1]<<" "<<coordinatesTmp[2]<<endl;
		}

		//Write faces into file, the index is 1-based
		Face_iterator fib = facets_begin();
		Face_iterator fie = facets_end();
		Face_iterator fitemp = fib;
		int faceIndexTmp = 1;
		for (; fitemp!=fie; fitemp++)
		{
			//Format: Face 1 100 101 102

			int vertexIndex[3];
			Halfedge_iterator tmpedge = fitemp->halfedge();
			Vertex_iterator v = tmpedge->vertex();
			vertexIndex[0] = v->GetIndex();
			tmpedge=tmpedge->next();
			v = tmpedge->vertex();
			vertexIndex[1] = v->GetIndex();
			tmpedge=tmpedge->next();
			v = tmpedge->vertex();
			vertexIndex[2] = v->GetIndex();

			stream<<"Face "<<faceIndexTmp<<" "<<vertexIndex[0]<<" "<<vertexIndex[1]<<" "<<vertexIndex[2]<<endl;

			faceIndexTmp++;
		}

		stream.close();
		return true;
	}

	/************************************************************************/
	/* /brife                                                                     */
	/************************************************************************/
	bool Another_HDS_model::save_m_file(string filename)
	{
		//Open a file and check the statu of file
		ofstream stream;
		stream.open(filename.c_str());
		if (stream.fail()||!stream.good())
		{
			cerr<<"Error Another_HDS_model::save_m_file: Can not open file: "<<filename<<endl;
			return false;
		}

		//Write vertices into file. the index is 1-based
		Vertex_iterator vib = vertices_begin();
		Vertex_iterator vie = vertices_end();
		Vertex_iterator vitmp = vib;
		for (; vitmp!=vie; vitmp++ )
		{
			//Vertex 1  -0.998566 -0.998968 0.882 
			int vertexIndexTmp = vitmp->GetIndex();
			Point_3 coordinatesTmp = vitmp->point();
			stream<<"Vertex "<<vertexIndexTmp<<" "<<coordinatesTmp[0]<<" "<<coordinatesTmp[1]<<" "<<coordinatesTmp[2]<<endl;
		}

		//Write faces into file, the index is 1-based
		Face_iterator fib = facets_begin();
		Face_iterator fie = facets_end();
		Face_iterator fitemp = fib;
		int faceIndexTmp = 1;
		for (; fitemp!=fie; fitemp++)
		{
			//Format: Face 1 100 101 102

			int vertexIndex[3];
			Halfedge_iterator tmpedge = fitemp->halfedge();
			Vertex_iterator v = tmpedge->vertex();
			vertexIndex[0] = v->GetIndex();
			tmpedge=tmpedge->next();
			v = tmpedge->vertex();
			vertexIndex[1] = v->GetIndex();
			tmpedge=tmpedge->next();
			v = tmpedge->vertex();
			vertexIndex[2] = v->GetIndex();

			stream<<"Face "<<faceIndexTmp<<" "<<vertexIndex[0]<<" "<<vertexIndex[1]<<" "<<vertexIndex[2]<<endl;

			faceIndexTmp++;
		}
	
		stream.close();
		return true;
	}
}