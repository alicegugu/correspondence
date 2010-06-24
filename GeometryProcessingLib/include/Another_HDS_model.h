#pragma once

#include <CGAL/basic.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream_ture.h>
#include <iostream>

#include <nvMath.h>

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <stdio.h>
#include <stdlib.h>

#include "Traits.h"
#include "My_Vertex.h"
#include "My_Halfedge.h"
#include "My_face.h"
#include "Build_triangle.h"



using namespace std;

namespace Another
{
	class Surface_point;

	enum RenderOption
	{
		DISCRETE_HARMONIC_MAP = 0,//show harminic map
		PLANAR_EMBEDDING = 1, //show embedded
		MESH = 2,
		CURVATURE =3,
		TRIANGULAR = 4, //show triangular mesh
		FEATURE_POINT = 5,//select 
		CORRESPONDENCE = 6,
		DEFAULT = 7,
		LAPLACE = 8,
		CUTFACESELECTION = 9,
		ALIGNMENT_SELECT = 10,
		ALIGNMENT_RESULT = 11
	};
    //An item wrapper using My_vertex, My_halfedge, My_face
	struct My_items: public CGAL::Polyhedron_items_3
    {
    
        template<class Refs, class Traits>
        struct Vertex_wrapper
        {
			typedef struct 
			{
			public:
				typedef typename Traits::Point_3 Point_3;
				typedef typename Traits::Point_2 Point_2;
				typedef typename Traits::Vector_3 Vector_3;
			} FGeomTraits;
			typedef typename Traits::Point_3 Point_3;
			typedef My_vertex < Refs, CGAL::Tag_true, Point_3, FGeomTraits > Vertex;
        };
        
        template<class Refs, class Traits>
        struct Halfedge_wrapper
        {
			typedef struct 
			{
			public:
				typedef typename Traits::Point_3 Point_3;
				typedef typename Traits::Point_2 Point_2;
				typedef typename Traits::Vector_3 Vector_3;
			} FGeomTraits;
            typedef My_halfedge<Refs, CGAL::Tag_true,CGAL::Tag_true,CGAL::Tag_true,FGeomTraits >   Halfedge;
        };
        
        template<class Refs, class Traits>
        struct Face_wrapper
        {
			//typedef typename Traits::Vector_3 Vector_3;
			//all types needed by the facet...
			typedef struct {
			public:
				typedef typename Traits::Vector_3 Vector_3;
			} FGeomTraits;
			//custom type instantiated...
            typedef My_face<Refs, FGeomTraits>       Face;
        };
        
   };


    typedef CGAL::Polyhedron_3< Another_Traits, My_items>      Another_HDS;
    typedef CGAL::HalfedgeDS_decorator<Another_HDS>  Decorator;
    typedef Another_HDS::Vertex	                 Vertex;
    typedef Another_HDS::Face	         Face;
    typedef Another_HDS::Face_handle     Face_handle;
    typedef Another_HDS::Vertex_handle	Vertex_handle;
    typedef Another_HDS::Halfedge		Edge;
    typedef Another_HDS::Halfedge_handle    Halfedge_handle;
    typedef Another_HDS::Vertex_iterator	Vertex_iterator;
    typedef Another_HDS::Halfedge_iterator 	Halfedge_iterator;
    typedef Another_HDS::Face_iterator	Face_iterator;
    typedef Another_HDS::HalfedgeDS      HalfedgeDS ;
	typedef Another_HDS::Point_3 Point_3;
	typedef Another_HDS::Facet Facet;

	typedef Sort<Another::Vertex_iterator> SORT_PAINTMODEL_VERTEX;


	class Another_HDS_model:public Another_HDS
    {
        private:
        
            //Color 
            float m_color[4];
           //Material 
            float m_ambient_matrix[4];
            float m_diffuse_matrix[4];
            float m_shininess;
            float m_specular_matrix[4];
            
           //Model's geometry info 
            int m_corner_count;
            int m_edge_count;
            int m_facet_count;
            
            int m_mode;
            
            //Object space coordinates
            Vector_3 m_axis_x; 
            Vector_3 m_axis_y ;
            Vector_3 m_axis_z;
            Point m_original;
            
           //View mode control 
            bool m_draw_axis;
            int m_draw_type; 
			
			My_One2OneMap<unsigned int,Another::Vertex_iterator,SORT_PAINTMODEL_VERTEX> m_vertexMap;

			//Feature points
			vector<Vertex_handle> m_features;
			vector<bool> m_features_selected;
			vector<int> m_index;
			vector<bool> m_enlarge;
			vector<complex<double>> m_middle_edge_features;
			vector<bool> m_middle_edge_features_select;

            //Load mesh
            T_Build_triangle<HalfedgeDS> m_triangle;

			//correspondence
			map<int, int> m_correspondence;
	public:
        //    Another_HDS m_original_model; 

            bool m_has_UV;
			bool m_has_gaussian_curvature;
			bool m_has_1D_u;
			int m_color_correspoendence;

	public:
			int m_cutface; //Face index of cut face

			vector<int> m_three_points;

        public:
        
            Another_HDS_model(void);
            ~Another_HDS_model(void);
             
  			bool set_curvature() 
			{
				m_has_gaussian_curvature = true; 
				return m_has_gaussian_curvature;
			}
			bool un_set_curvature()
			{
				m_has_gaussian_curvature = false; 
				return m_has_gaussian_curvature;
			}

			vector<Vertex_handle>& get_samples()
			{
				return m_features;
			}

		//	void insert_correspondence( int name);
           
           /******************************************Render Model and Polycube *********************************************/ 
           //main draw function, render content based on mode and type control
            bool draw(int type,int PolyOrOriginal,int debugmode);

            //Draw object space coordinates
            void draw_XYZ_Normals();
			bool draw_midedge();
			void draw_HamonicMap();
			bool draw(RenderOption renderOption);
			void draw_feature_points(int mode);
			void set_enlarge(GLint name);
			void draw_mid_edge_features(int mode);
			void draw_laplacian();
           /*************************************************************************************************************************************/
            
           
           
           
           
            
            /***************************************Set and get attributes*********************************************************************/
           
            bool set_draw_type(int t) { m_draw_type = t;return true; }
            int get_draw_type() { return m_draw_type; }
            void set_draw_axis( bool dr ) { m_draw_axis = dr; } 
            bool get_draw_axis( ) { return m_draw_axis; }
            bool set_color(float* color);
            float* get_color( ) { return m_color; }
            void enable_material();
			void get_guassian_curvature(vector<double>& gc);
			bool set_1D_texture_coordinates(vector<double>& texture_coordinates);
		//	void set_translation_matrix_screen_to_objspace(QMatrix translation_matrix);
			void get_sample_indexes(vector<int>& features) 
			{
				features.clear();
				for (int i = 0; i<m_features.size(); ++i)
				{
					Vertex_handle Vtmp = m_features[i];
					features.push_back(Vtmp->GetIndex());
				} 
			}
            
           /****************************************************************************************************************************************/
           
           
           /************************************************File access************************************************************************/
            //Load model and polycube from m file
            bool load_m_file(FILE * mf);
			bool load_off_file(const char * filename);
			bool load_obj_file(FILE* mf);
			bool load_laplacian_mesh(ifstream& file);
            
           //Save all the setting and model data to a file 
            bool Save();
           
            // Read UV for each vertex from a file
            bool load_UV_file(FILE* uvfile);
            
            //Set up a one to one map between index and vertex in original model
            bool setup_OneToOneVertexMap();

			bool Init_Faces_Unvisited();

			//Export the model into a m file
			bool save_m_file(string filename);
			bool Another_HDS_model::save_m_file(string filename, int type);
           
           /**************************************************************************************************************************************/ 

            //Calculate plane type
            //0: x-positive, 1:x-negative, 2: y-positivem 3:y-negative, 4:z-positive, 5:z-negative
            int get_root_facet_type();

			bool set_U(double* u, int m, int n);
			bool PlanarEmbedding();
			void caculate_conjuct_facet(Halfedge_handle seed);
			void BFS_faces(Halfedge_handle root,bool visited);
             
			/************************************************Geometry Processing*********************************************/
			void Sampling();

			int closest_vertices(Surface_point& p, std::vector<Vertex_handle>& storage);

			void set_featurePoint(const std::vector<Vertex_handle>& features);

			void calculate_edge_length();

			//Another_HDS::Iso_cuboid_3 get_bounding_box();

			void set_correspondence(const vector<int>& indexs);
			static 	void init_color_map(int mapSize = 100, int min_bright=4, int max_dark=12);
			static GLfloat get_color_map(int index, int component);
			static float random(float min, float max);

			int get_cut_face() { return m_cutface; }
			void set_cut_face(int face) { m_cutface = face; }
			void set_mid_edge_features(const std::vector<complex<double>>& features);

			//Feature selections
			void init_features_unselected();
			void set_feature_visited(int name);
			void set_feature_unvisited(int name);

			//Get the feature point vertex index by name
			int get_feature_index(int name);

			bool add_three_points(int featureindex)
			{
				m_three_points.push_back(featureindex);
				if (m_three_points.size()>3)
				{
					return false;
				}
				return true;
			}
        };

    }
