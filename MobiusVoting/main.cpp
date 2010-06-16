///////////////////////////////////////////////////////////////////////////////////
//
// Author: GU Huiqin Date: 2010 Feb Title: Mobius Voting for Surface Correspondence
//
///////////////////////////////////////////////////////////////////////////////////

#include "headers.h"


namespace po = boost::program_options;


#define CHECKER_BOX "./textures/checker_medium.gif"

#define ITEM_1   1
#define ITEM_2   2


GLuint v,f,f2,p;
float lpos[4] = {1,0.5,1,0};
vector<Another::Another_HDS_model*> g_model;
int sourceModel = 0;
int targetModel = 0;
int g_current_model = 0;
int sample_size = 300;

//vertex list for triangulation
vector<int> sourceFeatures;
vector<int> targetFeatures;

void drawPickSource(void)
{
	/* Name-stack manipulation for the purpose of
       selection hit processing when mouse button
       is pressed.  Names are ignored in normal
       OpenGL rendering mode.       */
//	glViewport(0, 0, 1920/2, 1080);
	g_model[sourceModel]->draw_feature_points(GL_SELECT);

}
void drawPickTarget(void)
{
	/* Name-stack manipulation for the purpose of
       selection hit processing when mouse button
       is pressed.  Names are ignored in normal
       OpenGL rendering mode.       */
//	glViewport(0, 0, 1920/2, 1080);
	g_model[targetModel]->draw_feature_points(GL_SELECT);
}
void ProcessPickUpSource(GLint name)
{
// 	if (name >-1)
// 	{
// 		g_model[0]->set_enlarge(name);
// 		g_model[1]->set_enlarge(name);
// 	}
// 	else
// 	{
// 		g_model[0]->set_enlarge(name);
// 		g_model[1]->set_enlarge(name);
// 	}
	if (name>-1)
	{
		int indexTmp = g_model[sourceModel]->get_feature_index(name) ;
		cout<<"Pick up source model's feature point: "<<indexTmp<<endl;
		//g_model[sourceModel]->insert_correspondence(name);
		sourceFeatures.push_back( indexTmp);
	}
}

void ProcessPickUpTarget(GLint name)
{
// 	if (name >-1)
// 	{
// 		g_model[0]->set_enlarge(name);
// 		g_model[1]->set_enlarge(name);
// 	}
// 	else
// 	{
// 		g_model[0]->set_enlarge(name);
// 		g_model[1]->set_enlarge(name);
// 	}
	if (name>-1)
	{
		int indexTmp = g_model[targetModel]->get_feature_index(name);
		cout<<"Pick up target model's feature point: "<<indexTmp<<endl;
		//g_model[sourceModel]->insert_correspondence(name);
		targetFeatures.push_back( indexTmp );
	}
}

void init_color_map(int mapSize = 6, int min_bright=4, int max_dark=12)
{
	int step = 45;
	int m_color_correspoendence = 50;
	for (int i =0; i<mapSize; i++)
	{
		int t1,t2,t3,tmp;
		for(;;)
		{
			if (m_color_correspoendence+ step>215)
				m_color_correspoendence=m_color_correspoendence - 216 + step;
			else
				m_color_correspoendence+=step;
			tmp = m_color_correspoendence;
			t1 = (int)tmp/36;
			t2 = (int)(tmp%36) /6;
			t3 = (tmp%36) %6;

			if ((t1+t2+t3 >min_bright) && (t1+t2+t3 < max_dark))
			{                    
				break;
			}
		}

		const int k_max = 255;

		float m1 = (float)((k_max -t1*51) / 256.0);
		tmp = tmp % 36;
		float m2 = (float)((k_max -t2*51) / 256.0);
		float m3 = (float)((k_max -t3*51) / 256.0);

		//glColor4f(m1,m2,m3,.7f);
		glutSetColor(i,m1, m2, m3);
	}
}
double fix_sine(double sine)
{
	if(sine >= 1)
		return 1;
	else if(sine <= -1)
		return -1;
	else
		return sine;
}

double compute_angle_rad(Another::Point P, Another::Point Q, Another::Point R)
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




double computeCot(Another::Point position_v_i, Another::Point position_v_j, Another::Point position_v_k, Another::Point position_v_l )
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
	CGAL_surface_mesh_parameterization_assertion(len != 0.0);    // two points are identical!
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
/* /brief  Caculate the discrete HarmonicMap                            */
/************************************************************************/
double* DiscreteHarmonicMap(Another::Another_HDS_model* model, int m, int n, int Length)
{
	//Allocate matrix A and B, A is a length*length matrix, the last colume is B
	double* matrixU = (double*)malloc(sizeof(double)*Length*(Length+1));

	for (int i = 0; i<Length; i++)
	{
		for (int j = 0; j<Length+1; j++)
		{
			matrixU[i*Length+j] = 0;
		}
	}
	//For each vertex, get a row of matrix A and B, according to
	//sigma ( j in neignbor i ){(ui - uj)(cotai +cotbi) }= 0
	//for i!= m and n and j!= m and n ( sigma(j in neighbor i) { (cotai+ cotbi) } ) ui - sigma(j in neighbor i) { (cotai + cot bi) uj }
	//
	for( Another::Vertex_iterator v = model->vertices_begin(); v != model->vertices_end();++v)
	{
		int tmpIndex  = v->GetIndex();
		Another::Point  Vi_position = v->point();
		if ( tmpIndex != m && tmpIndex != n ) //skip m and n rows 
		{
			Another::Halfedge_handle tmpEdge;
			Another::Halfedge_handle firstEdge = v->halfedge();
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

				int nIndex = vtmpj->GetIndex();

				double delta = 0.0;
				delta = computeCot(Vi_position, Vj_position, Vk_position, Vl_position);
				sum = sum + delta;
				if ( (nIndex != m) && (nIndex!=n) && (nIndex!=tmpIndex) ) //skip m and n columes  caculate uj
				{
					matrixU[ (tmpIndex -1 )*Length + nIndex -1 ] = -delta;

				}
				else if (nIndex == m) // for m um = 1 un = -1  
				{
					matrixU[(tmpIndex-1)*Length + Length] = matrixU[(tmpIndex-1)*Length + Length] + delta;
				}
				else if (nIndex == n)
				{
					matrixU[(tmpIndex-1)*Length + Length] = matrixU[(tmpIndex-1)*Length + Length] - delta;
				}
				else
				{

					cout<<"Error"<<endl;
				}
				//visite vl
				tmpEdge = Ejl->next();
			}

			matrixU[(tmpIndex-1)*Length + tmpIndex - 1 ] = sum;
		}
		else //Set m and n rows 0
		{
			for (int i = 0; i<Length+1; i++)
			{
				matrixU[(tmpIndex-1)*Length + i] = 0.0;
			}
		}
	}
	return matrixU;
}

double* SolveLinearEquation(double* matrixA, double* matrixB, int length)
{
	cout << "calling matlab" << endl;
	Engine *m_pEngine;
	m_pEngine = engOpen(NULL);
	if (m_pEngine == NULL)
	{
		//Error! Fail to connect to MATLAB engine.
		// The plot function will be disabled!
		printf("Fail to open MATLAB Engine!\n");
		exit(0);
	}	

	engSetVisible(m_pEngine, 0);
	mxArray *A = NULL;
	mxArray *B = NULL;

	A = mxCreateNumericMatrix( length -2 , length -2 , mxDOUBLE_CLASS, mxREAL) ;
	B = mxCreateNumericMatrix( length -2 , 1 , mxDOUBLE_CLASS , mxREAL) ;

	mxSetPr (A, matrixA) ;
	mxSetPr (B, matrixB) ;

	engPutVariable(m_pEngine, "A", A);
	engPutVariable(m_pEngine, "B", B);


	engEvalString(m_pEngine, "X = linsolve(A,B)"); 
	mxArray *X = NULL;
	X = engGetVariable(m_pEngine, "X");

	double* re = mxGetPr(X);


	FILE* refile;	
	refile = fopen("D:\\matrixX.txt","w+");
	for (int i = 0; i<length-2;i++)
	{
		fprintf(refile,"%f ",re[i]);
	}
	fclose(refile);

	return re;
}

double* GetA(double*matrix, int m, int n, int length )
{
	double* A = (double*)malloc(sizeof(double)*(length-2)*(length-2));
	int ai = 0;
	int aj = 0 ;
	//FILE* Afile;	
	//Afile = fopen("D:\\matrixA.txt","w+");

	for (int i = 0; i< length ; i++)
	{
		aj =0;
		if ( i!=(m-1) && i!=(n-1) )
		{
			for (int j = 0; j<length; j++)
			{
				if(j!= (m-1) &&j!= (n-1))
				{
					A[(length-2)*ai + aj ] = matrix[i*length+j];
					//fprintf(Afile,"%f ",A[(length-2)*ai + aj ]);
					aj++;
				}
				else //skip
					continue;
			}
			ai++;
			//	fprintf(Afile,"\n");
		}
		else //skip all the item
		{
			continue;
		}
	}
	//fclose(Afile);

	return A;
}
double* GetB(double*matrix, int m, int n, int length )
{
	double* B = (double*)malloc(sizeof(double)*(length-2));
	int bi = 0;
	//FILE* Bfile;	
	//Bfile = fopen("D:\\matrixB.txt","w+");
	for (int i = 0; i< length; i++)
	{
		if ( (i!= (m-1) ) && (i!= (n-1)) )
		{
			B[bi] = matrix[i*length+length];
			//fprintf(Bfile,"%f ",B[bi]);
			bi++;
		}
		else
		{
			continue;
		}
	}
	//fclose(Bfile);
	return B;
}



class mobiusWindow;
class mobiusApplication:public cwc::glApplication
{
private:
	vector<Another::Another_HDS_model*> m_models;
	mobiusWindow* m_view;
public:
	mobiusWindow* AttachWindow(mobiusWindow* view) 
	{
		m_view = view;
		return m_view;
	}

	/************************************************************************/
	/* /brief Insert a model in controller                                  */
	/************************************************************************/
	Another::Another_HDS_model* InsertModel(Another::Another_HDS_model* model) 
	{
		if (model)
		{
			m_models.push_back(model);
		}
		else
		{
			cerr<<"Error:: mobiusApplication::InsertModel Insert a NULL model"<<endl;
		}
		return model;
	}
	//Load model, Init render option
	virtual void OnInit()
	{
		// 		SetAppName("Mobius Voting");
		// 		ShowConsole();	

#pragma region Import Models

		//Load model

/************************************************************************/
/* M file load                                                          */
/************************************************************************/
	//	string sourceFilename("D:\\Research\\Correspondence\\3D correspondence\\Project\\MobiusVoting\\source_model.m");
		string sourceFilename("D:\\Models\\polycube data\\cat\\cat0_no_hole.obj2m.m");
		string targetFilename("D:\\Models\\polycube data\\cat\\cat_polycube_boundary_xyz.pp.m");

		FILE* f1 = fopen(sourceFilename.c_str(),"r");
		FILE* f2 = fopen(targetFilename.c_str(),"r");

		Another::Another_HDS_model* model1 = new Another::Another_HDS_model();
		model1->load_m_file(f1);
		g_model.push_back(model1);
		sourceModel = g_model.size()-1;
		fclose(f1);

		Another::Another_HDS_model* model2 = new Another::Another_HDS_model();
		model2->load_m_file(f2);
		g_model.push_back(model2);
		targetModel = g_model.size()-1;
		fclose(f2);


		/************************************************************************/
		/* OFF file load                                                        */
		/************************************************************************/

// 		string sourceFilename("D:\\Models\\OFFfiles\\cat0_no_hole.off");
// 		string targetFilename("D:\\Models\\OFFfiles\\cat1_no_hole.off");
//  		sourceModel = ImportModel(sourceFilename);
//  		targetModel = ImportModel(targetFilename);

#pragma endregion

#pragma region Embedding

		LoadLaplace(sourceModel);
		LoadLaplace(targetModel);

 		LaplacianMatrix(sourceModel);
 		LaplacianMatrix(targetModel);

#pragma endregion

#pragma region Harmonic map

		HarmonicMap(sourceModel);
		HarmonicMap(targetModel);

#pragma endregion

#pragma region Calculate gussain curvature sampling

		GuassianCurvature(sourceModel);
		GuassianCurvature(targetModel);

#pragma endregion

#pragma region Voting


		
#pragma endregion

	}

	void LoadLaplace(int model)
	{
		cout<<"Loading Laplacian coordinates for model: "<< model <<endl;

		//Get the model
		Another::Another_HDS_model* tmpModel = g_model[model];
		ifstream laplaceFile;
		string name = "LaplaceCoordinates";
		string pofix = ".txt";
		string s = boost::lexical_cast<string>( model );
		name += s;
		name += pofix;
		laplaceFile.open(name.c_str());
		tmpModel->load_laplacian_mesh(laplaceFile);
		laplaceFile.close();

		string mfileposix = ".m";
		string savefilename = "model";
		savefilename +=s;
		savefilename +=mfileposix;
		tmpModel->save_m_file(savefilename,1);
		cout<<"Finish Loading Laplacian matrix " <<endl;
	}


	void LaplacianMatrix(int model)
	{
		cout<<"Calculate Laplacian matrix for model: "<< model <<endl;

		//Get the model
		Another::Another_HDS_model* tmpModel = g_model[model];
		if (tmpModel)
		{
			Another::GeometryProcessingAlgorithm al;
			double* lpMatrix = al.LaplacianMatrix(*tmpModel,Another::LapalacianOperatorType::COTANGENT);
			
			if (lpMatrix == NULL)
			{
				return;
			}

			//Save the matrix to file
			int n = tmpModel->size_of_vertices();
			ofstream matrixFile;
			string name = "LaplacaianMatrix";
			string pofix = ".txt";
			string s = boost::lexical_cast<string>( model );
			name += s;
			name += pofix;
			matrixFile.open(name.c_str());
			if (matrixFile.good())
			{
				for (int i = 0; i<n; ++i)
				{
					for (int j = 0; j<n; ++j)
					{
						matrixFile<<lpMatrix[i*n+j]<<" ";
					}
					matrixFile<<"\n";
				}
			}
			matrixFile.close();
		}
		cout<<"Finish Calculating Laplacian matrix " <<endl;
	}

	/************************************************************************/
	/* /brief Load model from m file                                        */
	/* /param [in] filename													*/
	/* /return the index of the loaded model in g_model						*/
	/************************************************************************/
	Another::Another_HDS_model* ImportModel_OFF(string filename)
	{
		cout<<"Importing model "<<filename<<endl;
		if ( filename.c_str() != NULL )
		{
			Another::Another_HDS_model* tmpModel = new Another::Another_HDS_model();
			ifstream stream;
			stream.open(filename.c_str());
			if (stream.good())
			{
				stream >> *tmpModel;
				stream.close();
				tmpModel->setup_OneToOneVertexMap();
				cout<<"Imported model:"<<filename<<endl;
				InsertModel(tmpModel);
				return tmpModel;
			}
			else
			{
				cerr << "open file:" << filename << "Failed"<<endl;
				stream.close();
				return NULL;
			}

		}
		else
		{
			cerr<<"file name is null"<<endl;
			return NULL;
		}
	}

	/************************************************************************/
	/* /brief Load model from file                                          */
	/* /param [in] filename													*/
	/* /return the index of the loaded model in g_model						*/
	/************************************************************************/
	int ImportModel(string filename)
	{
		cout<<"Importing model "<<filename<<endl;
		Another::Another_HDS_model* tmpModel;
		if ( filename.c_str() != NULL )
		{
			tmpModel = new Another::Another_HDS_model();
			ifstream stream;
			stream.open(filename.c_str());
			if (stream.good())
			{
				stream >> *tmpModel;
				stream.close();
				g_model.push_back(tmpModel);
				tmpModel->setup_OneToOneVertexMap();
				cout<<"Finish importing model"<<endl;
				return g_model.size()-1;
			}
			else
			{
				cerr << "open file:" << filename << "Failed"<<endl;
				stream.close();
				return -1;
			}

		}
		else
		{
			cerr<<"file name is null"<<endl;
			return -1;
		}
	}

	/************************************************************************/
	/* /brief Calculate the harmonic map and mid-edge embedding of a model  */
	/************************************************************************/
	void HarmonicMap(int model)
	{
		cout<<"Harmonic mapping for mesh, the boundary is cut from a triangle"<<endl;
		cout<<"Calculate harmonic mapping for model<<<" << model <<endl;

		//Get the model
		Another::Another_HDS_model* tmpModel = g_model[model];
		
		if (tmpModel)
		{
			//Specify the cut face
			int m , n;
			cout<<"cut edge:";
			cin >> m;
			cin >> n;
			//Harmonic map
			int length = tmpModel->size_of_vertices();
			double* matrix = DiscreteHarmonicMap(tmpModel, m, n, length);
			double* matrixA = GetA(matrix, m, n , length);
			double* matrixB = GetB(matrix, m , n, length);
			free(matrix);
			matrix = NULL;
			double* result = SolveLinearEquation(matrixA, matrixB,length);
			tmpModel->set_U(result, m,n);
			//free(result);
			result = NULL;
			free(matrixA);
			matrixA = NULL;
			free(matrixB);
			matrixB = NULL;
			free(result);
			cout<<"Finish calculating harmonic map"<<endl;

			//Planar embedding
			cout<<"Start planar embedding"<<endl;
			tmpModel->PlanarEmbedding();
			cout<<"Finish planar embedding"<<endl;
		}
		else
		{
			cerr <<"Model:"<<model<<" dose not exist, can not calculate the harmonic map"<<endl;
			return;
		}
	}

	/************************************************************************/
	/* /brief Calculate the Gaussian curvature of model                      */
	/************************************************************************/
	void GuassianCurvature(int model)
	{
		Another::Another_HDS_model* tmpModel = g_model[model];
		if (tmpModel)
		{
			//Options
			unsigned int d_fitting = 2 ;
			unsigned int d_monge = 2;
			unsigned int nb_rings = 0;//seek min # of rings to get the required #pts
			unsigned int nb_points_to_use = 0;//
			bool verbose = false;
			unsigned int min_nb_points = (d_fitting+1)*(d_fitting+2)/2;
			unsigned int neihborPointsToUse = 0;

			//modify global variables which are fact of options:
			min_nb_points = (d_fitting + 1) * (d_fitting + 2) / 2;
			if (nb_points_to_use < min_nb_points && nb_points_to_use != 0)
			{std::cerr << "the nb of points asked is not enough to perform the fitting" << std::endl; exit(0);}

			//Calculate Gaussian curvature

			cout<<"Calculating Gaussian curvature for model:"<<model<<endl;
			Another::GeometryProcessingAlgorithm al;
			Another::Vertex2Monge_map_type vertex2Mongemap;
			Another::Vertex2Monge_PM_type mongePM(vertex2Mongemap);
			al.Monge_jet_fitting(*tmpModel,min_nb_points,d_fitting,d_monge,nb_rings,neihborPointsToUse, mongePM);
			Another::Vertex2GuassianCurvature_map_type vertex2Guassainmp;
			Another::Vertex2GuassianCurvature_PM_type gaussianCurvature(vertex2Guassainmp);
			al.GaussianCurvature(mongePM,*tmpModel,gaussianCurvature);
			Another_CommonUtil::CommonUtil* util = new Another_CommonUtil::CommonUtil();
			vector<double> gc;
			tmpModel->get_guassian_curvature(gc);
			vector<double> gcn;
			util->NormalizeScalar(gc, gcn);
			if(tmpModel->set_1D_texture_coordinates(gcn))
			{
				cout<<"Finish calculating Gaussian curvature..."<<endl;
			}
			else
			{
				cerr<<"assign 1D texture coordinates failed"<<endl;
			}

			cout<<"Calculate the local minum and maxium of guassian curvature..."<<endl;
			std::vector<Another::Vertex_handle> features;
			al.LocalMaxandMinGuassian(*tmpModel,5, features);
			cout<<"There are " <<features.size()<<" feature points of model"<<model<<endl;
			tmpModel->set_featurePoint(features);
			cout<<"Finish calculating the local minum and maxium of guassian curvature..."<<endl;
			delete util;

			cout<<"Start sampling"<<endl;

			int SamplePointNumber = 100;
			cout<<"Enter sample point number:";
			cin >> SamplePointNumber;
			sample_size = SamplePointNumber;
			Another::Sampling_algorithm samplingAl;
			samplingAl.FPS(features,*tmpModel,SamplePointNumber);
			tmpModel->set_featurePoint(features);

			cout<<"End sampling"<<endl;
		}
		else
		{
			cerr <<"Model:"<<model<<" is not exist, can not calculat the harmonic map"<<endl;
			return;
		}
	}
	
};

enum DualView
{
	LEFT = 0,
	RIGHT
};
class mobiusWindow:public cwc::glutWindow
{
private:
	mobiusApplication* m_application;
	Another::Another_HDS_model* m_model;
protected:
	cwc::glShaderManager m_shaderManager;
	cwc::glShader *m_shader;
	GLuint m_programObject;
	map<Another::RenderOption , cwc::SmartPtr<cwc::TextureBase> > m_textures; //textures
	cwc::SmartPtr<cwc::TextureBase> m_texture;
	Another::RenderOption m_renderOption; //render option
	bool m_show_texture;

	//Feature option
	bool m_set_correspondence;
	

	//status
	bool m_correspoendence; //has correspondence or not
	bool m_planar; //has planar embedding or not

	DualView m_active_window;
	bool m_dual_window;


public:
	mobiusWindow() { m_renderOption = Another::RenderOption::DEFAULT;  /*  OnAddMenu();*/ m_correspoendence = false; m_planar = false;}
	mobiusWindow(mobiusApplication* app):m_renderOption(Another::RenderOption::DEFAULT),
		m_correspoendence(false),
		m_planar(false),
		m_application(app),
		m_texture(NULL),
		m_show_texture(false),
		m_set_correspondence(false),
		m_dual_window(false)
	{
		// Create textures
		// Todo: change the m_textures vector type to map type
		m_textures.clear();
		cwc::SmartPtr<cwc::TextureBase> hue_texture = cwc::TextureFactory::CreateTextureFromFile(".\\textures\\hue.png");
		pair<Another::RenderOption,cwc::SmartPtr<cwc::TextureBase>> curvature_texture(Another::RenderOption::CURVATURE,hue_texture);
		m_textures.insert(curvature_texture);
		cwc::SmartPtr<cwc::TextureBase> white_black_check_box = cwc::TextureFactory::CreateCheckBoxTexture(cwc::CheckBoxTextureType::WHITE_BLACK_CHECK_BOX, 64,64);
		pair<Another::RenderOption,cwc::SmartPtr<cwc::TextureBase>> harmonic_texture(Another::RenderOption::DISCRETE_HARMONIC_MAP,white_black_check_box);
		m_textures.insert(harmonic_texture);
		if (m_textures.size() ==0)
			std::cerr << "***WARNING: Failed loading texture!!"<<endl;
	}
	Another::Another_HDS_model* AttachModel(Another::Another_HDS_model* model) 
	{
		m_model = model;
		return m_model;
	}
	virtual void OnRender(void)
	{
		/************************************************************************/
		/* Clear color, Set up texture, shader, light and viewpoints             */
		/************************************************************************/
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

		//Preparing texture and bind def if has a texture
		if (m_texture)
		{
			m_texture->bind(0);
			glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);	
			glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST); // Disable Filtering!
		}

		//Get the window size
		int l_wsW = GetWidth();
		int l_wsH= GetHeight();

		//Prepare variables for the shader
		GLfloat baseColor[3] = {0.4, 0.4, 1.0};
		GLfloat lightPos[3] = {0.0, 0.0, 4.0};

		switch(m_renderOption)
		{
		case Another::RenderOption::CORRESPONDENCE:

			//Left view for source
			glViewport(0, 0, l_wsW/2, l_wsH);

			//Preparing shader;
			if(m_shader) 
			{
				glEnable(GL_TEXTURE_2D);
				m_shader->setUniform1i("EnvMap", 0);
				m_shader->setUniform3fv("BaseColor", 3, baseColor);
				m_shader->setUniform3fv("LightPos", 3, lightPos);

				m_shader->begin();
				g_model[0]->draw(SMOOTH,1,1);
				m_shader->end();
				glDisable(GL_TEXTURE_2D);
			}
			g_model[sourceModel]->draw_feature_points(GL_SMOOTH);
			//	set_curvature();

			//Right view for target
			glViewport(l_wsW/2, 0, l_wsW/2, l_wsH);
			if (m_shader) 
			{
				glEnable(GL_TEXTURE_2D);
				m_shader->setUniform1i("EnvMap", 0);
				m_shader->setUniform3fv("BaseColor", 3, baseColor);
				m_shader->setUniform3fv("LightPos", 3, lightPos);
				m_shader->begin();
				g_model[targetModel]->draw(SMOOTH,1,1);
				m_shader->end();
				glDisable(GL_TEXTURE_2D);
			}
			g_model[targetModel]->draw_feature_points(GL_SMOOTH);
			break;

		case Another::RenderOption::FEATURE_POINT:

			//left view for source
			glViewport(0, 0, l_wsW/2, l_wsH);
			g_model[sourceModel]->draw(m_renderOption);

			//right view for target
			glViewport(l_wsW/2, 0, l_wsW/2, l_wsH);
			g_model[targetModel]->draw(m_renderOption);
			break;
		case Another::RenderOption::PLANAR_EMBEDDING: 
			//left view for source
			glViewport(0, 0, l_wsW/2, l_wsH);
			glDisable(GL_TEXTURE_2D);
			g_model[sourceModel]->draw(m_renderOption);

			//right view for target
			glViewport(l_wsW/2, 0, l_wsW/2, l_wsH);
			g_model[targetModel]->draw(m_renderOption);
			break;
		case Another::RenderOption::CURVATURE: 
		case Another::RenderOption::TRIANGULAR: 
		case Another::RenderOption::DISCRETE_HARMONIC_MAP:
			glViewport(0, 0, l_wsW, l_wsH);
			glEnable(GL_TEXTURE_2D);
			g_model[sourceModel]->draw(m_renderOption);
			glDisable(GL_TEXTURE_2D);
			break;
		case Another::RenderOption::LAPLACE:
			//left view for source
			glViewport(0, 0, l_wsW/2, l_wsH);
			g_model[sourceModel]->draw(m_renderOption);

			//right view for target
			glViewport(l_wsW/2, 0, l_wsW/2, l_wsH);
			g_model[targetModel]->draw(m_renderOption);
			break;
		case Another::RenderOption::DEFAULT:
			glViewport(0, 0, l_wsW, l_wsH);
			g_model[sourceModel]->draw(m_renderOption);
		}

	
// 		glMatrixMode(GL_PROJECTION);
// 		glLoadIdentity();
// 		glOrtho(-l_wsW/4, l_wsW/4, l_wsH/2, -l_wsH/2, 4.0, 10000.0);
// 		glMatrixMode(GL_MODELVIEW);

		if (m_show_texture)
		{
			glEnable(GL_TEXTURE_2D);
			m_texture->drawToScreen(0,0,10,10);
			glDisable(GL_TEXTURE_2D);
		}
		glutSwapBuffers();
	}

	virtual void OnIdle() {}

	// When OnInit is called, a render context (in this case GLUT-Window) 
	// is already available!
	virtual void OnInit()
	{
		//OpenGL Init
		glutInitDisplayMode(GLUT_SINGLE | GLUT_INDEX);
		glClearColor(0.5f, 0.5f, 0.5f, 0.0f);
		glShadeModel(GL_SMOOTH);
		glEnable(GL_DEPTH_TEST);
		
		//Load default glass shader
		m_shader = m_shaderManager.loadfromFile(".\\shaders\\Glass.vert",".\\shaders\\Glass.frag"); // load (and compile, link) from file
		if (m_shader==0) 
			std::cerr<< "Error Loading, compiling or linking shader"<<endl;
		else
		{
			m_programObject = m_shader->GetProgramObject();
		}

	}

	virtual void OnResize(int w, int h) 
	{
		// 		if(h == 0) h = 1;
		// 		float ratio = 1.0f * (float)w / (float)h;
		// 
		// 		glMatrixMode(GL_PROJECTION);
		// 		glLoadIdentity();
		// 
		// 		glViewport(0, 0, w, h);
		// 
		// 		 gluPerspective(45,ratio,1,100);
		// 		//glOrtho(-810,810,-540,540,1,10000);
		// 		glMatrixMode(GL_MODELVIEW);
		// 		glLoadIdentity();
		// 		gluLookAt(0.0f,0.0f,15.0f, 
		// 			0.0,0.0,-1.0,
		// 			0.0f,1.0f,0.0f);

			zprReshape(w,h);

	}
	virtual void OnClose(void){}
	virtual void OnMouseDown(int button, int x, int y) 
	{
		//If the rendering statue in setting correspondence, redraw the active window
		//If the rendering statue in setting correspondence, track the active window
		if (m_dual_window)
		{
			//Get the window size
			int l_wsW = GetWidth();

			//set the active window
			if ( x< l_wsW/2 +1 )
			{
				m_active_window = DualView::LEFT;
			}
			else
			{
				m_active_window = DualView::RIGHT;
			}
		}
		if ( m_set_correspondence&&m_dual_window )
		{
			int l_wsW = GetWidth();
			int l_wsH = GetHeight();
			if (m_active_window == DualView::LEFT)
			{
				//left view for source
				glViewport(0, 0, l_wsW/2, l_wsH);
				g_model[sourceModel]->draw(m_renderOption);
				zprSelectionFunc(drawPickSource);     /* Selection mode draw function */
				zprPickFunc(ProcessPickUpSource);              /* Pick event client callback   */
			}
			else
			{
				//right view for target
				glViewport(l_wsW/2, 0, l_wsW/2, l_wsH);
				g_model[targetModel]->draw(m_renderOption);
				zprSelectionFunc(drawPickTarget);     /* Selection mode draw function */
				zprPickFunc(ProcessPickUpTarget);              /* Pick event client callback   */
			}
		}
				
		zprMouse(button, GLUT_DOWN, x, y);

	}    
	virtual void OnMouseUp(int button, int x, int y)
	{
		zprMouse(button, GLUT_UP, x, y);
	}
	virtual void OnLeftMouseDrag(int x, int y)
	{
		zprMotion(x, y);
	}
	virtual void OnMouseMove(int x, int y)
	{
		zprMotion(x,y);
	}
	virtual void OnMouseWheel(int nWheelNumber, int nDirection, int x, int y){}

	virtual void OnKeyDown(int nKey, char cAscii)
	{       
		//m_application->OnKeyDown(nKey, cAscii);
		switch(cAscii)
		{
		case 27:	// 0x1b = ESC
			Close();// Close Window!
			break;
		case 'a':
			m_renderOption = Another::RenderOption::LAPLACE;
			break;
		case 'm': //show original models todo
			m_renderOption = Another::RenderOption::MESH;
			break;
		case 'c': //show correspondence todo
			m_renderOption = Another::RenderOption::CURVATURE;
			break;
		case 'p'://planar embedding option
			if (!m_planar)
			{
				//project the feature point on nearst vertex of mid-edge mesh
				Embed(sourceModel,targetModel);
				m_planar = true;
			}
			m_renderOption = Another::RenderOption::PLANAR_EMBEDDING;
			break;
		case 'h':
			m_renderOption = Another::RenderOption::DISCRETE_HARMONIC_MAP;
			m_show_texture = false;
			m_texture = m_textures[m_renderOption];
			break;
		case 't':
			m_renderOption = Another::RenderOption::TRIANGULAR;
			break;
		case 'f':
			m_renderOption = Another::RenderOption::FEATURE_POINT;
			m_set_correspondence = false;
			m_dual_window = true;
			break;
		case 'v':
			if (!m_correspoendence)
			{
				vector<pair<int,int>> pairs;
				Voting(sourceModel,targetModel,pairs);
				SetCorrespoendence(sourceModel,targetModel,pairs);
				Another::Another_HDS_model::init_color_map();
				m_correspoendence = true;
				m_renderOption = Another::RenderOption::CORRESPONDENCE;
			}
			else
			{
				cout<<"Already construct correspondence"<<endl;
			}
			break;
		case 's':
			{
				switch(m_renderOption)
				{
				//Discrete harmonic map mode, 's' key to toggle 2D texture display
				case Another::RenderOption::DISCRETE_HARMONIC_MAP:
					m_show_texture = !m_show_texture;
					break;

				//Feature point mode, 's' key to toggle correspondence selection
				case Another::RenderOption::FEATURE_POINT:
					m_set_correspondence = !m_set_correspondence;
					if (m_set_correspondence)
					{
						sourceFeatures.clear();
						targetFeatures.clear();
						cout<<"<<<<<<<<<<<<<<"<<endl;
						cout<<"<<<<<<<<<<<<<<"<<endl;
						cout<<"<<<<<<<<<<<<<<"<<endl;
						cout<<"Please Select Feature Point..."<<endl;
					}
					else
					{
						if (sourceFeatures.size()!=targetFeatures.size())
						{
							cout<<"Number of feature points on source model is not equal to feature points on target model"<<endl;
						}
						else
						{
							ofstream sourceFeatureStream;
							sourceFeatureStream.open("source_vetices_list.txt");
							for (int i = 0; i<sourceFeatures.size();i++)
							{
								sourceFeatureStream<<sourceFeatures[i]<<" ";
							}
							sourceFeatureStream.close();

							ofstream targetFeatureStream;
							targetFeatureStream.open("target_vertices_list.txt");
							for (int i = 0; i<targetFeatures.size();i++)
							{
								targetFeatureStream<<targetFeatures[i]<<" ";
							}
							targetFeatureStream.close();

							cout<<"Output: source_vetices_list.txt and target_vertices_list.txt"<<endl;

							g_model[sourceModel]->save_m_file("source_model.m");
							g_model[targetModel]->save_m_file("target_model.m");

							cout<<"Output: source_model.m and target_model.m"<<endl;
						}
					}
					break;
				}
			}
		break;

		case 'd':
			Delaunay_triangulation( g_model[sourceModel]);
		default:
			break;
		} 
		
		Repaint();
	}

	virtual void OnKeyUp(int nKey, char cAscii)
	{
		if (cAscii == 's')      // s: Shader
			m_shader->enable();
// 		else if (cAscii == 'f') // f: Fixed Function
// 			m_shader->disable();
	}

	virtual void OnMenu(int itemi)
	{
		switch (itemi) 
		{
		case ITEM_1:
			open_file();

		}

	}


	virtual void OnAddMenu()
	{
		glutAddMenuEntry("Open ",ITEM_1);
		glutAttachMenu(GLUT_RIGHT_BUTTON);
	}

	void Delaunay_triangulation(Another::Another_HDS_model* model)
	{
		Another::Triangulation_algorithm al;
		Another::Another_triangulation_data_structure tds;
		vector<int> features;
		model->get_sample_indexes(features);
		al.triangulate(*model, features, tds);
		//al.read_triangulation_file("result.txt");
		vector<Another::Another_HDS_model> patches;
		al.segment_mesh(*model,tds, patches);
	}

	void open_file()
	{
		static char szFilter[]= "Ply files (*.ply)\0*.ply\0M files (*.m)\0*.m\0OFF Files (*.off)\0*.off\0All Files (*.obj)\0*.obj\0All Files (*.*)\0*.*\0";
		OPENFILENAME ofn;
		char pszFileLocn[256] = {'\0'};

		// Set up OPENFILENAME struct to use commond dialog box for open

		ZeroMemory(&ofn, sizeof(OPENFILENAME));
		ofn.lStructSize = sizeof(ofn);	// size of struct
		ofn.hwndOwner = NULL;			// window that owns Dlg
		ofn.lpstrFilter = szFilter;		// Filter text
		ofn.lpstrFile = pszFileLocn;		// File name string
		ofn.nMaxFile = sizeof(pszFileLocn); // size of file name
		ofn.Flags = OFN_HIDEREADONLY;	// don't display "Read Only"
		ofn.lpstrDefExt = "ply";		// extension name
		ofn.lpstrTitle = "Open Mesh File"; // title of dlg box

		// call common dlg control for file open
		if (!GetOpenFileName(&ofn)) 
		{
			return ;	
		}

		// see if file exists
		// 		struct stat fileStat;
		// 		if (stat(ofn.lpstrFile, &fileStat))
		// 		{
		// 			char errormsg[1024];
		// 			sprintf(errormsg, "%s not found.", ofn.lpstrFile);
		// 			MessageBox(NULL,errormsg,"File Not Found Error",MB_OK | MB_ICONINFORMATION);
		// 			return ;	
		// 		}


		Another::Another_HDS_model* model = new Another::Another_HDS_model();
		switch(ofn.nFilterIndex)
		{
		case 1:
			//g_model->load_ply_file(mfile);
			break;
		case 2:

			load_m_file(ofn.lpstrFile, model);
			g_model.push_back(model);
			g_current_model = g_model.size()-1;
			break;
		case 3:

			load_off_file(ofn.lpstrFile, model);
			g_model.push_back(model);
			g_current_model = g_model.size()-1;
			break;
			// 		case 4:
			// 			load_obj_file(ofn.lpstrFile, model);
			// 			g_model.push_back(model);
			// 			g_current_model = g_model.size()-1;
		default:
			break;
		}

	}
	void load_m_file(const char* fileName, Another::Another_HDS_model* model)
	{
		FILE *mfile = fopen(fileName,"r");
		if (!mfile)
		{
			//MessageBox("File can not open",MB_OK | MB_ICONINFORMATION);
			return;
		}
		model->load_m_file(mfile);
		fclose(mfile);
	}

	void load_off_file(const char* fileName, Another::Another_HDS_model* model)
	{
		ifstream stream;
		stream.open(fileName);
		stream >> *model;
		stream.close();
	}

	// 	void load_obj_file(const char* fileName, Another::Another_HDS_model* model)
	// 	{
	// 		FILE *mfile = fopen(fileName,"r");
	// 		if (!mfile)
	// 		{
	// 			//MessageBox("File can not open",MB_OK | MB_ICONINFORMATION);
	// 			return;
	// 		}
	// 		model->load_obj_file(mfile);
	// 		fclose(mfile);
	// 	}

	void Voting(int source, int target, vector<pair<int,int>>&pairs)
	{
		cout<<"Start voting.."<<endl;
		Another::Another_HDS_model* src = g_model[source];
		Another::Another_HDS_model* targ = g_model[target];
		Another::GeometryProcessingAlgorithm geoAl;
		double* correspoendenceMatrix = geoAl.MobiusVoting(*src, *targ);
		pairs.clear();
		ProcessingMatrix(correspoendenceMatrix, pairs);
		cout<<"End voting, there are "<<pairs.size()<<" points are matched"<<endl;
	}

	/************************************************************************/
	/* /brief embed source and target model in planar domain                */
	/* /param source[in] souce model index									*/
	/* /param target[in] target model index									*/
	/************************************************************************/
	void Embed(int source, int target)
	{
		cout<<"Start Embed.."<<endl;
		Another::Another_HDS_model* src = g_model[source];
		Another::Another_HDS_model* targ = g_model[target];
		Another::GeometryProcessingAlgorithm geoAl;
		geoAl.Embed(*src, *targ);
		cout<<"End Embed.."<<endl;
	}

	/************************************************************************/
	/* /brief                                                               */
	/************************************************************************/
	void ProcessingMatrix(double* correspoendenceMatrix, vector<pair<int,int>>& pairs)
	{
		pair<int, int> tmpPair;
		double tmpEntry = correspoendenceMatrix[0];
		double largest = 0.0;
		do
		{
			largest = 0.0;
			for ( int k = 0; k<sample_size; k++ )
			{
				for ( int l = 0; l<sample_size; l++ )
				{
					tmpEntry = correspoendenceMatrix[k*sample_size+l];
					if (tmpEntry > largest)
					{
						tmpPair.first = k;
						tmpPair.second = l;
						largest = tmpEntry;
					}
				}

			}
			if (largest >0) //if find a entry is larger than zero and is the largest
			{
				pairs.push_back(tmpPair);

				for (int i = 0; i<sample_size; i++)
				{
					correspoendenceMatrix[tmpPair.first*sample_size+i] = 0.0;
					correspoendenceMatrix[i*sample_size+tmpPair.second] = 0.0;
				}
			}

		}while(largest>0);

	}
	void SetCorrespoendence(int sourceModel, int targetModel,vector<pair<int,int>>& pairs)
	{
		Another::Another_HDS_model* src = g_model[sourceModel];
		Another::Another_HDS_model* targ = g_model[targetModel];
		vector<int> IndexSrc;
		IndexSrc.clear();
		vector<int> IndexTag;
		IndexTag.clear();
		for (int i = 0; i<pairs.size(); i++)
		{
			IndexSrc.push_back(pairs[i].first);
			IndexTag.push_back(pairs[i].second);
		}
		src->set_correspondence(IndexSrc);
		targ->set_correspondence(IndexTag);
	}
};


bool Setup2DdrawingMode(double windowLeft, double windowRight, double windowBottom, double windowTop)
{
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	glOrtho(windowLeft, windowRight, windowBottom, windowTop, 0, 1);
	glMatrixMode (GL_MODELVIEW);
	return true;
}
// 
// void renderScene(void) 
// {
// 
// 	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
// 
// 	glLoadIdentity();
// 	gluLookAt(0.0,0.0,10.0, 
// 		0.0,0.0,-1.0,
// 		0.0f,1.0f,0.0f);
// 
// 	glLightfv(GL_LIGHT0, GL_POSITION, lpos);
// 	if ( g_model == NULL )
// 	{
// 		glutSolidTeapot(1);
// 	}
// 	else
// 	{
// 		switch ( g_option )
// 		{
// 		case RenderOption::DISCRETE_HARMONIC_MAP:
// 			g_model->draw_HamonicMap();
// 		case RenderOption::PLANAR_EMBEDDING:
// 			g_model->draw_midedge( );
// 		case RenderOption::MESH:
// 			g_model->draw(GL_SMOOTH,1,1);
// 			// 		case RenderOption::CURVATURE:
// 			// 			g_model->draw(CURVATURE);
// 		}
// 	}
// 	//a+=0.1;
// 	glutSwapBuffers();
// }

void processNormalKeys(unsigned char key, int x, int y) {

	if (key == 27) 
		exit(0);
}

#define printOpenGLError() printOglError(__FILE__, __LINE__)

int printOglError(char *file, int line)
{
	//
	// Returns 1 if an OpenGL error occurred, 0 otherwise.
	//
	GLenum glErr;
	int    retCode = 0;

	glErr = glGetError();
	while (glErr != GL_NO_ERROR)
	{
		printf("glError in file %s @ line %d: %s\n", file, line, gluErrorString(glErr));
		retCode = 1;
		glErr = glGetError();
	}
	return retCode;
}


void printShaderInfoLog(GLuint obj)
{
	int infologLength = 0;
	int charsWritten  = 0;
	char *infoLog;

	glGetShaderiv(obj, GL_INFO_LOG_LENGTH,&infologLength);

	if (infologLength > 0)
	{
		infoLog = (char *)malloc(infologLength);
		glGetShaderInfoLog(obj, infologLength, &charsWritten, infoLog);
		printf("%s\n",infoLog);
		free(infoLog);
	}
}

void printProgramInfoLog(GLuint obj)
{
	int infologLength = 0;
	int charsWritten  = 0;
	char *infoLog;

	glGetProgramiv(obj, GL_INFO_LOG_LENGTH,&infologLength);

	if (infologLength > 0)
	{
		infoLog = (char *)malloc(infologLength);
		glGetProgramInfoLog(obj, infologLength, &charsWritten, infoLog);
		printf("%s\n",infoLog);
		free(infoLog);
	}
}



void setShaders() 
{

	char *vs = NULL,*fs = NULL,*fs2 = NULL;

	v = glCreateShader(GL_VERTEX_SHADER);
	f = glCreateShader(GL_FRAGMENT_SHADER);
	f2 = glCreateShader(GL_FRAGMENT_SHADER);

	vs = textFileRead(".\\shaders\\minimal.vert");
	fs = textFileRead(".\\shaders\\minimal.frag");

	const char * vv = vs;
	const char * ff = fs;

	glShaderSource(v, 1, &vv,NULL);
	glShaderSource(f, 1, &ff,NULL);

	free(vs);free(fs);

	glCompileShader(v);
	glCompileShader(f);

	printShaderInfoLog(v);
	printShaderInfoLog(f);
	printShaderInfoLog(f2);

	p = glCreateProgram();
	glAttachShader(p,v);
	glAttachShader(p,f);

	glLinkProgram(p);
	printProgramInfoLog(p);

	glUseProgram(p);

}



bool DeleteModel(Another::Another_HDS_model* model)
{
	try
	{
		delete model;
	}
	catch (exception&)
	{
		return false;
	}
	return true;
}


int main(int argc, char **argv)
{
	//Init freeImage
#ifdef FREEIMAGE_LIB
	FreeImage_Initialise();
#endif

	//Create a default model
	Another::Another_HDS_model* model = new Another::Another_HDS_model();

	//Create Controller
	mobiusApplication* app = new mobiusApplication();

	//Set up the MVC links: one-to-one multrual link between view and controller, insert models into controllers and views
	mobiusWindow* win = new mobiusWindow(app);
	app->AttachWindow(win);
	app->InsertModel(model);

	app->run();
	delete app;

#ifdef FREEIMAGE_LIB 
	FreeImage_DeInitialise();
#endif 
	return 0;
}