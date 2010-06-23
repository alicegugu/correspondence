//STD
#include <stdio.h>
#include <stdlib.h>
#include <exception>
#include <vector>
#include <iostream>

//GLEW
#include <GL/glew.h>

//GLUT
#include <GL/glut.h>

//BOOST
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>

//CGAL
#include <CGAL/Parameterizer_traits_3.h>
#include <CGAL/IO/Polyhedron_iostream_ture.h> // this head file is edited by Gu Huiqin and is not a standard CGAL file
#include <CGAL/Taucs_solver_traits.h>

//MATLAB
#include <engine.h>

//CWC
#include "glApplication.h"
#include "glutwindow.h"
#include "texture.h"
#include "glsl.h"
#include "smartptr.h"

//LOCAL
#include "Textfile.h"
#include "Another_HDS_model.h"
#include "Geometry_processing_algorithm.h"
#include "CommonUtil.h"
#include "zpr.h"
#include "Another_CGAL_cast.h"
#include "geodesic_algorithm_dijkstra.h"
#include "Sampling_algorithm.h"
#include "Triangulation_algorithm.h"

//FreeImage
#include "FreeImage.h"

//STL
#include <algorithm>

//Windows
#include <windows.h>
#include <Commdlg.h>
