#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Color.h>
//#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_data_structure_2.h>

#include "Another_HDS_model.h"


namespace Another
{
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef CGAL::Triangulation_vertex_base_with_info_2<int,K> Vb;
	typedef CGAL::Triangulation_face_base_with_info_2<Another::Face_handle,K> Fb;
	typedef CGAL::Triangulation_data_structure_2<Vb,Fb> TDS;

	class Another_triangulation_data_structure:public TDS
	{
	public:
		Another_triangulation_data_structure(void);
		~Another_triangulation_data_structure(void);
		//friend std::ifstream& operator <<(std::istream& in, Another_triangulation_data_structure & tds);
	};

}