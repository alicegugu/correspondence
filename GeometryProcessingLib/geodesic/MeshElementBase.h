#pragma once

namespace Another
{


	enum PointType
	{
		VERTEX_TYPE = 0,
		EDGE_TYPE = 1,
		FACE_TYPE = 2,
		UNDEFINED_POINT_TYPE = 3
	};

	class Mesh_element_base	//prototype of vertices, edges and faces
	{
	public:

		Mesh_element_base():m_type(UNDEFINED_POINT_TYPE){};
		PointType type(){return m_type;};
		PointType m_type;							//vertex, edge or face
	};
	typedef  Mesh_element_base* Base_handle;

}