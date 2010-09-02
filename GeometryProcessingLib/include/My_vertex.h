#pragma once

#include "MeshElementBase.h"
 
namespace Another
{

	/****************************************************************************************************************/
	/*************An Vertex type with following additional information:****************************************/
	/**************1. normal (The normal of vertex, assigned when is read from m file)**********************/
	/**************2. polyPoint (Polycube vertex coordinates, which is assigned when reading m file)*****/ 
	/**************3. cornerID (The corner ID of the vertex, if it is not a corner, the value is UNDEFINED)*/
	/**************4. planePoint: the coordinates on the 2D plane facet **************************************/
	/**************5. index (Index of the vertex)******************************************************************/
	/**************6. vertex_type (CORNER EDGE or FACE, firstly initialized with UNDEFINED)******************/
	/**************7. edgeCount (The edge connected to the vertex, first initialized with UNDEFINED)****/
	/****************************************************************************************************************/

	template<class Refs, class Tag, class Pt, class GeoTraits >
	struct My_vertex:public CGAL::HalfedgeDS_vertex_base<Refs,Tag ,Pt>, public Mesh_element_base
	{
		typedef typename GeoTraits::Point_3 Point_3;
		typedef typename GeoTraits::Point_2 PlanePoint;
		typedef typename GeoTraits::Vector_3 Normal;
		typedef CGAL::HalfedgeDS_vertex_base<Refs,CGAL::Tag_true,Point> BASECLASS;

	private:
		Normal normal;
		PlanePoint planePoint;
		int cornerID;
		int index;
		int vertex_type;
		int edgeCount;
		double m_u;
		double m_v;
		double m_gaussian_curvature;
		double m_1D_u;

		//this flag speeds up exact geodesic algorithm
		bool m_saddle_or_boundary;		//it is true if total adjacent angle is larger than 2*PI or this vertex belongs to the mesh boundary
	public:
		double m_distance;
		int m_node_index;
		Point m_laplace_position;
		bool visited;
	

	public:

		My_vertex() 
		{
			Normal n(0,0,0);
			normal = n;
			Point pop(0,0,0);
			PlanePoint plp(0,0);
			planePoint = plp;
			index = UNDEFINED; 
			vertex_type = UNDEFINED; 
			edgeCount = UNDEFINED; 
			cornerID = UNDEFINED; 
			m_u = 0.0f;
			m_v = 0.0f;
			m_type = VERTEX_TYPE;
		}

		bool& saddle_or_boundary(){return m_saddle_or_boundary;};

		~My_vertex() {}

		My_vertex(const My_vertex& v) :BASECLASS(v.point()) 
		{ 
			normal = v.normal; 
			index = v.index;
			vertex_type = v.vertex_type;
			cornerID = v.cornerID;
			edgeCount = v.edgeCount;
			planePoint = v.planePoint;
			m_u = v.m_u;
			m_v = v.m_v;       
			m_type = v.m_type;
		}

		My_vertex(const Point_3 & pt):
		CGAL::HalfedgeDS_vertex_base < Refs, Tag, Point_3 > (pt)
		{
			Normal n(0,0,0);
			normal = n;
			Point pop(0,0,0);
			PlanePoint plp(0,0);
			planePoint = plp;
			vertex_type = UNDEFINED;
			index = UNDEFINED; 
			edgeCount = UNDEFINED; 
			cornerID = UNDEFINED; 
			m_type = VERTEX_TYPE;
		}
		bool SetNormal( Normal& n ) { normal = n; return true; }
		Normal& GetNormal() { return normal; }
		bool SetPlanePoint( PlanePoint& plp ) { planePoint = plp; return true; }
		PlanePoint& GetPlanePoint() { return planePoint; }
		bool SetIndex( int in ) { index = in; return true; }
		int GetIndex() { return index; }
		bool SetCornerID( int id ) { cornerID = id; return true; }
		int GetCornerID() { return cornerID; }
		bool SetType( int t ) { vertex_type = t; return true; }
		int GetType() { return  vertex_type; }
		bool SetedgeCount( int ec ) { edgeCount = ec; return true; }
		int GetedgeCount() { return edgeCount; }
		bool SetUV(double u, double v) { m_u = u; m_v = v;  return true;}
		double GetU()
		{
			return m_u;	
		}
		double GetV()
		{
			return m_v;
		}
		bool SetU(double u)
		{
			m_u = u;
			return true;
		}

		double Get1DU()
		{
			return m_1D_u;
		}
		bool Set1DU(double u)
		{
			m_1D_u = u;
			return true;
		}
		double distance(Point_3 p)
		{
			Point_3 cp = point();
			double dx = cp[0] - p[0];
			double dy = cp[1] - p[1];
			double dz = cp[2] - p[2];
			return sqrt(dx*dx + dy*dy + dz*dz);
		}

		double GetGaussianCurvature() {return m_gaussian_curvature;}
		bool SetGaussianCurvature(double gc) { m_gaussian_curvature = gc; return true;}
	};

}