/*
PIC-C Blog: Surface Interaction Test
https://www.particleincell.com/2017/line-triangle-intersection/
Written by: Lubos Brieda, 12/20/2017
*/

#ifndef _SURFACE_H
#define _SURFACE_H

#include "surf-int.h"
#include <string>
#include <vector>

//node
struct Node
{
	double pos[3];
	Node (double x, double y, double z) {pos[0]=x;pos[1]=y;pos[2]=z;}
};

namespace SurfUtils
{
	double getLPIt(double normal[3], double x0[3],double x1[3], double x2[3]);
	const double EPS = 1e-7;	/*tolerance*/
	const double LINE_INTERSECT_EPS = 1e-7;	/*tolerance for line intersection, [+EPS:1+EPS] is valid*/	
};

class Surface;

/*surface element*/
class SurfaceElement
{
public:
	SurfaceElement () : SurfaceElement(-1,-1,-1,-1) {}
	SurfaceElement (int n1, int n2, int n3, int n4=-1) {
		connect[0]=n1;connect[1]=n2;connect[2]=n3;connect[3]=n4;
		num_nodes = 4 ? connect[3]>=0 : 3;
		normal[0] = 0;
		normal[1] = 0;
		normal[2] = 0;
		num_hits = 0;
	}

	//computes normal vectors
	void init(Surface &surf);

	//computes random point on this element, 
	//returns true if primary, false if secondary
	bool randomPos(double x[3]);
	
	///computes random Lambertian vector, uses point location to determine if secondary normal should be used on quads (ignored on tris)*/
	void lambertianVector(double v[3], bool primary);

	//checks if point is inside the triangle. For quads can specify to check against secondary triangle
	inline bool containsPoint(double pos[3], bool primary);

	//return t=[0,1] or t=[2,3] if segment from x1 to x2 intersect the element
	//values >=2 indicate impact with secondary triangle in quads
	double lineIntersect(double x1[3], double x2[3], double t_min=99);

	//recomputes point location so it's not on an edge
	void pushOffEdge(double pos[3], bool primary);

	//computes basis functions from fractional areas
	void computeBasis(double pos[3], double T[3], bool primary);

	//element data
	int connect[4];
	int num_nodes;
	int id;
	double area;

	//simulation data
	int num_hits;
	
	//geometry	
	double normal[3];
	double tang1[3];
	double tang2[3];

	double sec_normal[3];
	double sec_tang1[3];
	double sec_tang2[3];

	double vertex_angle[2][3];		//vertex angle 
	double vertex_vm[2][3][3];		//normalized vm vector 
	double vertex_vp[2][3][3];		//normalized vp vector 
	double vertex_pos[2][3][3];		//node positions;
	double vertex_normal[2][3][3];	//normal vector at each vertex
	
	/*** various is point in triangle methods ***/
	inline bool isPointInTriAngle(double pos[3], bool primary, double eps);
	inline bool isPointInTriAngleSum(double pos[3], bool primary, double eps);
	inline bool isPointInTriNormals(double pos[3], bool primary, double eps);
	inline bool isPointInTriBasis(double pos[3], bool primary, double eps);
	inline bool isPointInTriArea(double pos[3], bool primary, double eps);
};

struct SurfaceZone
{
public:
	SurfaceZone():SurfaceZone("") {}
	SurfaceZone(std::string name):name(name) { /* */}
	std::vector<SurfaceElement> elements;

	void addElement(SurfaceElement ele) {elements.emplace_back(ele);}
	void init(Surface &surf) {for (SurfaceElement &ele:elements) {ele.init(surf);area+=ele.area;}}
	
	std::string name;
	double area = 0;
};

/*definition of a surface*/
class Surface
{
public:
	std::vector <Node> nodes;
	std::vector <SurfaceZone> surface_zones;

	void addNode(double x, double y, double z) {nodes.emplace_back(x,y,z);}

	//adds zone if not already present
	SurfaceZone& addZone(std::string name){
		SurfaceZone *sz = findZone(name);
		if (sz) return *sz;
		
		surface_zones.push_back(SurfaceZone(name));
		return surface_zones.back();
	}
	
	SurfaceZone& addZone(SurfaceZone &sz){
		surface_zones.push_back(sz);
		return surface_zones.back();
	}
	
	SurfaceZone* findZone(std::string name) {
		for (SurfaceZone &z:surface_zones) {
			if (!name.compare(z.name)) return &z; 
		}
		return nullptr;
	}

	//returns true if point inside box given some tolerance
	bool inBoundingBox(double pos[3], double tol = 1e-6)
	{
		for (int i=0;i<3;i++) 
			if (pos[i]<(bbox_min[i]-tol) || pos[i]>(bbox_max[i]+tol)) return false;
		return true;
	}

	//placeholder to retrieve elements within a bounding box, returns all elements
	std::vector<SurfaceElement*> getElements(double x0[3], double x1[3]) { return all_elements;}

	//initializes all zones;
	void init() {
		for (SurfaceZone &sz:surface_zones)
			sz.init(*this);	

		//build a list of all elements
		all_elements.clear();
		for (SurfaceZone &sz:surface_zones)
			for (SurfaceElement &ele:sz.elements)
				all_elements.push_back(&ele);

		//set bounding box
		vec::copy(bbox_min,nodes[0].pos);
		vec::copy(bbox_max,nodes[0].pos);
		for (Node &node:nodes)
		{
			for (int i=0;i<3;i++)
			{
				if (node.pos[i]<bbox_min[i]) bbox_min[i] = node.pos[i];
				if (node.pos[i]>bbox_max[i]) bbox_max[i] = node.pos[i];
			}
		}
	}
protected:
	std::vector<SurfaceElement *> all_elements;
	double bbox_min[3], bbox_max[3];
};


#endif