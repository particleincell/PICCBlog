/*
PIC-C Blog: Surface Interaction Test
https://www.particleincell.com/2017/line-triangle-intersection/
Written by: Lubos Brieda, 12/20/2017
*/

#include <math.h>
#include "Surface.h"
#include "surf-int.h"


/** Computes normal vectors*/
void SurfaceElement::init(Surface &surf)
{
	/*compute vertex data, the first "k" index indicates primary[0] or secondary[1] triangle (on quads)*/
	int vlut[2][3] = {{0,1,2},{0,2,3}};
	
	//grab node positions;
	double V[4][3];
	for (int i=0;i<num_nodes;i++)
		vec::copy(V[i],surf.nodes[connect[i]].pos);

	for (int k=0;k<2;k++)
	{
		if (num_nodes == 3 && k==1) break;
		for (int i=0;i<3;i++)
		{
			int ip=i+1; if (ip>2) ip=0;
			int im=i-1; if (im<0) im=2;
		
			double *xc = V[vlut[k][i]];
			double *xm = V[vlut[k][im]];
			double *xp = V[vlut[k][ip]];
			vec::copy(vertex_pos[k][i],xc);

			vec::subtract(vertex_vm[k][i],xm,xc);
			vec::subtract(vertex_vp[k][i],xp,xc);
			vec::unit(vertex_vm[k][i],vertex_vm[k][i]);
			vec::unit(vertex_vp[k][i],vertex_vp[k][i]);
			
			double cos_vertex_angle = vec::dot(vertex_vp[k][i],vertex_vm[k][i]);
			vertex_angle[k][i] = acos(cos_vertex_angle);

			//vertex normal
			vec::cross(vertex_normal[k][i],vertex_vp[k][i],vertex_vm[k][i]);
		}
	}

	//set normal vector using 0-1 and 0-2 edges
	vec::clear(normal);
	double v1[3];
	double v2[3];
	vec::subtract(v1,V[1],V[0]);
	vec::subtract(v2,V[2],V[0]);
	vec::cross(normal,v1,v2);
	vec::unit(normal,normal);
		
	/*tang1 is given by the first edge (one used to compute normal)*/
	vec::unit(tang1,v1);

	/*tangent2, normal cross tang1*/
	vec::cross(tang2,normal,tang1);
			
	/*sanity check*/
	double dot1 = vec::dot(tang1,normal);
	double dot2 = vec::dot(tang2,normal);
	double dot3 = vec::dot(tang1,tang2);
	if (fabs(dot1)>1e-4 || fabs(dot2)>1e-4 || fabs(dot3)>1e-4)
		warn("Incorrect tangent vectors computed for el "+std::to_string(id));
	
	/*quads can be slightly warped, need to store secondary triangle data to prevent particle leaks*/
	if (num_nodes==4)
	{
		//compute secondary normal for 023
		vec::subtract(v1,V[2],V[0]);
		vec::subtract(v2,V[3],V[0]);
		vec::cross(sec_normal,v1,v2);
		vec::unit(sec_normal,sec_normal);
		
		//set secondary tanget
		/*tang1 is given by the first edge (one used to compute normal)*/
		vec::unit(sec_tang1,v1);

		/*tangent2, normal cross tang1*/
		vec::cross(sec_tang2,sec_normal,sec_tang1);
		
		/*sanity check*/
		double dot1 = vec::dot(sec_tang1,sec_normal);
		double dot2 = vec::dot(sec_tang2,sec_normal);
		double dot3 = vec::dot(sec_tang1,sec_tang2);
		if (fabs(dot1)>1e-4 || fabs(dot2)>1e-4 || fabs(dot3)>1e-4)
			warn("Incorrect tangent vectors computed for el "+std::to_string(id));			
	}
		
	/*compute area*/
	double vx[3];
	vec::subtract(v1,V[1],V[0]);
	vec::subtract(v2,V[2],V[0]);
	vec::cross(vx,v1,v2);
	area = 0.5*vec::mag(vx);
	if (num_nodes==4)
	{
		vec::subtract(v1,V[2],V[0]);
		vec::subtract(v2,V[3],V[0]);
		vec::cross(vx,v1,v2);
		area += 0.5*vec::mag(vx);
	}

	//check for zero area elements, will be removed by calling loop
	if (area<=0.0)
	{		
		warn("Found element id = "+std::to_string(id)+" with zero area");		
		return;
	}
}


/** random position on a element
per http://mathworld.wolfram.com/TrianglePointPicking.html need to pick
random point in a quad, then discard points that are outside the triangle
returns true if primary, false if secondary
*/
bool SurfaceElement::randomPos(double *pos)
{
	//set k=0 if primary triangle, 1 if secondary
	int k = 0;
	if (num_nodes==4 && rnd()>0.5)
		k = 1;
	
	bool primary = (k==0);

	double (*verts)[3] = vertex_pos[k];
	
	double v1[3], v2[3];
	vec::subtract(v1,verts[1],verts[0]);
	vec::subtract(v2,verts[2],verts[0]);
	
	while (true)
	{	
		/*pick two random numbers*/
		double a1 = rnd();
		double a2 = rnd();

		for (int dim=0;dim<3;dim++) pos[dim] = verts[0][dim] + a1*v1[dim] + a2*v2[dim];
	
		/*check if point is in triangle, we have 50/50 chance*/
		if (containsPoint(pos,primary))	break;
	};

	return primary;
}

/*samples a random Lambertian vector off this surface element
based on https://www.particleincell.com/2015/cosine-distribution/*/
void SurfaceElement::lambertianVector(double v[3], bool primary)
{
	bool pick_another = false;

	//don't have secondary data on a triangle
	if (num_nodes==3) primary = true;

	double *normal = primary ? this->normal : this->sec_normal;
	double *tang1 = primary ? this->tang1 : this->sec_tang1;
	double *tang2 = primary ? this->tang2 : this->sec_tang2;

	do {	
		//this (0.999) is roughly 87 degrees, need this offset otherwise skewed quads can get self-intersections
		//can probably relax on triangles
		//todo: eventually need to handle cases of particle self-intersecting with originating quad other half triangle
		double sin_theta_max  = 0.999;		

		/*sample angle from cosine distribution*/
		double sin_theta = sqrt(rnd());	/*positive half*/
		if (sin_theta>sin_theta_max) sin_theta=sin_theta_max;	
		double cos_theta = sqrt(1-sin_theta*sin_theta);
    
		double vn[3],vt1[3],vt2[3];
		vec::mult(vn,normal,cos_theta);

		/*random rotation about surface normal*/
		double phi = 2*PI*rnd();
		vec::mult(vt1,tang1,sin_theta*cos(phi));
		vec::mult(vt2,tang2,sin_theta*sin(phi));

		/*add components and mulitply by magnitude*/
		for (int i=0;i<3;i++)
			v[i] = vn[i]+vt1[i]+vt2[i];

		/*make sure we don't self-intersect with either normal*/	
		pick_another = false;
		if (vec::dot(v,this->normal)<=0 || (num_nodes==4 && vec::dot(v,this->sec_normal)<=0))
			pick_another = true;
	} while (pick_another);	
}

//wrapper for isPointInTriangle check
bool SurfaceElement::containsPoint(double p[3], bool primary)
{
	double eps = 1e-4;
	bool check1 = isPointInTriNormals(p,primary,eps);
	bool check2;
	
	if (!check1) check2 = isPointInTriAngle(p,primary,eps);

	return check1 || check2;
}

/** line intersection from x1 to x2
Returns -1 if no intersection
On triangle, returns [0,1] if intersection
On quad, returns [0,1] if intersection using primary normal, or [2,3] if with secondary normal (skewed quads)

*/
double SurfaceElement::lineIntersect(double x1[3], double x2[3], double t_min)
{
	double x21[3],p[3];
	double t;

	//first check 
	double *node_pos = vertex_pos[0][0];
	t=SurfUtils::getLPIt(normal,node_pos,x1,x2);

	//check if inside this triangle
	if (t>=0 && t<t_min) 
	{	
		/*compute xyz coordinate of the intersection point*/
		vec::subtract(x21,x2,x1);
		vec::mult_and_add(p,x1,x21,t);
	
		bool in = containsPoint(p,true);
		if (!in) t = -1;
	}		
	
	//if this is a quad, check for another intersection using secondary triangle 0 2 3	
	if (num_nodes==4) 
	{
		double t2 = SurfUtils::getLPIt(sec_normal,node_pos,x1,x2);
		if (t2>=0 && t2<t_min && (t<0 || (t>0 && t2<t)) )
		{
			//check for inside the triangle
			bool in = false;
			/*compute xyz coordinate of the intersection point*/
			vec::subtract(x21,x2,x1);
			vec::mult_and_add(p,x1,x21,t2);
				
			if (containsPoint(p,false))
				return 2.0+t2;  //add 2 to indicate secondary normal
		}		
	}	//if t2<t_min
	
	return t;
}


/*===========================================================================
determines if point is in a triangle by checking angles for edges
angle between edges and line connecting edge to test point is compared to the
angle between the edges

--note same algorithm could be applied for quads but since they are not 
always planar this is now hardcoded for triangles
===========================================================================*/
bool SurfaceElement::isPointInTriAngle(double p[3], bool primary, double eps)
{
	int k = primary?0:1;

	if (num_nodes==3) k=0;

	double (*vm)[3] = this->vertex_vm[k];
	double (*vp)[3] = this->vertex_vp[k];
	double (*xc)[3] = this->vertex_pos[k];		
	double *vertex_angle = this->vertex_angle[k];

	/*loop through vertices*/
	for (int i=0;i<3;i++)
	{		
		//test for point on plane

		double v1[3];
		vec::subtract(v1,p,xc[i]);
		vec::unit(v1,v1);

		/*now compute the angle including the point*/
		double cos_angle1=vec::dot(vm[i],v1);
		if (cos_angle1>1.0) cos_angle1=1.0;
		if (cos_angle1<-1.0) cos_angle1=-1.0;
		double angle1 =acos(cos_angle1);
		if (angle1>vertex_angle[i]+eps) return false;		//need fairly large tolerance here but that gives many stuck particles??

		double cos_angle2=vec::dot(vp[i],v1);
		if (cos_angle2>1.0) cos_angle2=1.0;
		if (cos_angle2<-1.0) cos_angle2=-1.0;
		double angle2 = acos(cos_angle2);
		if (angle2>vertex_angle[i]+eps) return false;
	}
	return true;
}


bool SurfaceElement::isPointInTriAngleSum(double pos[3], bool primary, double eps)
{

}


/*checks for point in triangle by computing normal vector using point at each vertex and comparing to triangle normal*/
bool SurfaceElement::isPointInTriNormals(double pos[3], bool primary, double eps)
{
	int k = primary?0:1;

	if (num_nodes==3) k=0;

	double (*vp)[3] = this->vertex_vp[k];
	double (*xc)[3] = this->vertex_pos[k];	
	double (*normal) = primary?this->normal:this->sec_normal;
	
	for (int i=0;i<3;i++)
	{
		double r[3];
		double norm_p[3];
	
		vec::subtract(r,pos,xc[i]);
		vec::unit(r,r);		//not needed, just to standarize error

		vec::cross(norm_p,vp[i],r);
		double d = vec::dot(norm_p,normal);

		if (d<-eps) return false;		
	}

	return true;

}

/*pushes off point away from edges*/
void SurfaceElement::pushOffEdge(double pos[3], bool primary)
{
	//0 on primary, 1 on secondary
	int k = primary?0:1;	

	double T[3];
	computeBasis(pos,T,primary);
					
	for (int i=0;i<3;i++)
		if (T[i]<0.01)
		{						
			//ray from vertex to point, p = v + r
			double r[3];
			vec::subtract(r,pos,vertex_pos[k][i]);

			//hopefully particle is not so far outside the triangle that this still leaves it out, should check
			vec::mult_and_add(pos,vertex_pos[k][i],r,0.999);	
			break;
		}
}

/*computes basis functions from fractional areas*/
void SurfaceElement::computeBasis(double pos[3], double T[3], bool primary)
{
	int k=primary?0:1;

	//areas are given by 0:(1p)x(2p), 1:(2p)x(0p), 2:(0p)x(1p)
	for (int i=0;i<3;i++)
	{
		int ip = i+1;
		if (ip>2) ip=0;
		int im = i-1;
		if (im<0) im=2;

		double r1[3],r2[3];
		vec::subtract(r1,vertex_pos[k][ip],pos);
		vec::subtract(r2,vertex_pos[k][im],pos);

		double c[3];
		vec::cross(c,r2,r1);
		T[i] = 0.5*vec::mag(c)/area;	
	}

}


/*---------------------------
getLPIt

Performs a line-plane intersections and returns parametric position of the 
intersection point.
x0 and normal define the plane, x1 and x2 are line end points

\return
t:	parametric t for the intersection point, x=x1+t(x2-x1) [0:1)
	t<0 if no intersection
*/
double SurfUtils::getLPIt(double normal[3], double x0[3],double x1[3], double x2[3])
{
	double v[3];
	double t;
	double x21[3];
	double mag;
	double den;
	
	/*get vector connecting the two line end points,x21=x2-x1*/
	vec::subtract(x21,x2,x1);
	mag = vec::mag(x21);

	/*calculate denominator, n*(x2-x1)*/
	den=vec::dot(normal,x21);

	/*den~0 indicates the line is parallel to the plane, so no intersection*/
	if (fabs(den)/mag<SurfUtils::EPS) return -1;
 
	/*get parametric t of the intersection point, t=[n*(x0-x1)]/[n*(x2-x1)]*/
	vec::subtract(v,x0,x1);
	t=vec::dot(normal,v)/den;
	
	/*check bounds, t=[0,1] if the intersection is within the two end points
	using a small delta to account for numerical errors, this needs to be very small like 1e-20 otherwise
	can get fake intersections if particle moving only a tiny distance in time step*/
	if (t<-SurfUtils::LINE_INTERSECT_EPS || t>1+SurfUtils::LINE_INTERSECT_EPS)
		return -1;
 
	/*adjust floating point errors*/
	if (t<0) t=0;
	if (t>1) t=1.0;

	return t;
}
