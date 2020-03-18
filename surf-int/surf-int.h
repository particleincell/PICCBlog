/*
PIC-C Blog: Surface Interaction Test
https://www.particleincell.com/2017/line-triangle-intersection/
Written by: Lubos Brieda, 12/20/2017
*/

#ifndef _SURF_INT_H
#define _SURF_INT_H

#include <iostream>
#include <string>

//returns random number in [0,1)
double rnd();
extern const double PI;

inline void warn(std::string msg) {std::cout<<"WARNING: "<<msg;}

/*vector utils*/
namespace vec
{
	inline void clear(double r[3]) {r[0]=0;r[1]=0;r[2]=0;}
	inline void copy(double r[3], double s[3]) {r[0]=s[0];r[1]=s[1];r[2]=s[2];}
	
	inline void add(double r[3],double v1[3], double v2[3])
	{
		for (int i=0;i<3;i++) r[i]=v1[i]+v2[i];
	}

	inline void subtract(double r[3],double v1[3], double v2[3])
	{
		for (int i=0;i<3;i++) r[i]=v1[i]-v2[i];
	}

	inline void mult(double r[3], double v1[3], double s) 
	{
		for (int i=0;i<3;i++) r[i] = v1[i]*s;
	}

	inline void mult_and_add(double r[3], double v1[3], double v2[3], double s) 
	{
		for (int i=0;i<3;i++) r[i] = v1[i] + v2[i]*s;
	}

	inline void unit(double r[3], double v1[3])
	{
		double v_mag = sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
		for (int i=0;i<3;i++) r[i]=v1[i]/v_mag;
	}	
	
	inline void cross(double r[3], double v1[], double v2[])
	{
		r[0] = v1[1]*v2[2]-v1[2]*v2[1];
		r[1] = -v1[0]*v2[2]+v1[2]*v2[0];
		r[2] = v1[0]*v2[1]-v1[1]*v2[0];
	}

	inline double dot(double a[3], double b[3]) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
	inline double mag(double a[3]) {return sqrt(dot(a,a));}

};


#endif