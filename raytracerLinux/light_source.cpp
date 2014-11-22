/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

void PointLight::shade( Ray3D& ray ) {
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.
	Vector3D n = ray.intersection.normal;
    Vector3D s = _pos - ray.intersection.point;
    Vector3D r = -s+ 2 *(n.dot(s)) * n;
	Vector3D b = -ray.dir;
	
    s.normalize();
    n.normalize();
    r.normalize();
	b.normalize();
	Material * mat  = ray.intersection.mat;
/*
	Colour ambient = mat->ambient;
	Colour difuse = std::max(0.0, n.dot(s)) * mat->diffuse;
	Colour specular = pow((std::max(0.0, (-d).dot(m))), mat->specular_exp ) * mat->specular;
	ray.col = ambient + difuse + specular;
*/	
	ray.col = mat->ambient * _col_ambient + 
			std::max(0.0, n.dot(s)) * mat->diffuse * _col_diffuse +
			pow(std::max(0.0, r.dot(b)), mat->specular_exp) * mat->specular * _col_specular ;
			
	//clipping
	if (ray.col[0] > 1)
	{
			ray.col[0] = 1;
	}
	if (ray.col[1] > 1)
	{
			ray.col[1] = 1;
	}
	if (ray.col[2] > 1)
	{
			ray.col[2] = 1;
	}


}
