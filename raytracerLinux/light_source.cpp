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
    Vector3D d = -ray.dir;
    Vector3D r = -s+ 2 *(n.dot(s)) * n;
    
    s.normalize();
	Material * mat  = ray.intersection.mat;
	//ray.col[0] = mat->ambient[0] + mat->diffuse[0] * std::max(0, n * s) + mat->specular[0] * mat->specular_exp * std::max(0, d * r);
	ray.col[0] = mat->ambient[0] + mat->diffuse[0] * std::max(0.0, n.dot(s)) + mat->specular[0] * mat->specular_exp * std::max(0.0, d.dot(r));
	ray.col[1] = mat->ambient[1] + mat->diffuse[1] + mat->specular[1] * mat->specular_exp;
	ray.col[2] = mat->ambient[2] + mat->diffuse[2] + mat->specular[2] * mat->specular_exp;

}

