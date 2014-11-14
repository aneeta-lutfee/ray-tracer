/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"

bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0), 
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
	
	// ray_objectSpace = worldToModel * (*ray);
	
	//Vector4D p = Vector4D(ray.origin[0], ray.origin[1], ray.origin[2], 1);
	//Vector4D v = Vector4D(ray.dir[0], ray.dir[1], ray.dir[1], 0);
	
	// transform the ray into object space
	ray.dir = worldToModel * ray.dir;
	ray.origin = worldToModel * ray.origin;
	
	ray.intersection.t_value = - ray.origin[2] / ray.dir[2];
	double a = ray.origin[0] + ray.intersection.t_value * ray.dir[0];
	double b = ray.origin[1] + ray.intersection.t_value * ray.dir[1];
	
	if (
		((a <= 0.5) && (a >= -0.5)) &&
		((b <= 0.5) && (b >= -0.5))
		)	
	{
		ray.intersection.point[0] = a; 
		ray.intersection.point[1] = b;
		ray.intersection.point[2] = 0;
		
		ray.intersection.normal[0] = 0;
		ray.intersection.normal[1] = 0;
		ray.intersection.normal[2] = 1;
		
		ray.intersection.none = false;
	}
	
	// put the ray back in world space
	ray.dir = modelToWorld * ray.dir;
	ray.origin = modelToWorld * ray.origin;
	ray.intersection.normal = modelToWorld * ray.intersection.normal;
	ray.intersection.point = modelToWorld * ray.intersection.point;
	
	if (ray.intersection.none)
		return false;
	return true;
}

bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSphere, which is centred 
	// on the origin.  
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
	
	return false;
}

