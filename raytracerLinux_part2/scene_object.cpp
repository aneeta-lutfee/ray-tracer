/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"
#if 1
// alanwu
bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0), 
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).

	// ray_objectSpace = worldToModel * (*ray);

	// transform the ray into object space
	Vector3D ray_dir = worldToModel * ray.dir;
	Point3D ray_origin = worldToModel * ray.origin;
	double lambda = - ray_origin[2] / ray_dir[2];
	

	if ((!ray.intersection.none )&&  (lambda> ray.intersection.t_value))
	{
		return false;
	}

	if (lambda >= 0)
	{
		double x = ray_origin[0] + lambda * ray_dir[0];
		double y = ray_origin[1] + lambda * ray_dir[1];
		// determine it there is an intersection
		if (
			((x <= 0.5) && (x >= -0.5)) &&
			((y <= 0.5) && (y >= -0.5))
			)	
		{
			// Set intersection point
			ray.intersection.point[0] = x; 
			ray.intersection.point[1] = y;
			ray.intersection.point[2] = 0;
			ray.intersection.point = modelToWorld * ray.intersection.point;
			
			// Set intersection normal
			ray.intersection.normal = transNorm(worldToModel, Vector3D(0.0, 0.0, 1.0));
			
			// set t_value
			ray.intersection.t_value = lambda;
			ray.intersection.none = false;
		}
	}
	
	if (ray.intersection.none)
		return false;
	return true;
}
#else
#endif

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
	

	/* This derivation is taken from textbook page: 76 and
	 * http://ray-tracer-concept.blogspot.ca/2011/11/ray-sphere-intersection.html
	 * 
	 * r: sphere radius = 1
	 * c: sphere origin (0,0,0)
	 * Sphere equation: (x-xc)^2 + (y-yc)^2 + (x-xc)^2 = r^2
	 * In parametric form: (p-c).(p-c) = r^2 
	 * 		where p is a point on the sphere and c is the center
	 * Ray equation in parametric form: p(t) = e + td
	 * 		e: camera position
	 * 		d: ray direction
	 * 		t: unknown
	 * 
	 * Replacing ray equation in shpere equation
	 * 		((e+td)-c).((e+td)-c) - r^2 = 0
	 * 	==> ((e-c)+td).((e-c)+td) - r^2 = 0
	 * 	==>	((e-c).(e-c)) + ((e-c)td.(e-c)td) + (td.td) - r^2 = 0
	 * 	==>	((e-c).(e-c)) + 2d.(e-c)t + t^2(d.d) - r^2 = 0
	 * 	
	 * Rewrite the equation as a quadratic formula: at^2 + bt + c = 0
	 * a = d.d
	 * b = 2d.(e-c)
	 * c = (e-c).(e-c) - r^2
	 * To determine if the ray intersects with the shpere or no
	 * 		Calculate discriminant q = b^2 - 4ac:
	 * 		if q < 0 ==> No intersection
	 * 		if q > 0 ==> 2 intersections
	 *  	if q = 0 ==> 1 intersection
	 * 
	 * Roots are: (-b +- sqrt(d))/a
	 * 
	 * */		
		
	Point3D sphereOrigin(0,0,0);
	float lambda;
	
	// Transform direction vector to object space	
	Vector3D t_dir = worldToModel * ray.dir;
	
	// Transform origin point to object space
	Point3D t_origin = worldToModel * ray.origin;
	Vector3D t_origin_v = Vector3D(t_origin[0] - sphereOrigin[0], 
								   t_origin[1] - sphereOrigin[1],
								   t_origin[2] - sphereOrigin[2]);

	float a = t_dir.dot(t_dir); 
	float b = (t_dir.dot(t_origin_v));
	float c = (t_origin_v.dot(t_origin_v)) - 1;
	float q = pow(b, 2) - (a * c);

	if (q >= 0) // There is intersection
	{
		// Intersection points 
		float front_lambda = (-b + sqrt(q)) / a; 
		float back_lambda = (-b - sqrt(q)) / a;
		
		lambda = front_lambda;
		if (back_lambda >= 0)
			lambda = back_lambda;
		else 
			return false;

		if ((!ray.intersection.none) && (lambda > ray.intersection.t_value))
		{
			// put the ray back to world coordinates
			ray.dir = modelToWorld * ray.dir;
			ray.origin = modelToWorld * ray.origin;
			return false;
		}	

		// intersection point p(lambda) = c + lambda(pw - c)
		Vector3D intersection_v = t_origin_v + (lambda * t_dir);
		Point3D intersection_p = Point3D(intersection_v[0], 
										intersection_v[1],
										intersection_v[2]);
										
		ray.intersection.point = modelToWorld * intersection_p;
		
		// The unit normal is (p-c)/R = intersection point
		ray.intersection.normal = transNorm(worldToModel, intersection_v);
		ray.intersection.t_value  = lambda;
		ray.intersection.none = false;
	
	}
	
	if (ray.intersection.none)
		return false;
	return true;
}

