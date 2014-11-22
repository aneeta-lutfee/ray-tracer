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
	double lambda = - ray.origin[2] / ray.dir[2];
	
	if ((!ray.intersection.none )&&  (lambda> ray.intersection.t_value))
	{
		// put the ray back to world coordinates
		ray.dir = modelToWorld * ray.dir;
		ray.origin = modelToWorld * ray.origin;
		return false;
	}
    

	ray.intersection.t_value = - ray.origin[2] / ray.dir[2];
	double a = ray.origin[0] + ray.intersection.t_value * ray.dir[0];
	double b = ray.origin[1] + ray.intersection.t_value * ray.dir[1];
	
	// determine it there is an intersection
	if (
		((a <= 0.5) && (a >= -0.5)) &&
		((b <= 0.5) && (b >= -0.5))
		)	
	{
		// Set intersection point
		ray.intersection.point[0] = a; 
		ray.intersection.point[1] = b;
		ray.intersection.point[2] = 0;
		
		// Set intersection normal
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
	

	/* This derivation is taken from 
	 * http://ray-tracer-concept.blogspot.ca/2011/11/ray-sphere-intersection.html
	 * 
	 * r: sphere radius = 1
	 * c: sphere origin (0,0,0)
	 * Sphere equation: (x-xc)^2 + (y-yx)^2 + (x-xc)^2 = r^2
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
	 * 		Calculate q = b^2 - 4ac:
	 * 		if q < 0 ==> No intersection
	 * 		if q > 0 ==> 2 intersections
	 *  	if q = 0 ==> 1 intersection
	 * 
	 * Roots are: (-b +- sqrt(d))/2a
	 * 
	 * */		
		
	// Transform direction vector to object space	
	Vector3D t_dir = worldToModel * ray.dir;
	

	// Transform origin point to object space
	Point3D t_origin = worldToModel * ray.origin;
	Vector3D t_origin_v = Vector3D(t_origin[0], t_origin[1], t_origin[2]);

	float a = t_dir.dot(t_dir); 
	float b = (t_dir.dot(t_origin_v));
	float c = (t_origin_v.dot(t_origin_v)) - 1;
	float q = pow(b, 2) - (a * c);
	
	float lambda;
	
	if (q >= 0) // There is intersection
	{
		// Intersection points 
		float front_lambda = (-b + sqrt(q)) / (2*a); 
		float back_lambda = (-b - sqrt(q)) / (2*a);
		
		if (front_lambda > 0)
			lambda = front_lambda;
		else if (back_lambda > 0)
			lambda = back_lambda;

		if ((!ray.intersection.none )&&  (lambda> ray.intersection.t_value))
		{
			// put the ray back to world coordinates
			ray.dir = modelToWorld * ray.dir;
			ray.origin = modelToWorld * ray.origin;
			return false;
		}	
		ray.intersection.t_value  = lambda;
		ray.intersection.none = false;
		
		// intersection point p(lambda) = c + lambda(pw - c)
		Vector3D intersection_v = t_origin_v + (lambda * t_dir);
		Point3D intersection_p = Point3D(intersection_v[0], 
										intersection_v[1],
										intersection_v[2]);
										
		ray.intersection.point = modelToWorld * intersection_p;
		ray.intersection.normal = transNorm(modelToWorld, intersection_v);
	}
	
	if (ray.intersection.none)
		return false;
	return true;
}

