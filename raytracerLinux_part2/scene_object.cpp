/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "scene_object.h"

// alanwu
bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld, double time) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0), 
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).

	// ray_objectSpace = worldToModel * (*ray);
	
	// transform the ray into object space
	Vector3D ray_dir = worldToModel * ray.dir;
	Point3D ray_origin = worldToModel * ray.origin;
	double lambda = - ray_origin[2] / ray_dir[2];

	//double time = double (rand())/(RAND_MAX) * 0.5;
	time = 0;

	if (lambda >= 0)
	{
		double x = ray_origin[0] + lambda * ray_dir[0];
		double y = ray_origin[1] + lambda * ray_dir[1];
		// determine it there is an intersection
		if (
			((x <= 0.5 + time) && (x >= -0.5 + time)) &&
			((y <= 0.5) && (y >= -0.5))
			)	
		{
			if (ray.intersection.t_value > lambda || ray.intersection.none)
			{
				// Set intersection point
				ray.intersection.point[0] = x; 
				ray.intersection.point[1] = y;
				ray.intersection.point[2] = 0;
				ray.intersection.point = modelToWorld * ray.intersection.point;
				
				// Set intersection normal
				ray.intersection.normal = transNorm(worldToModel, Vector3D(0.0, 0.0, 1.0));
				ray.intersection.normal.normalize();
				// set t_value
				ray.intersection.t_value = lambda;
				// set intersection
				ray.intersection.none = false;
				return true;
			}
		}
	}
	return false;
}

#if 1
//Anita's code
bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld, double time ) {

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
	
	Point3D sphereOrigin(time, 0, 0);
	float lambda;
	
	// Transform direction vector to object space	
	Vector3D t_dir = worldToModel * ray.dir;
	
	// Transform origin point to object space
	Point3D t_origin = worldToModel * ray.origin;
	Vector3D t_origin_v = t_origin - sphereOrigin;

	// Calculate the values of a, b, c, q to substitute in the equation
	float a = t_dir.dot(t_dir); 
	float b = t_dir.dot(t_origin_v);
	float c = t_origin_v.dot(t_origin_v) - 1;
	float q = pow(b, 2) - (a * c);
	
	if (q < 0) // There is no intersection
	{
		return false;
	} else { // There is intersection
		
		// Find lambda at the Intersection points 
		float front_lambda = (-b + sqrt(q)) / a; 
		float back_lambda = (-b - sqrt(q)) / a;
		
		// Set lambda values
		if (front_lambda < 0 && back_lambda < 0) // hits are behind the view plane
			return false;
		else if (front_lambda > 0 && back_lambda < 0) // front_lambda is a valid hit
			lambda = front_lambda;
		else  // 2 valid hits, smallest lambda gives intersection closest to camera
			lambda = std::min(front_lambda, back_lambda);
		}
		if (lambda < 0.01){
			return false;
		}
		// There is already a closer intersection
		if ((!ray.intersection.none) && (lambda > ray.intersection.t_value))
		{
			return false;
		}	

		// intersection point p(lambda) = c + lambda(pw - c)
		Point3D intersection_p = t_origin + lambda * t_dir;										
		ray.intersection.point = modelToWorld * intersection_p; 
		
		// The unit normal is (p-c)/R = intersection point
		Vector3D normal = intersection_p - sphereOrigin;
		ray.intersection.normal = modelToWorld * normal;
		ray.intersection.normal.normalize();
		ray.intersection.t_value  = lambda;
		ray.intersection.none = false;
		return true;
}

#else
bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld, double time ) {
	// TODO: implement intersection code for UnitSphere, which is centred 
	// on the origin.  
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
	
	Point3D rayOrigin = worldToModel * ray.origin;

		Vector3D rayDir = worldToModel * ray.dir;
		Point3D sphereOrigin(time,0,0);
		double lambdaStar;
		//std::cout << "ray origin is " << rayOrigin << "\n";
		//std::cout << "ray direction is " << rayDir << "\n";

		double A = rayDir.dot(rayDir);
		double B = (rayOrigin-sphereOrigin).dot(rayDir);
		double C = (rayOrigin-sphereOrigin).dot(rayOrigin-sphereOrigin) - 1;
		double D = B*B-A*C;

		//std::cout << "A, B, C, D are "<< A << "," << B << "," << C << "," << D << "\n";
		if (D<0)
		{
			return false;
		}else if(D*D < 0.00001){
			double lambda = -B/A;
			if (((lambda-0)*(lambda-0)) < 0.00001){
				return false;
			}else{
				lambdaStar = lambda;
			}
		}
		else
		{
			double lambda_1 = -B/A + sqrt(D) / A;
			double lambda_2 = -B/A - sqrt(D) / A;
			//std::cout << "lambda1  "<< lambda_1 << " lambda2 " << lambda_2 << "\n";
			if (lambda_1 < 0 && lambda_2 < 0)
			{
				return false;
			}
			else if (lambda_1 > 0 && lambda_2 < 0)
			{
				lambdaStar = lambda_1;
			}
			else
			{
				lambdaStar = lambda_2;
			}
		}
		if (lambdaStar < 0.01){
			return false;
		}
		Point3D intersectionPoint = rayOrigin + lambdaStar * rayDir;
		Vector3D normal = intersectionPoint - sphereOrigin;

		if (!ray.intersection.none && lambdaStar > ray.intersection.t_value){
			return false;
		}
		ray.intersection.point = modelToWorld * intersectionPoint;
		//std::cout << "ray origin " << ray.origin << "\n";
		//std::cout << "intersection point " << ray.intersection.point << "\n";
		/*if ((ray.origin - ray.intersection.point).length() < 0.0001){
			return false;
		}*/
		ray.intersection.normal = modelToWorld * normal;
		ray.intersection.normal.normalize();
		ray.intersection.t_value = lambdaStar;
		ray.intersection.none = false;
		return true;

	//return false;
}
#endif
bool UnitSphereStatic::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld, double time ) {
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

	Point3D sphereOrigin(0, 0, 0);
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
		return true;
	}
	return false;
}

