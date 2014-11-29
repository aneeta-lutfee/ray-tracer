/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

#include <cstdlib>
#include <time.h>

// Enums for rendering mode
enum RenderMode {Phong, Signature, DiffuseAndAmibient};

RenderMode my_render_mode = Phong;

void PointLight::shade( Ray3D& ray ) {
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

	Colour ambient = mat->ambient * _col_ambient;
	Colour diffuse = std::max(0.0, n.dot(s)) * mat->diffuse * _col_diffuse;
	Colour specular = pow(std::max(0.0, r.dot(b)), mat->specular_exp) * mat->specular * _col_specular;
	
	switch (my_render_mode)
	{
		default:
		case Phong:
			ray.col = ambient + diffuse + specular;
			break;
			
		case DiffuseAndAmibient:
			ray.col = ambient + diffuse;
			break;
			
		case Signature:
			
			break;
	}
	
	//clampping
	ray.col.clamp();
}

void AreaLight::shade( Ray3D& ray ) {
	Point3D p = get_position();
	Point3D light_pos;
	Vector3D n = get_normal();
	Vector3D v;
	double r = get_radius();
	double nHigh = r; 
	double nLow = 0;
	
	// (point - p).dot(n) = 0
	// |point - p| < r
	
	// v = point - p
	// v.x * n.x + v.y * n.y + v.z * n*z = 0
	// randomize x and y [0, r], solve for z
	// check sqrt(x^2 + y^2 + z^2) < r
	// point = p + v
	// launch new light at point
	Colour final_colour;
	int num_lights = 20;
	for (int i = 0; i < num_lights; i ++)
	{
		do
		{
			v[0] = ((double) rand() / (RAND_MAX)) * r;
			v[1] = ((double) rand() / (RAND_MAX)) * r;
			v[2] = -(v[0]*n[0] + v[1]*n[1])/n[2];
		} while (v.length() > r);
		light_pos = p + v;
		// create a new point light at light_pos
		PointLight pt_light(light_pos, _col_ambient, _col_diffuse, _col_specular);
		pt_light.shade(ray);
		final_colour = final_colour + ray.col;
	}
	ray.col = (1.0/num_lights) * final_colour;
}
