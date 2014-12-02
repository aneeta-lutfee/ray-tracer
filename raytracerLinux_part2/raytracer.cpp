/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		Implementations of functions in raytracer.h, 
		and the main function which specifies the 
		scene to be rendered.	

***********************************************************/


#include "raytracer.h"
#include "bmp_io.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include "util.h"

Raytracer::Raytracer() : _lightSource(NULL) {
	_root = new SceneDagNode();
}

Raytracer::~Raytracer() {
	delete _root;
}

SceneDagNode* Raytracer::addObject( SceneDagNode* parent, 
		SceneObject* obj, Material* mat, Colour ID) {
	SceneDagNode* node = new SceneDagNode( obj, mat );
	node->parent = parent;
	node->next = NULL;
	node->child = NULL;
	node->ID = ID;
	
	// Add the object to the parent's child list, this means
	// whatever transformation applied to the parent will also
	// be applied to the child.
	if (parent->child == NULL) {
		parent->child = node;
	}
	else {
		parent = parent->child;
		while (parent->next != NULL) {
			parent = parent->next;
		}
		parent->next = node;
	}
	
	return node;;
}

LightListNode* Raytracer::addLightSource( LightSource* light) {
	LightListNode* tmp = _lightSource;
	_lightSource = new LightListNode( light, tmp );
	return _lightSource;
}

void Raytracer::rotate( SceneDagNode* node, char axis, double angle ) {
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;
	
	for (i = 0; i < 2; i++) {
		switch(axis) {
			case 'x':
				rotation[0][0] = 1;
				rotation[1][1] = cos(angle*toRadian);
				rotation[1][2] = -sin(angle*toRadian);
				rotation[2][1] = sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'y':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][2] = sin(angle*toRadian);
				rotation[1][1] = 1;
				rotation[2][0] = -sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'z':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][1] = -sin(angle*toRadian);
				rotation[1][0] = sin(angle*toRadian);
				rotation[1][1] = cos(angle*toRadian);
				rotation[2][2] = 1;
				rotation[3][3] = 1;
			break;
		}
		if (i == 0) {
		    node->trans = node->trans*rotation; 	
			angle = -angle;
		} 
		else {
			node->invtrans = rotation*node->invtrans; 
		}	
	}
}

void Raytracer::translate( SceneDagNode* node, Vector3D trans ) {
	Matrix4x4 translation;
	
	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	node->trans = node->trans*translation; 	
	translation[0][3] = -trans[0];
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	node->invtrans = translation*node->invtrans; 
}

void Raytracer::scale( SceneDagNode* node, Point3D origin, double factor[3] ) {
	Matrix4x4 scale;
	
	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];
	node->trans = node->trans*scale; 	
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	node->invtrans = scale*node->invtrans; 
}

Matrix4x4 Raytracer::initInvViewMatrix( Point3D eye, Vector3D view, 
		Vector3D up ) {
	Matrix4x4 mat; 
	Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);

	mat[0][0] = w[0];
	mat[1][0] = w[1];
	mat[2][0] = w[2];
	mat[0][1] = up[0];
	mat[1][1] = up[1];
	mat[2][1] = up[2];
	mat[0][2] = -view[0];
	mat[1][2] = -view[1];
	mat[2][2] = -view[2];
	mat[0][3] = eye[0];
	mat[1][3] = eye[1];
	mat[2][3] = eye[2];

	return mat; 
}

void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray, double time ) {
	SceneDagNode *childPtr;

	// Applies transformation of the current node to the global
	// transformation matrices.
	_modelToWorld = _modelToWorld*node->trans;
	_worldToModel = node->invtrans*_worldToModel; 
	if (node->obj) {
		// Perform intersection.
		if (node->obj->intersect(ray, _worldToModel, _modelToWorld, time)) {
			ray.intersection.mat = node->mat;
			
#ifdef SCENE_SIGNATURE
			ray.col = node->ID;
#endif
		}
		
	}
	// Traverse the children.
	childPtr = node->child;
	while (childPtr != NULL) {
		traverseScene(childPtr, ray, time);
		childPtr = childPtr->next;
	}

	// Removes transformation of the current node from the global
	// transformation matrices.
	_worldToModel = node->trans*_worldToModel;
	_modelToWorld = _modelToWorld*node->invtrans;
}

void Raytracer::computeShading( Ray3D& ray, double time ) {
	LightListNode* curLight = _lightSource;
	
	// light ray used to collect color from each light source
	Ray3D lightRay = ray;
	Colour final_col;
	
	// counter for number of lights
	int counter = 0;
	
	for (;;) {
		if (curLight == NULL) break;

		// alanwu: implement shadows
		// for hard shadows, if object block the ray from intersection to light source, 
		// only set the color to ambient color of the object.
#ifdef SHADOWS
		Ray3D shadowRay;
		shadowRay.dir = curLight->light->get_position() - lightRay.intersection.point ;
		shadowRay.origin = lightRay.intersection.point + 0.0001 * (shadowRay.dir);
		shadowRay.intersection.none = true;
		
		traverseScene(_root, shadowRay, time); 
		if (shadowRay.intersection.none || shadowRay.intersection.mat->refractivity != 0.0)
		{
			// shade only if the light shines the intersection
			curLight->light->shade(lightRay);
			final_col = final_col + lightRay.col;
		}
		else
		{
			// set color to the ambient colour of the object
			//ray.col = ray.intersection.mat->ambient;
			final_col = final_col + lightRay.intersection.mat->ambient;
		}
#else
	curLight->light->shade(lightRay);
	final_col = final_col + lightRay.col;
	
#endif
	curLight = curLight->next;
	counter ++;
	}
	
	// average out the final colour
	ray.col = (1.0/counter) * final_col;
	ray.col.clamp();
}

void Raytracer::initPixelBuffer() {
	int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
	_rbuffer = new unsigned char[numbytes];
	_gbuffer = new unsigned char[numbytes];
	_bbuffer = new unsigned char[numbytes];
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			_rbuffer[i*_scrWidth+j] = 0;
			_gbuffer[i*_scrWidth+j] = 0;
			_bbuffer[i*_scrWidth+j] = 0;
		}
	}
}

void Raytracer::flushPixelBuffer( char *file_name ) {
	bmp_write( file_name, _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer );
	delete _rbuffer;
	delete _gbuffer;
	delete _bbuffer;
}

Colour Raytracer::shadeRay( Ray3D& ray , int depth) {
	Colour col(0.0, 0.0, 0.0); 
	
#ifndef SCENE_SIGNATURE
	double time = 0.0;
	int counter = 0;
	
#ifdef MOTION_BLUR
	time = TIME_DURATION;
#endif // END_MOTION_BLUR

	for (; time >= 0; time-=TIME_SLICE)
	{
		traverseScene(_root, ray, time); 
		// Don't bother shading if the ray didn't hit 
		// anything.
		if (!ray.intersection.none) {
			computeShading(ray, time); 
			col = col + ray.col; 
	
			if (depth < 2)
			{
#ifdef REFLECTIONS	
				// You'll want to call shadeRay recursively (with a different ray, 
				// of course) here to implement reflection/refraction effects. 
		#ifndef GLOSSY
				// ---------------------- Reflection -------------------
				Vector3D d = ray.dir;
				Vector3D n = ray.intersection.normal;
				n.normalize();
				Vector3D reflected_dir =  d - (2.0) * (d.dot(n)) * n;
				reflected_dir.normalize(); 
				Ray3D reflection_ray(ray.intersection.point + 0.00001* reflected_dir, reflected_dir);

				col = col + ray.intersection.mat->reflectivity * shadeRay(reflection_ray ,depth + 1);	
				
		#else	
				// ---------------------- Glossy reflection ------------- 
				Vector3D d = ray.dir;
				Vector3D n = ray.intersection.normal;
				n.normalize();
				Vector3D reflected_dir =  d - (2.0) * (d.dot(n)) * n;
				reflected_dir.normalize(); 
				Ray3D reflection_ray(ray.intersection.point + 0.00001* reflected_dir, reflected_dir);
				Colour main_ref_colour = ray.intersection.mat->reflectivity * shadeRay(reflection_ray ,depth + 1);
				
				//launch rays randomly, similar code to creating the disk light
				//reflected_dir is the disk normal
				double radius = 0.5;
				Vector3D v, new_ref_dir;
				Point3D ref_pos, disk_pos;
				disk_pos = ray.intersection.point + 1 * reflected_dir;
				int num_rays = 200 * radius * radius;
				Colour gloss_color(0,0,0);
				for (int i = 0; i < num_rays; i ++)
				{
					v[0] = ((double) rand() / (RAND_MAX)) * radius;
					v[1] = ((double) rand() / (RAND_MAX)) * radius;
					v[2] = -(v[0]*reflected_dir[0] + v[1]*reflected_dir[1])/reflected_dir[2];
					v.normalize();
					// set the length of v less than radius
					v = ((double) rand() / (RAND_MAX)) * radius * v;
				
					ref_pos = disk_pos + v;
					
					// fire a new reflection ray
					new_ref_dir = ref_pos - ray.intersection.point;
					
					new_ref_dir.normalize();
					Ray3D new_reflection_ray(ray.intersection.point + 0.00001* new_ref_dir, new_ref_dir);
					gloss_color = gloss_color + ray.intersection.mat->reflectivity * shadeRay(new_reflection_ray ,depth + 1);
				}
				gloss_color = (1.0/(1+num_rays)) * (main_ref_colour + gloss_color);
				col = col + gloss_color;

		#endif	//END_GLOSSY
	#endif	// END_REFLECTIONS		
	
		// ------------- REFRACTION ----------------	
		#ifdef REFRACTION 
			if (ray.intersection.mat->refractivity)
			{
				Vector3D inDir = ray.dir;
				inDir.normalize();
				Vector3D normal = ray.intersection.normal;
				normal.normalize();
				
				double cosTheta1 = inDir.dot(normal); // incident theta
				double c1, c2;
				double nr;
				Vector3D T, Tout;
				
				// ray transmisison from air to object 
				if (cosTheta1 < 0) {
					c1 = 1.0; // 
					c2 = ray.intersection.mat->lightSpeed;
					nr = c1/c2;
					float ndotNegv = normal.dot(-inDir);
					float rootContent = 1.0 - pow(nr, 2) * (1-pow(ndotNegv, 2));
					if (rootContent >= 0.0) {
						// Find refracted ray
						T = (nr * ndotNegv - sqrt(rootContent)) * normal - (nr * (-inDir));
						T.normalize();
						Ray3D refractive(ray.intersection.point + 0.009 * T, T);	
						Colour refractiveCol = shadeRay(refractive, depth+1);
						col = col + ray.intersection.mat->refractivity * refractiveCol;
					} 
				} else { // ray transmisison from object to air 
					c1 = ray.intersection.mat->lightSpeed;
					c2 = 1.0;
					nr = c1/c2;
					Vector3D matDir = ray.dir;
					matDir.normalize();
					Vector3D negN = -normal;
					negN.normalize();
					float negNDotNegV = negN.dot(-matDir);
					float rootContent = 1.0 - pow(nr, 2) * (1 - pow(negNDotNegV, 2));
					if (rootContent >= 0.0) {
						// Find refracted ray
						Tout = (nr * negNDotNegV - sqrt(rootContent)) * negN - (nr * (-matDir));
						Tout.normalize();
						Ray3D refractive1(ray.intersection.point + 0.009 * Tout, Tout);	
						Colour refractiveCol1 = shadeRay(refractive1, depth+1);
						col = col + ray.intersection.mat->refractivity * refractiveCol1;
					}
				}
			}
		#endif	// End refraction
			}
		
		}
		counter++;
	}// end of time loop
	col = (1.0/counter) * col;
	col.clamp();

#else // ELSE_SCENCE_SIGNATURE
	traverseScene(_root, ray, 0.0);
	col = ray.col;
#endif

	return col; 
}	

void Raytracer::render( int width, int height, Point3D eye, Vector3D view, 
		Vector3D up, double fov, char* fileName ) {
	Matrix4x4 viewToWorld;
	_scrWidth = width;
	_scrHeight = height;
	double subPixel = 2;
	double factor = (double(height)/2)/tan(fov*M_PI/360.0);
	double accumulate_red; // Colour accumulator for RGB
	double accumulate_green; // Colour accumulator for RGB
	double accumulate_blue; // Colour accumulator for RGB
	
	// alanwu: Depth of Field parameters
	int num_rays = 2;
	double focal_length = 5;
	double jitter = 0.2; // large jitter (~1) cause artifacts
	
	initPixelBuffer();
	viewToWorld = initInvViewMatrix(eye, view, up);


	// Construct a ray for each pixel.
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			// Shoot subPixed^2 rays for each pixel
			for (double sub_y = 0; sub_y < subPixel; sub_y++) {
				for (double sub_x = 0; sub_x < subPixel; sub_x++) {
					// Sets up ray origin and direction in view space, 
					// image plane is at z = -1.
					Point3D origin(0, 0, 0);
					Point3D imagePlane;
					
					double sub_xx = double(1)/double(subPixel*2) + sub_x * double(1)/double(subPixel);
					double sub_yy = double(1)/double(subPixel*2) + sub_y * double(1)/double(subPixel);	
					imagePlane[0] = (-double(width)/2 + sub_xx + j)/factor;
					imagePlane[1] = (-double(height)/2 + sub_yy + i)/factor;
					imagePlane[2] = -1;

					Ray3D ray;
					//alanwu: filling the origin and the direction of the ray
					ray.origin = viewToWorld * (imagePlane);
					ray.dir = viewToWorld * (imagePlane - origin);

					// Get back a color
					Colour col = shadeRay(ray, 0); 
					//printf("col: %f, %f, %f\n", col[0], col[1], col[2]);
					
#ifdef DOF
					// alanwu: Depth of field, initiate multiple rays around this point, and 
					// average them out
					Vector3D f_dir = viewToWorld * (imagePlane - origin);
					f_dir.normalize();
					
					// get the focus_point
					Point3D focus_point = viewToWorld * imagePlane + focal_length * f_dir;
					focus_point;
					int counter = 0;
					for (int m = -num_rays; m <= num_rays; m++)
					{
						for (int n = -num_rays; n <= num_rays; n++)
						{
							Ray3D new_ray;
							Vector3D origin_shift(m*jitter, n*jitter, 0);
							new_ray.origin = origin + origin_shift;
							new_ray.dir = focus_point - viewToWorld * new_ray.origin;
							new_ray.origin = viewToWorld * new_ray.origin;
							
							Colour new_col = shadeRay(new_ray, 0);
							counter += 1;
							col = col + new_col;
						}
					}
					// alanwu: average the color
					col = (1.0/(1+counter)) * col;
#endif
					
					// accumulate separately for each color
					if (sub_y == 0 && sub_x == 0) {
						accumulate_red = col[0];
						accumulate_green = col[1];
						accumulate_blue = col[2];
					}
					else {
						accumulate_red = accumulate_red + col[0];
						accumulate_green = accumulate_green + col[1];
						accumulate_blue = accumulate_blue + col[2];
					}
				}
			}
			
			// Divide the colors by the number of subpixels
			double col_factor = double(1)/double(pow(subPixel, 2));
			
			// Fill up the rgb buffer
			_rbuffer[i*width+j] = int(accumulate_red*col_factor*255);
			_gbuffer[i*width+j] = int(accumulate_green*col_factor*255);
			_bbuffer[i*width+j] = int(accumulate_blue*col_factor*255);
			
			// show progress
			printf("Percentage Complete: %.f\%%           \r", ((i * _scrWidth + j ) * 100.0)/(_scrHeight * _scrWidth));
		}
		
	}
	printf("\nFile Outputed!\n");
	flushPixelBuffer(fileName);
}

int main(int argc, char* argv[])
{	
	// Build your scene and setup your camera here, by calling 
	// functions from Raytracer.  The code here sets up an example
	// scene and renders it from two different view points, DO NOT
	// change this if you're just implementing part one of the 
	// assignment.  
	Raytracer raytracer;
#ifdef TEST_IMAGE
	int width = 32 *4; 
	int height = 24 *4; 
#else
	int width = 320; 
	int height = 240; 
#endif

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}

    // Set random
    srand( time(0) );
	// Camera parameters.
	Point3D eye(0, 0, 1);
	Vector3D view(0, 0, -1);
	Vector3D up(0, 1, 0);
	double fov = 60; // This is changed from 60

	// Defines a material for shading.
	Material gold( Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648), 
			Colour(0.628281, 0.555802, 0.366065), 
			51.2 , 0.5, 0.0, 0.5);
	Material jade( Colour(0, 0, 0), Colour(0.54, 0.89, 0.63), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8 ,0.3, 0.0, 0.5);
	Material blue( Colour(0, 0, 0), Colour(100.0/255, 149.0/255, 237.0/255), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8 , 0.9, 0.0, 0.5);
	Material glass( Colour(0.05, 0.05, 0.05), Colour(0.05, 0.05, 0.05),
			Colour(0.5, 0.5, 0.5),
			32, 0.0, 0.9, 1.2); 

	
#ifdef AREA_LIGHT
	// Defines an area light source by randomly creating many point
	// light sources around a center location on a plane to simulate
	// a disk light panel.
	Point3D area_light_position(0, 0, 5);
	Vector3D area_light_normal(0, 0, -1);
	Colour area_light_col(0.9, 0.9, 0.9);
	double radius = 5.0;
	
	Point3D light_pos;
	Vector3D v;

	// scales with the radius
	int num_lights = 6 * radius*radius;
	for (int i = 0; i < num_lights; i ++)
	{
		do
		{
			v[0] = ((double) rand() / (RAND_MAX)) * radius;
			v[1] = ((double) rand() / (RAND_MAX)) * radius;
			v[2] = -(v[0]*area_light_normal[0] + v[1]*area_light_normal[1])/area_light_normal[2];
		} while (v.length() > radius);
		light_pos = area_light_position + v;
		// create a new point light at light_pos
		PointLight * newlight = new PointLight(light_pos, area_light_col);
		raytracer.addLightSource(newlight);
	}
#else
	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(0, 5, 5), 
				Colour(0.9, 0.9, 0.9) ) );
#endif

#if 1

	// Add a unit square into the scene with material mat.
	SceneDagNode* sphere = raytracer.addObject( new UnitSphere(), &jade, Colour(1,1,1));
	
#ifdef REFRACTION
	SceneDagNode* sphere1 = raytracer.addObject( new UnitSphere(), &glass, Colour(1,1,1));
#endif
	SceneDagNode* sphere2 = raytracer.addObject( new UnitSphere(), &blue, Colour(1,1,1));
	SceneDagNode* plane = raytracer.addObject( new UnitSquare(), &gold, Colour(0,1,0));
//	SceneDagNode* plane1 = raytracer.addObject( new UnitSquare(), &jade );  // Remove: new plane for testing
	
	// Apply some transformations to the unit square.
	double factor1[3] = { 1.0, 1.0, 1.0 };
	double factor2[3] = { 5.0, 5.0, 5.0 };
	double factor3[3] = { 1.5, 1.5, 1.5 };
	raytracer.translate(sphere, Vector3D(-1, 0, -4));	
	raytracer.rotate(sphere, 'x', -45); 
	raytracer.rotate(sphere, 'z', 45); 
	raytracer.scale(sphere, Point3D(0, 0, 0), factor1);
	
	raytracer.translate(sphere2, Vector3D(-1, 2, -4));	
	raytracer.rotate(sphere2, 'x', -45); 
	raytracer.rotate(sphere2, 'z', 45); 
	raytracer.scale(sphere2, Point3D(0, 0, 0), factor1);
	
#ifdef REFRACTION
	raytracer.translate(sphere1, Vector3D(0, 1, -2));	
	raytracer.rotate(sphere1, 'x', -45); 
	raytracer.rotate(sphere1, 'z', 45); 
	raytracer.scale(sphere1, Point3D(0, 0, 0), factor1);
#endif
	raytracer.translate(plane, Vector3D(0, 0, -5));	
	raytracer.rotate(plane, 'z', 45); 
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);
/*
	// Apply some transformations to the unit square.
	double factor1[3] = { 1.0, 1.0, 1.0 };
	double factor2[3] = { 1.5, 1.5, 1.5 };
	double factor3[3] = { 15, 15, 15};
	
	// 3 spheres
	SceneDagNode* plane = raytracer.addObject( new UnitSphere(), &jade, Colour(1,1,1) );
	SceneDagNode* sphere = raytracer.addObject( new UnitSphere(), &gold, Colour(0,1,1) );
	SceneDagNode* sphere2 = raytracer.addObject( new UnitSphere(), &gold, Colour(1,0,1) );

	
	raytracer.translate(sphere, Vector3D(-3.5, 0, -3.5));	
	raytracer.rotate(sphere, 'x', -45); 
	raytracer.rotate(sphere, 'z', 45); 
	raytracer.scale(sphere, Point3D(0, 0, 0), factor1);

	raytracer.translate(plane, Vector3D(0, 0, -5));	
	raytracer.rotate(plane, 'z', 45); 
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);
	
	raytracer.translate(sphere2, Vector3D(-6, -7, -20));
	raytracer.scale(sphere2, Point3D(0, 0, 0), factor3);
	*/
#else
	// Add a unit square into the scene with material mat.
	SceneDagNode* sphere = raytracer.addObject( new UnitSphere(), &jade, Colour(1,0,1));
	SceneDagNode* plane = raytracer.addObject( new UnitSphereStatic(), &blue, Colour(1,1,0) );

//	SceneDagNode* plane1 = raytracer.addObject( new UnitSquare(), &jade );  // Remove: new plane for testing
	
	// Apply some transformations to the unit square.
	double factor1[3] = { 1.0, 1.0, 1.0 };
	double factor2[3] = { 17.0, 17.0, 17.0 };
	raytracer.translate(sphere, Vector3D(0, 0, -5));	
	raytracer.rotate(sphere, 'x', -45); 
	raytracer.rotate(sphere, 'z', 45); 
	raytracer.scale(sphere, Point3D(0, 0, 0), factor1);

	raytracer.translate(plane, Vector3D(0, 0, -25));	
	raytracer.rotate(plane, 'z', 45); 
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);
#endif	
	// Render the scene, feel free to make the image smaller for
	// testing purposes.	
	raytracer.render(width, height, eye, view, up, fov, "view1.bmp");
	
	// Render it from a different point of view.
	Point3D eye2(4, 2, 1);
	Vector3D view2(-4, -2, -6);
	raytracer.render(width, height, eye2, view2, up, fov, "view2.bmp");
	
	return 0;
}

