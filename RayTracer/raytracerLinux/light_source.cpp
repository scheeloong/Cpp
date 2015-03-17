/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/
// Defines the basic light class.  You could define different types of 
// lights, which shades the ray differently.  Point lights are sufficient 
// for most scenes.  

#define RADIUS 1 // assume radius of sphere is 1 
#include <cmath>
#include "light_source.h"

// define in raytracer.cpp 
extern int width; 
extern int height; 
extern	unsigned char* rbuffer; 
extern	unsigned char* gbuffer; 
extern	unsigned char* bbuffer; 

void PointLight::shade( Ray3D& ray, bool isInShadow, int onlyAmbient, int textureMap  ) {
	// Implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  

	if(!ray.intersection.none) {
		// Intersection point to point light source position
		Vector3D s = Vector3D(this->_pos - ray.intersection.point);
		s.normalize();

		// Reflection of s across n
		Vector3D r = Vector3D(-1 * s + 2 * ray.intersection.normal.dot(s) * ray.intersection.normal);

		r.normalize();
		// Intersection point to camera position
		Vector3D b = ray.origin - ray.intersection.point;
		b.normalize();

		// Material properties, easier access
		Colour mat_a =  Colour(ray.intersection.mat->ambient);
		Colour mat_d =  Colour(ray.intersection.mat->diffuse);
		Colour mat_s =  Colour(ray.intersection.mat->specular);
		
		// Light source properties
		Colour I_a =   this->_col_ambient;
		Colour I_d =  this->_col_diffuse;
		Colour I_s =  this->_col_specular;

		double alpha = ray.intersection.mat->specular_exp;
		// if don't do image map 
		//std::cout <<" IMAGE MAP IS " <<ray.intersection.mat->imageMap <<std::endl; 
		// Note: Code always assume image map is applied to spheres only 
		if(ray.intersection.mat->imageMap == 1)
		{
			Point3D IntersectPoint = ray.intersection.point; // get the x, y, z coordinates of the sphere
			// Need center of sphere
			Point3D centerOfSphere = ray.intersection.CenterPoint; 
			// in radians 
			// theta = [0, pi]
			// phi = [-pi, pi]
			double theta = acos((IntersectPoint[2]-centerOfSphere[2])/RADIUS); 
			double phi = atan2(IntersectPoint[1]-centerOfSphere[1],IntersectPoint[0]-centerOfSphere[0]); 
			double u = phi/(double)(2.0*M_PI);
			double v = (M_PI - theta)/(double)M_PI; 
			u += 0.5; 
			// Now both u and v goes from 0 to 1 
			u *= width; 
			v *= height; 
			int temp = floor(v)*width+floor(u); 
			ray.col[0] += (rbuffer[temp])/(double)255 ;
			ray.col[1] += (gbuffer[temp])/(double)255 ;
			ray.col[2] += (bbuffer[temp])/(double)255 ;
		}
		// Texture Mapping for plane 
		else if(ray.intersection.mat->imageMap == 2)
		{
			Point3D IntersectPoint = ray.intersection.point; // get the x, y, z coordinates of the sphere
			// Need center of sphere
			Point3D centerOfPlane = ray.intersection.CenterPoint; 
			// Determine if colour black or white based on the magnitude of the point
			double magnitude = 0;
			magnitude += IntersectPoint[0];
			magnitude = (magnitude * M_PI); // 
			double sineValue = sin(magnitude); 

			double color = 0; 
			if (sineValue >= 0) color = 1; // set to 1 half the wave 
			ray.col[0] = (color);// /(double)255 ;
			magnitude = IntersectPoint[1];
			magnitude = (magnitude * M_PI); // /
			 sineValue = sin(magnitude); 

			 color = 0; 
			if (sineValue >= 0) color = 1; // set to 1 ha
			ray.col[1] = (color);// /(double)255 ;
			magnitude = IntersectPoint[2];
			magnitude = (magnitude * M_PI); // 
			 sineValue = sin(magnitude); 

			 color = 0; 
			if (sineValue >= 0) color = 1; // set to 1 ha
			ray.col[2] = (color);// /(double)255 ;
		}
		else 
		{
			// Compute ray colours using Phong model
			ray.col[0] += mat_a[0] * I_a[0]; 
			ray.col[1] += mat_a[1] * I_a[1];
			ray.col[2] += mat_a[2] * I_a[2];					
			if(isInShadow == false && onlyAmbient==0 ) 
			{
				ray.col[0] += mat_d[0] * I_d[0] * fmax(0, s.dot(ray.intersection.normal))  + mat_s[0] * I_s[0] * pow(fmax(0, r.dot(b)), alpha);
				ray.col[1] += mat_d[1] * I_d[1] * fmax(0, s.dot(ray.intersection.normal))  + mat_s[1] * I_s[1] * pow(fmax(0, r.dot(b)), alpha);
				ray.col[2] += mat_d[2] * I_d[2] * fmax(0, s.dot(ray.intersection.normal))  + mat_s[2] * I_s[2] * pow(fmax(0, r.dot(b)), alpha);
			}
		}
		// Clamp values down to 1.0
		ray.col.clamp();
	}
}

