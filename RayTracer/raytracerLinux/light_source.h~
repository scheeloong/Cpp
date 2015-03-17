/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		   light source classes

***********************************************************/
// This file contains information about: 
//---------------------
// LightSource
// PointLight (child of LightSource) 
//---------------------
#include "util.h"
#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif


// Base class for a light source.  You could define different types
// of lights here, but point light is sufficient for most scenes you
// might want to render.  Different light sources shade the ray differently.
class LightSource {
public:
	// Apply Shading 
	virtual void shade( Ray3D&, bool isInShadow, int onlyAmbient, int textureMap) = 0;
	// Get position 
	virtual Point3D get_position() const = 0; 
};

// A point light is defined by its position in world space and its colour.
class PointLight : public LightSource {
private:
	Point3D _pos;
	Colour _col_ambient;
	Colour _col_diffuse; 
	Colour _col_specular; 

public:
	// Constructor (positition and same color for all 3 components of Phong Shading)
	PointLight( Point3D pos, Colour col ) : _pos(pos), _col_ambient(col), 
	_col_diffuse(col), _col_specular(col) {}
	// Constructor, position and different colors for all 3 components of Phong Shading 
	PointLight( Point3D pos, Colour ambient, Colour diffuse, Colour specular ) 
	: _pos(pos), _col_ambient(ambient), _col_diffuse(diffuse), 
	_col_specular(specular) {}
	void shade( Ray3D& ray, bool isInShadow, int onlyAmbient, int textureMap );
	Point3D get_position() const { return _pos; }	
};
