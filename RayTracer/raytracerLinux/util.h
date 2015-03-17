/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		utility functions and structures
		(based on code from CGL, University of Waterloo), 
		modify this file as you see fit.

***********************************************************/

// This file contains information about: 
//---------------------
// point3D
// Vector3D
// Vector4D
// Matrix 
//---------------------
// Colour
// Material
// Intersection
// Ray3D
//---------------------

#ifndef _UTIL_
#define _UTIL_

#include <iostream>
#include <cmath>

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

class Point3D {
private:
	// The 3 coordinates of the point (x,y,z) 
	double m_data[3];

public:
	// Constructors 
	Point3D(); 
	Point3D(double x, double y, double z);  
	Point3D(const Point3D& other); 
	// Overloading (=) operator 
	Point3D& operator =(const Point3D& other); 
	// Overloading ([]) operator with different return values
	double& operator[](int i); 
	double operator[](int i) const; 
};

// This class represents a Vector in 3D
class Vector3D {
private:
	// The 3 coordinates of the vector 
	double m_data[3];

public:
	// Constructors
	Vector3D(); 
	Vector3D(double x, double y, double z); 
	Vector3D(const Vector3D& other); 

	// Overloading Operators 
	Vector3D& operator =(const Vector3D& other); 
	double& operator[](int i);  
	double operator[](int i) const;  
	// Public Methods 
	
	// Find the length of the vector 
	double length() const; 
	// Normalize the vector 
	double normalize();
	// Do a dot product 
	double dot(const Vector3D& other) const; 
	// Do a cross product 
	Vector3D cross(const Vector3D& other) const; 


};
// Global Functions
// standard operators on points and vectors
// Pointwise operation on points and vectors 
Vector3D operator *(double s, const Vector3D& v); 
Vector3D operator +(const Vector3D& u, const Vector3D& v); 
Point3D operator +(const Point3D& u, const Vector3D& v); 
Vector3D operator -(const Point3D& u, const Point3D& v); 
Vector3D operator -(const Vector3D& u, const Vector3D& v); 
Vector3D operator -(const Vector3D& u); 
Point3D operator -(const Point3D& u, const Vector3D& v); 
Vector3D cross(const Vector3D& u, const Vector3D& v); 

// Output the points on the output stream 
std::ostream& operator <<(std::ostream& o, const Point3D& p); 
std::ostream& operator <<(std::ostream& o, const Vector3D& v); 

// This class represents a 4D Vector 
class Vector4D {
private:
	double m_data[4];

public:
	Vector4D(); 
	Vector4D(double w, double x, double y, double z); 
	Vector4D(const Vector4D& other); 

	Vector4D& operator =(const Vector4D& other); 
	double& operator[](int i);  
	double operator[](int i) const;  
};


// This is a 4by4 matrix 
// To set up transformation matrices 
class Matrix4x4 {
private:
  // (w,x,y,z) 
  double m_data[16];


public:
  Matrix4x4(); 
  Matrix4x4(const Matrix4x4& other); 
  Matrix4x4& operator=(const Matrix4x4& other); 
  // Methods 
  Vector4D getRow(int row) const; 
  double *getRow(int row); 
  Vector4D getColumn(int col) const; 

  Vector4D operator[](int row) const; 
  double *operator[](int row); 
  // Get the transpose of a Matrix 
  Matrix4x4 transpose() const; 
};

Matrix4x4 operator *(const Matrix4x4& M, const Matrix4x4& N); 
Vector3D operator *(const Matrix4x4& M, const Vector3D& v); 
Point3D operator *(const Matrix4x4& M, const Point3D& p);
// Multiply n by the transpose of M, useful for transforming normals.  
// Recall that normals should be transformed by the inverse transpose 
// of the matrix.  
Vector3D transNorm(const Matrix4x4& M, const Vector3D& n); 

// Overload the "<<" Operator to output properly onto the screen 
std::ostream& operator <<(std::ostream& os, const Matrix4x4& M); 

// This class represents the color 
class Colour {
private:
	// (r,g,b) 
	double m_data[3];
public:
	Colour(); 
	Colour(double r, double g, double b); 
	Colour(const Colour& other); 

	// Operator Overloading for (= , *, []) 
	Colour& operator =(const Colour& other); 
	Colour operator *(const Colour& other); 
	double& operator[](int i);  
	double operator[](int i) const; 
	// This function clamps the output 
	void clamp(); 	
};

// Pointwise multiplier 
Colour operator *(double s, const Colour& c); 
// Pointwise adder 
Colour operator +(const Colour& u, const Colour& v); 
std::ostream& operator <<(std::ostream& o, const Colour& c); 

// Material Structure for Phong Shading 
struct Material {
	// Ambient components for Phong shading.
	Colour ambient; 
	// Diffuse components for Phong shading.
	Colour diffuse;
	// Specular components for Phong shading.
	Colour specular;
	// Specular exponent.
	double specular_exp;

	// ImageMapping
	int imageMap; // 0 => no image mapping
	               // 1 =>image map to worldmap
			// 2 =>something else 

	// Constructor 
	Material( Colour ambient, Colour diffuse, Colour specular, double exp , int _imageMap) :
		ambient(ambient), diffuse(diffuse), specular(specular), specular_exp(exp), imageMap(_imageMap) {}
};

struct Intersection {
// Public 
	// Location of intersection (point)
	// Note:
	    // Position of the intersection point on your ray.
	    // (i.e. point = ray.origin + t_value * ray.dir)
	// Location of intersection.
	Point3D point;

	Point3D CenterPoint; // added center point of object for intersection 

	// Normal at the intersection.
	Vector3D normal;
	// Material at the intersection.
	Material* mat;
	// Position of the intersection point on your ray.
	// (i.e. point = ray.origin + t_value * ray.dir)
	// This is used when you need to intersect multiply objects and
	// only want to keep the nearest intersection.
	double t_value;	
	// Set to true when no intersection has occured.
	bool none; // true if no intersection at all
		   // false otherwise
};

// Ray structure. 
struct Ray3D {
	// Origin and direction of the ray.
	Point3D origin; // The point of origin of the ray Tracer
	Vector3D dir;	// the direction of the ray tracer 
	// Intersection status, should be computed by the intersection function.
	Intersection intersection;
	
	Matrix4x4 viewToWorldRay3D;

	// Current colour of the ray, should be computed by the shading function.
	Colour col;
	int depth; // number of depth 

	// Constructor without point or intersection 
	Ray3D() {
		intersection.none = true;  // Initialize with no intersection for the ray 
	}
	// Constructor with point and vector 
	Ray3D( Point3D p, Vector3D v ) : origin(p), dir(v) {
		intersection.none = true; // Initialize with no intersection for the ray 
	}

};
#endif





