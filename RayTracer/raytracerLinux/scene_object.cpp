/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"

//----------------------------------------------------------------------------------------------------------------
bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel, const Matrix4x4& modelToWorld )
{
	// Implement intersection code for UnitSquare, which is
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
	// The original point is at 0,0,0 at its own model coordinates 

	// Get ray direction in model coordinates
	Ray3D ObjectRay = Ray3D(worldToModel*ray.origin, worldToModel * ray.dir);
	Point3D origin(0,0,0);
	// always assume unitsquare is on x,y plane 
	Vector3D squareNormal(0,0,1);
	double t_value;
	Point3D intersectionPoint;

	// If the ray is not moving towards the square, it means there is no intersection since
	// the square lies on the x and y axis, ray must move at some value on the z direction 
	if (ObjectRay.dir[2] != 0)
	{
		t_value  = -ObjectRay.origin[2]/ObjectRay.dir[2];
		// Check for no self intersection 
		if ((t_value )< 0.001)
		{
			return false;
		}
	}

	//the ray and the square are on the same plane
	else if (ObjectRay.dir[2] == 0 && ObjectRay.origin[2] == 0)
	{
		double temp1 = (-1.5 - ObjectRay.origin[0])/ObjectRay.dir[0];
		double temp2 = (1.5 - ObjectRay.origin[0])/ObjectRay.dir[0];
		double temp3 = (-1.5 - ObjectRay.origin[1])/ObjectRay.dir[1];
		double temp4 = (1.5 - ObjectRay.origin[1])/ObjectRay.dir[1];
		t_value  = temp1;
		if (temp2 < t_value )
		{
			t_value  = temp2;
		}
		if (temp3 < t_value )
		{
			t_value  =temp3;
		}
		if (temp4< t_value )
		{
			t_value  = temp4;
		}
	}
	else
	{

		return false;
	}
	
	// Compute t_value from z intersection at z = 0
	
	// If the ray already has an intersection at an earlier t_value, return but don't update
	// ray.intersection fields
	if(ray.intersection.none == false) 
	{
		if(t_value > ray.intersection.t_value) 
		{
			return false;
		}
	}

	// Compute coordinates
	Vector3D intersection(ObjectRay.origin[0] + ObjectRay.dir[0] * t_value, ObjectRay.origin[1] + ObjectRay.dir[1] * t_value, 0);
	double x = ObjectRay.origin[0] + ObjectRay.dir[0] * t_value;
	double y = ObjectRay.origin[1] + ObjectRay.dir[1] * t_value;

	// Test x
	if(x >= -0.5 && x <= 0.5) {
		// Text y
		if(y >= -0.5 && y <= 0.5) {
			// We have an intersection
			ray.intersection.t_value = t_value;
			ray.intersection.point = modelToWorld * Point3D(x, y, 0);
			ray.intersection.normal = Vector3D(0, 0, 1);

			// For Image Mapping 			
			ray.intersection.CenterPoint = modelToWorld * Point3D(0,0,0); // every model thinks it is in the center 
			ray.intersection.normal = transNorm(worldToModel, ray.intersection.normal);
			ray.intersection.normal.normalize();
			ray.intersection.none = false; // false means there is an intersection 

			return true;
		}
	}

	return false;
}

//----------------------------------------------------------------------------------------------------------------

bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel, const Matrix4x4& modelToWorld )
{
	// TODO: implement intersection code for UnitSphere, which is centred 
	// on the origin.  
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   

	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
	
	Ray3D ObjectRay = Ray3D(worldToModel * ray.origin, worldToModel * ray.dir);
	double a = pow(ObjectRay.dir[0], 2) + pow(ObjectRay.dir[1], 2) + pow(ObjectRay.dir[2], 2);
	double b = 2 * (ObjectRay.origin[0] * ObjectRay.dir[0] +
					ObjectRay.origin[1] * ObjectRay.dir[1] +
					ObjectRay.origin[2] * ObjectRay.dir[2]);

	double c = pow(ObjectRay.origin[0], 2) + 
			   pow(ObjectRay.origin[1], 2) + 
			   pow(ObjectRay.origin[2], 2) -1;

	double discriminant = pow(b, 2) - 4 * a * c;

	// No intersection since no real roots 
	if(discriminant < 0) 
	{ 
		return false;

	} 
	else 
	{
		double t_1 = (-1 * b + sqrt(discriminant)) / (2 * a);
		double t_2 = (-1 * b - sqrt(discriminant)) / (2 * a);
		double t_value = 0;

		if ( t_2 < 0 && t_1 < 0)
			return false;
		else if ( t_1 < 0 )
			t_value = t_2;
		else if ( t_2 < 0 )
			t_value = t_1;

		else
			t_value = fmin(t_1, t_2);

		// If the ray already has an intersection at an earlier t_value, return but don't update
		// ray.intersection fields

		if(ray.intersection.none == false) 
		{
			if(t_value > ray.intersection.t_value) 
			{
				return true;
			}
		}


		double x = ObjectRay.origin[0] + ObjectRay.dir[0]*t_value;
		double y = ObjectRay.origin[1] + ObjectRay.dir[1]*t_value;
		double z = ObjectRay.origin[2] + ObjectRay.dir[2]*t_value;
		
		Vector3D intersection = Vector3D(x, y, z);
		ray.intersection.t_value = t_value;		
		ray.intersection.point = modelToWorld * Point3D(x, y, z);
		ray.intersection.normal = transNorm(worldToModel, intersection);
		ray.intersection.normal.normalize();	
		// FOR TEXTURE MAPPING!
		ray.intersection.CenterPoint = modelToWorld * Point3D(0,0,0); // every model thinks it is in the center 
		ray.intersection.none = false;
		return true;
	}
}

//----------------------------------------------------------------------------------------------------------------

// Pseudocode from 
// Unit cylinder is whereby top and base are unit circle centered at (0,0,z), 
// z is between -0.5 to 0.5 => Height of cylinder is 1 
// Cylinder has equation x^2 + z^2 - r^2 = 0
// substituting line p + td
// Let q = intersectionPoint
// where t = (q-p)/d 
// Cylinder has equation: 
// (q-p-((q-p)/d)d)^2 - r^2 = 0 (Implicit form) 

// Basically, this is an intersection with 3 different connected surfaces 
// The top disk, the bottom disk and the cylinder body. 
// Test the intersection of all 3 bodies and figure out which has the lowest t value
bool UnitCylinder::intersect( Ray3D& ray, const Matrix4x4& worldToModel, const Matrix4x4& modelToWorld ) 
{
	// Get the ray in the object coordinates
	Ray3D ObjectRay = Ray3D(worldToModel * ray.origin, worldToModel * ray.dir);
	// A reference to the origin point in object/model coordinates
	Point3D centerOfCylinder(0,0,0); // Center of cylinder 
	// T Values 
	// One to check for intersection with disks, one to check for intersection with sides 
	double t_value; double t_valueTwo; 
	double tempOne, tempTwo;
	double radius = 1.0; // the radius of the disk of the cylinder

	Vector3D normal; // The normal of the cylinder Disk as well as wall of cylinder
			// Gets reused
	Point3D intersectionPoint; // intersection point with the cylinder 

//---------------------------------------------------------------------------------
// Checking Intersection with Cylinder Disks 
//---------------------------------------------------------------------------------

	// Note: Will only intersect the disks of the cylinder
	// if the direction value for z is not 0
	// Getting the tValue for the closer base
	if (ObjectRay.dir[2] != 0)
	{
		// Tvalue to intersect bottom disk of cylinder (negative z axis)
		tempOne = (-0.5-ObjectRay.origin[2])/ObjectRay.dir[2];
		// Tvalue to intersect top disk of cylinder (positive z axis)
		tempTwo = (0.5-ObjectRay.origin[2])/ObjectRay.dir[2];
		
		// Update the normal depending on which disk is closer 
		if (tempOne < tempTwo)
		{
			t_value = tempOne;
			// Construct the normal for the closer disk, which points on the negative z axis
			Point3D normal_temp(0,0,-1);
			normal = normal_temp - centerOfCylinder;
			normal.normalize();

		}
		else
		{
			t_value = tempTwo;
			Point3D normal_temp(0,0,1);
			normal = normal_temp - centerOfCylinder;
			normal.normalize();
		}
	}



//---------------------------------------------------------------------------------
//  Calculate the intersection point for the disk  
//---------------------------------------------------------------------------------
	intersectionPoint = ObjectRay.origin + t_value * ObjectRay.dir;
	
	// Prevent self intersection
	if (t_value < 0.001)
	{
		return false;
	}

	// Base or top of cylinder 
	// need make sure within the radius of the disk 
	if (intersectionPoint[0]*intersectionPoint[0] + intersectionPoint[1] * intersectionPoint[1] <= (radius*radius))
	{
		// If it intersection before and it's tvalue is lower then this, 
		// means there is no intersection, return false. 
		if (!ray.intersection.none && (t_value > ray.intersection.t_value))
		{
			return false;
		}
		// Else, there is intersection, update and return true 
		ray.intersection.point = intersectionPoint;
		
		ray.intersection.normal = normal;
		ray.intersection.t_value = t_value;
		ray.intersection.none = false;
		// Every model thinks it is in the center
		ray.intersection.CenterPoint = modelToWorld * Point3D(0,0,0); 
		return true;
	}

	// No intersection with cylinder disk
	// Check for intersection with cylinder body 
//---------------------------------------------------------------------------------
// Checking Intersection with Cylinder Body 
//---------------------------------------------------------------------------------

// Source: http://www.csie.ntu.edu.tw/~cyy/courses/rendering/pbrt-2.00/html/cylinder_8cpp_source.html

	// Compute Quadratic Cylinder Coefficients 
	double a = ObjectRay.dir[0]*ObjectRay.dir[0] + ObjectRay.dir[1]*ObjectRay.dir[1];
	double b = (ObjectRay.origin[0]*ObjectRay.dir[0] + ObjectRay.origin[1]* ObjectRay.dir[1]);
					
	double c = ObjectRay.origin[0]*ObjectRay.origin[0] + ObjectRay.origin[1]*ObjectRay.origin[1] - radius*radius;
													
	double discriminant = b*b-a*c;


	// If ray is aligned with body of the cylinder, there is no intersection 
	if (a == 0.0)
	{
		// return false
	}

	//If discrimant is < 0 => No real roots to the cylinder quadratic equation. 
	if (discriminant<0)
	{
		// It means no intersection 
		return false;
	}
	else
	{
		double rootOne = -b/a + sqrt(discriminant) / a;
		double rootTwo = -b/a - sqrt(discriminant) / a;
		if (rootOne< 0 && rootTwo < 0)
		{
			return false;
		}
		else if (rootOne> 0 && rootTwo < 0)
		{
			t_valueTwo = rootOne;
		}
		else
		{
			t_valueTwo = rootTwo;
		}
	}
//---------------------------------------------------------------------------------
//  Calculate the intersection point for the walls of the cylinder   
//---------------------------------------------------------------------------------
	// There is no intersection with the base of the cylinder on x,y plane
	// Check if there is intersection with sides. 

		// Get the intersection point for the sides 
		intersectionPoint = ObjectRay.origin+ t_valueTwo * ObjectRay.dir;
		// Check for no self intersection 
		if (t_valueTwo < 0.001)
		{
			return false;
		}
		// Get the coordinate of the normals for the sides of the cylinder 
		normal[0] = intersectionPoint[0];
		normal[1] = intersectionPoint[1];
		normal[2] = 0;
		normal.normalize();

		// If the intersection point is within the height of the cylinder
		if (intersectionPoint[2] < 0.5 && intersectionPoint[2] > -0.5)
		{
			// Ensure no previous intersection with lower t value 
			if (!ray.intersection.none && t_value > ray.intersection.t_value)
			{
				return false;
			}
			// Update the intersection to the ray 
			ray.intersection.point = modelToWorld * intersectionPoint;
			Point3D normalTemp;
			normalTemp[0] = intersectionPoint[0];
			normalTemp[1] = intersectionPoint[1];
			normalTemp[2] = 0;
			ray.intersection.normal = modelToWorld * (normalTemp - centerOfCylinder);
			ray.intersection.t_value = t_valueTwo;
			ray.intersection.none = false;
			// Every model thinks it is in the center
			ray.intersection.CenterPoint = modelToWorld * Point3D(0,0,0); 
			return true;
		}
		// If the intersection is not within the height of the cylinder,
		// It means there wasn't really a cylinder intersection to begin with  
		else
		{
			return false;
		}
}


// Assume disk is flat on the x,y plane with radius 1 
bool UnitDisk::intersect( Ray3D& ray, const Matrix4x4& worldToModel, const Matrix4x4& modelToWorld ) 
{
	// Get the ray in the object coordinates
	Ray3D ObjectRay = Ray3D(worldToModel * ray.origin, worldToModel * ray.dir);
	// A reference to the origin point in object/model coordinates
	Point3D centerofDisk(0,0,0); // Center of cylinder 
	// T Values 
	// One to check for intersection with disks, one to check for intersection with sides 
	double t_value; double t_valueTwo; 
	double tempOne, tempTwo;
	double radius = 1.0; // the radius of the disk of the cylinder

	Vector3D normal; // The normal of the cylinder Disk as well as wall of cylinder
			// Gets reused
	Point3D intersectionPoint; // intersection point with the cylinder 

//---------------------------------------------------------------------------------
// Checking Intersection with Cylinder Disks 
//---------------------------------------------------------------------------------

	// Note: Will only intersect the disks of the cylinder
	// if the direction value for z is not 0
	// Getting the tValue for the closer base
	if (ObjectRay.dir[2] != 0)
	{
		// To know if intersect from top or bottom of disk, and set normal accordingly 
		tempOne = (-0.0000001-ObjectRay.origin[2])/ObjectRay.dir[2];
		// Tvalue  
		tempTwo = (0.00000001-ObjectRay.origin[2])/ObjectRay.dir[2];
		
		// Update the normal depending on which disk is closer 
		if (tempOne < tempTwo)
		{
			t_value = tempOne;
			// Construct the normal for the closer disk, which points on the negative z axis
			Point3D normal_temp(0,0,-1);
			normal = normal_temp - centerofDisk;
			normal.normalize();

		}
		else
		{
			t_value = tempTwo;
			Point3D normal_temp(0,0,1);
			normal = normal_temp - centerofDisk;
			normal.normalize();
		}
	}



//---------------------------------------------------------------------------------
//  Calculate the intersection point for the disk  
//---------------------------------------------------------------------------------
	intersectionPoint = ObjectRay.origin + t_value * ObjectRay.dir;
	
	// Prevent self intersection
	if (t_value < 0.001)
	{
		return false;
	}

	// Base or top of cylinder 
	// need make sure within the radius of the disk 
	if (intersectionPoint[0]*intersectionPoint[0] + intersectionPoint[1] * intersectionPoint[1] <= (radius*radius))
	{
		// If it intersection before and it's tvalue is lower then this, 
		// means there is no intersection, return false. 
		if (!ray.intersection.none && (t_value > ray.intersection.t_value))
		{
			return false;
		}
		// Else, there is intersection, update and return true 
		ray.intersection.point = intersectionPoint;
		
		ray.intersection.normal = normal;
		ray.intersection.t_value = t_value;
		ray.intersection.none = false;
		// Every model thinks it is in the center
		ray.intersection.CenterPoint = modelToWorld * Point3D(0,0,0); 
		return true;
	}
	return false; // no intersection with disk 
}

