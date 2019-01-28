/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		Implementations of functions in raytracer.h, 
		and the main function which specifies the 
		scene to be rendered.	

***********************************************************/
//The main body of the raytracer, including the scene graph.  Simple 
// addition and traversal code to the graph are provided to you. 

#include "raytracer.h"
#include "bmp_io.h"
#include <cmath> // for fmax() 
#include <iostream>
#include <cstdlib>
#include <time.h> // for random 

// Dear T.A. for CSC418 Fall 2014, 
// you can simply run each of our states by uncommenting below. 
// Except for Depth of Field which compiles for about 2 minutes, everything else should compile the first image within 10 seconds. 
//----------------------------------------------------------------------------------------------------------
// MODES for this ray tracer 
//----------------------------------------------------------------------------------------------------------
// Note: Uncomment if want them, must only uncomment one of these at a time 
#define SIGNATURE 1 // only ambient components 
//#define PARTONEFULL 1 // has antialias, without shadows 
//#define SHADOW 1 // normal shadow (compulsory)
//#define REFLECTION 1 // specular perfect reflection (compulsory) 
//#define SOFTSHADOW 1 // demonstrates soft shadows
//#define DEPTHOFFIELD 1 // demonstrates use of depth of field with 30 antialiasing (2 minute compile per pic)
//#define TEXTUREMAPPINGSPHERE 1 // texture image map sphere and texture procedural map plane 
//#define MOTIONBLUR 1
//#define CYLINDERDRAW 1 // includes both blue cylinder and purple disk 
//#define GLOSSYREFLECTION 1 
//#define FULLEVERYTHING 1 // takes few hours compile time 
//#define ANIMATION 1 // for rendering animation , takes few days to compile 

// to execute, type on terminal: 
// >> make clean
// >> make
// >> ./raytracer 
//----------------------------------------------------------------------------------------------------------
// Part 1
//----------------------------------------------------------------------------------------------------------
#ifdef SIGNATURE // only define if SIGNATURE is defined 
#define SCALESPHERE 1 // 1 means scale the sphere, 0 means don't 
#define ONLYAMBIENT 1 // only ambient 1 means no other phong shading 
#define NOSHADOW 1 // 1 means shadow is always off 
#define NUMLIGHTS 6.0 // must stay 6 unless add more or less points to check 
#define OFFSET 0.3   // a reasonable offset of light 
#define NUMANTIALIASE 1
#define RAYDEPTH 1// no reflection 
#define HEADLIGHTS 0 // use point lights
#define DODEPTHOFFIELD 0 //means no depth of field 
#define FOCALLENGTH 5 // focal length for depth of field 
#define DOTEXTUREMAPSPHERE 0 // 1 means do texture map on sphere, 0 means not. 
#define DOMOTIONBLUR 0 //  means do motion blur, 0 means do not. 
#define NUMGLOSSY 5 // number of times to reshoot glossy for reflection 
#define NUMMOTIONBLUR 5 // length of motion blur to move
#define DRAWCYLINDER 0 // 0 means don't draw cylinder
#define ONLYPARTONE 1 // only do part 1, so 2 shapes only 
#define DOANIMATION 0 // don't do animation 
#endif 
//----------------------------------------------------------------------------------------------------------
// Full of part 1 
//----------------------------------------------------------------------------------------------------------
#ifdef PARTONEFULL 
#define SCALESPHERE 1 // 1 means scale the sphere, 0 means don't 
#define ONLYAMBIENT 0 // only ambient 0 means phong shading works as usual 
#define NOSHADOW 1 // 1 means shadow is always off 
#define NUMLIGHTS 6.0 // must stay 6 unless add more or less points to check 
#define OFFSET 0.3   // a reasonable offset of light 
#define NUMANTIALIASE 5
#define RAYDEPTH 1// no reflection  
#define HEADLIGHTS 0 // use point lights
#define DODEPTHOFFIELD 0 //means no depth of field  
#define FOCALLENGTH 5 // focal length for depth of field 
#define DOTEXTUREMAPSPHERE 0 // 1 means do texture map on sphere, 0 means not. 
#define DOMOTIONBLUR 0 //  means do motion blur, 0 means do not. 
#define NUMGLOSSY 5 // number of times to reshoot glossy for reflection 
#define NUMMOTIONBLUR 5 // length of motion blur to move 
#define DRAWCYLINDER 0 // 0 means don't draw cylinder
#define ONLYPARTONE 1 // only do part 1, so 2 shapes only 
#define DOANIMATION 0 // don't do animation 
#endif 
//----------------------------------------------------------------------------------------------------------
// Part 2:Compulsory One: Shadow
//----------------------------------------------------------------------------------------------------------
#ifdef SHADOW 
#define SCALESPHERE 1 // 1 means scale the sphere, 0 means don't 
#define ONLYAMBIENT 0 // only ambient 0 means phong shading works as usual  
#define NOSHADOW 0 // 0 means shadow works as usual  
#define NUMLIGHTS 6.0 // must stay 6 unless add more or less points to check 
#define OFFSET 0.3   // a reasonable offset of light 
#define NUMANTIALIASE 1 // no antialias for soft shadow 
#define RAYDEPTH 1// no reflection 
// Set HeadLights to 0 to do shadows without soft shadows 
#define HEADLIGHTS 0 // dont use headlights for normal shadows
#define DODEPTHOFFIELD 0 //means no depth of field  
#define FOCALLENGTH 5 // focal length for depth of field 
#define DOTEXTUREMAPSPHERE 0 // 1 means do texture map on sphere, 0 means not. 
#define DOMOTIONBLUR 0 //  means do motion blur, 0 means do not. 
#define NUMGLOSSY 5 // number of times to reshoot glossy for reflection 
#define NUMMOTIONBLUR 5// length of motion blur to move 
#define DRAWCYLINDER 0 // 0 means don't draw cylinder
#define ONLYPARTONE 0 // don't only do part 1, so  many shapes 
#define DOANIMATION 0 // don't do animation 
#endif 
//----------------------------------------------------------------------------------------------------------
// Part 2:Compulsory Two: Reflection  
//----------------------------------------------------------------------------------------------------------
#ifdef REFLECTION
#define SCALESPHERE 0 // 1 means scale the sphere, 0 means don't 
#define ONLYAMBIENT 0 // only ambient 0 means phong shading works as usual  
#define NOSHADOW 1 // Turn off shadows  
#define NUMLIGHTS 6.0 // useless for now 
#define OFFSET 0.3   // offset not used cause no head lights 
#define NUMANTIALIASE 1 // no antialias 
#define RAYDEPTH 2// reflection occurs 
#define HEADLIGHTS 0 // use point lights (no soft shadows) 
#define DODEPTHOFFIELD 0 //means no depth of field  
#define FOCALLENGTH 5 // focal length for depth of field 
#define DOTEXTUREMAPSPHERE 0 // 1 means do texture map on sphere, 0 means not. 
#define DOMOTIONBLUR 0 // 1 means do motion blur, 0 means do not. 
// Set NUMGLOSSY to 1 for no glossy reflection (perfect reflection)  
#define NUMGLOSSY 1 // do perfect reflection 
#define NUMMOTIONBLUR 5 // length of motion blur to move 
#define DRAWCYLINDER 1 // 0 means don't draw cylinder
#define ONLYPARTONE 0 // don't only do part 1, so  many shapes 
#define DOANIMATION 0 // don't do animation 
#endif 
//----------------------------------------------------------------------------------------------------------
// Part 2: a) Soft Shadows 
//----------------------------------------------------------------------------------------------------------
#ifdef SOFTSHADOW 
#define SCALESPHERE 1 // 1 means scale the sphere, 0 means don't 
#define ONLYAMBIENT 0 // only ambient 0 means phong shading works as usual  
#define NOSHADOW 0 // 0 means shadow works as usual  
#define NUMLIGHTS 6.0 // must stay 6 unless add more or less points to check 
#define OFFSET 0.3   // a reasonable offset of light 
#define NUMANTIALIASE 1 // no antialias for soft shadow 
#define RAYDEPTH 1// no reflection 
// Set HeadLights to 0 to do shadows without soft shadows 
#define HEADLIGHTS 1 // use head lights (convert soft shadows to shadow) 
#define DODEPTHOFFIELD 0 //means no depth of field  
#define FOCALLENGTH 5 // focal length for depth of field 
#define DOTEXTUREMAPSPHERE 0 // 1 means do texture map on sphere, 0 means not. 
#define DOMOTIONBLUR 0 //  means do motion blur, 0 means do not. 
#define NUMGLOSSY 5 // number of times to reshoot glossy for reflection 
#define NUMMOTIONBLUR 5// length of motion blur to move 
#define DRAWCYLINDER 0 // 0 means don't draw cylinder
#define ONLYPARTONE 0 // don't only do part 1, so  many shapes 
#define DOANIMATION 0 // don't do animation 
#endif 
//----------------------------------------------------------------------------------------------------------
// Part 2: b) Depth of Field  (note: This takes long to run due to needed high number of aliasing. 
//----------------------------------------------------------------------------------------------------------
#ifdef DEPTHOFFIELD 
#define SCALESPHERE 1 // 1 means scale the sphere, 0 means don't 
#define ONLYAMBIENT 0 // only ambient 0 means phong shading works as usual  
#define NOSHADOW 0 // 0 means shadow works as usual  
#define NUMLIGHTS 6.0 // useless for now 
#define OFFSET 0.3   // offset not used cause no head lights 
#define NUMANTIALIASE 30// The antialias number also affects depth of field!!
#define RAYDEPTH 1// no reflection 
#define HEADLIGHTS 0 // use point lights (no soft shadows) 
#define DODEPTHOFFIELD 1 //means do depth of field 
#define FOCALLENGTH 6 // focal length for depth of field 
#define DOTEXTUREMAPSPHERE 0 // 1 means do texture map on sphere, 0 means not.
#define DOMOTIONBLUR 0 //  1 means do motion blur, 0 means do not.  
#define NUMGLOSSY 5 // number of times to reshoot glossy for reflection 
#define NUMMOTIONBLUR 5 // length of motion blur to move 
#define DRAWCYLINDER 0 // 0 means don't draw cylinder
#define ONLYPARTONE 0 // don't only do part 1, so  many shapes 
#define DOANIMATION 0 // don't do animation 
#endif 
//----------------------------------------------------------------------------------------------------------
// Part 2: c) Texture Mapping on sphere
//----------------------------------------------------------------------------------------------------------
#ifdef TEXTUREMAPPINGSPHERE
#define SCALESPHERE 0 // 1 means scale the sphere, 0 means don't 
#define ONLYAMBIENT 0 // only ambient 0 means phong shading works as usual  
#define NOSHADOW 1 // Turn off shadows  
#define NUMLIGHTS 6.0 // useless for now 
#define OFFSET 0.3   // offset not used cause no head lights 
#define NUMANTIALIASE 1 // no antialias 
#define RAYDEPTH 1// no reflection  
#define HEADLIGHTS 0 // use point lights (no soft shadows) 
#define DODEPTHOFFIELD 0 //means no depth of field  
#define FOCALLENGTH 5 // focal length for depth of field 
#define DOTEXTUREMAPSPHERE 1 // 1 means do texture map on sphere, 0 means not. 
#define DOMOTIONBLUR 0 // 1 means do motion blur, 0 means do not. 
#define NUMGLOSSY 5 // number of times to reshoot glossy for reflection 
#define NUMMOTIONBLUR 5 // length of motion blur to move 
#define DRAWCYLINDER 0 // 0 means don't draw cylinder
#define ONLYPARTONE 0 // don't only do part 1, so  many shapes 
#define DOANIMATION 0 // don't do animation 
#endif 
// Global buffers for texture map, these are global and used in light_source.cpp as well  
	unsigned char* rbuffer; 
	unsigned char* gbuffer; 
	unsigned char* bbuffer; 
	int width = 320; 
	int height = 240; 

//----------------------------------------------------------------------------------------------------------
// Part 2: d) Motion Blur 
//----------------------------------------------------------------------------------------------------------
#ifdef MOTIONBLUR
#define SCALESPHERE 0 // 1 means scale the sphere, 0 means don't 
#define ONLYAMBIENT 0 // only ambient 0 means phong shading works as usual  
#define NOSHADOW 1 // Turn off shadows  
#define NUMLIGHTS 6.0 // useless for now 
#define OFFSET 0.3   // offset not used cause no head lights 
#define NUMANTIALIASE 1 // no antialias 
#define RAYDEPTH 1// no reflection  
#define HEADLIGHTS 0 // use point lights (no soft shadows) 
#define DODEPTHOFFIELD 0 //means no depth of field  
#define FOCALLENGTH 5 // focal length for depth of field 
#define DOTEXTUREMAPSPHERE 0 // 1 means do texture map on sphere, 0 means not. 
#define DOMOTIONBLUR 1 // 1 means do motion blur, 0 means do not. 
#define NUMGLOSSY 1 // number of times to reshoot glossy for reflection 
#define NUMMOTIONBLUR 5 // length of motion blur to move 
#define DRAWCYLINDER 0 // 0 means don't draw cylinder
#define ONLYPARTONE 0 // don't only do part 1, so  many shapes 
#define DOANIMATION 0 // don't do animation 
#endif 
//----------------------------------------------------------------------------------------------------------
// Part 2: e) Cylinder
//----------------------------------------------------------------------------------------------------------
#ifdef CYLINDERDRAW
#define SCALESPHERE 0 // 1 means scale the sphere, 0 means don't 
#define ONLYAMBIENT 0 // only ambient 0 means phong shading works as usual  
#define NOSHADOW 1 // Turn off shadows  
#define NUMLIGHTS 6.0 // useless for now 
#define OFFSET 0.3   // offset not used cause no head lights 
#define NUMANTIALIASE 1 // no antialias 
#define RAYDEPTH 1// no reflection 
#define HEADLIGHTS 0 // use point lights (no soft shadows) 
#define DODEPTHOFFIELD 0 //means no depth of field  
#define FOCALLENGTH 5 // focal length for depth of field 
#define DOTEXTUREMAPSPHERE 0 // 1 means do texture map on sphere, 0 means not. 
#define DOMOTIONBLUR 0 // 1 means do motion blur, 0 means do not. 
// Set NUMGLOSSY to 1 for no glossy reflection (perfect reflection)  
#define NUMGLOSSY 1 // number of times to reshoot glossy for reflection 
#define NUMMOTIONBLUR 5 // length of motion blur to move 
#define DRAWCYLINDER 1 // 0 means don't draw cylinder
#define ONLYPARTONE 0 // don't only do part 1, so  many shapes 
#define DOANIMATION 0 // don't do animation 
#endif 
//----------------------------------------------------------------------------------------------------------
//  Part 2: f) Glossy Reflection (Additional 1 from 5 required)  
//----------------------------------------------------------------------------------------------------------
#ifdef GLOSSYREFLECTION
#define SCALESPHERE 0 // 1 means scale the sphere, 0 means don't 
#define ONLYAMBIENT 0 // only ambient 0 means phong shading works as usual  
#define NOSHADOW 1 // Turn off shadows  
#define NUMLIGHTS 6.0 // useless for now 
#define OFFSET 0.3   // offset not used cause no head lights 
#define NUMANTIALIASE 1 // no antialias 
#define RAYDEPTH 2// reflection occurs 
#define HEADLIGHTS 0 // use point lights (no soft shadows) 
#define DODEPTHOFFIELD 0 //means no depth of field  
#define FOCALLENGTH 5 // focal length for depth of field 
#define DOTEXTUREMAPSPHERE 0 // 1 means do texture map on sphere, 0 means not. 
#define DOMOTIONBLUR 0 // 1 means do motion blur, 0 means do not. 
// Set NUMGLOSSY to 1 for no glossy reflection (perfect reflection)  
#define NUMGLOSSY 5 // number of times to reshoot glossy for reflection 
#define NUMMOTIONBLUR 5 // length of motion blur to move 
#define DRAWCYLINDER 1 // 0 means don't draw cylinder
#define ONLYPARTONE 0 // don't only do part 1, so  many shapes 
#define DOANIMATION 0 // don't do animation 
#endif 
//----------------------------------------------------------------------------------------------------------
// Everything 
//----------------------------------------------------------------------------------------------------------
// This executes everything except:
// depth of field cause that takes too long
// texture map sphere cause wanna see reflections 
// glossy reflections cause wanna see perfect scene 
#ifdef FULLEVERYTHING
#define SCALESPHERE 1 // 1 means scale the sphere, 0 means don't 
#define ONLYAMBIENT 0 // only ambient 0 means phong shading works as usual  
#define NOSHADOW 0 // Turn on shadows  
#define NUMLIGHTS 6.0 // 
#define OFFSET 0.3   // offset not used cause no head lights 
#define NUMANTIALIASE 20 // anti alias  20 times
#define RAYDEPTH 2// reflection occurs 
#define HEADLIGHTS 1 // use head lights (has soft shadows) 
#define DODEPTHOFFIELD 0 //means with depth of field  
#define FOCALLENGTH 5 // focal length for depth of field 
#define DOTEXTUREMAPSPHERE 0 // 1 means do texture map on sphere, 0 means not. 
#define DOMOTIONBLUR 1 // 1 means do motion blur, 0 means do not. 
#define NUMGLOSSY 1 // number of times to reshoot glossy for reflection 
#define NUMMOTIONBLUR 5 // length of motion blur to move 
#define DRAWCYLINDER 1 // 0 means don't draw cylinder
#define ONLYPARTONE 0 // don't only do part 1, so  many shapes 
#define DOANIMATION 0 // don't do animation 
#endif 
//----------------------------------------------------------------------------------------------------------
// ANIMATION
//----------------------------------------------------------------------------------------------------------
// This executes everything except:
// depth of field cause that takes too long
// texture map sphere cause wanna see reflections 
// glossy reflections cause wanna see perfect scene 
#ifdef ANIMATION
#define SCALESPHERE 1 // 1 means scale the sphere, 0 means don't 
#define ONLYAMBIENT 0 // only ambient 0 means phong shading works as usual  
#define NOSHADOW 0 // Turn on shadows  
#define NUMLIGHTS 6.0 // 
#define OFFSET 0.3   // offset not used cause no head lights 
#define NUMANTIALIASE 20 // anti alias  20 times
#define RAYDEPTH 2// reflection occurs 
#define HEADLIGHTS 1 // use head lights (has soft shadows) 
#define DODEPTHOFFIELD 0 //means with depth of field  
#define FOCALLENGTH 5 // focal length for depth of field 
#define DOTEXTUREMAPSPHERE 0 // 1 means do texture map on sphere, 0 means not. 
#define DOMOTIONBLUR 1 // 1 means do motion blur, 0 means do not. 
#define NUMGLOSSY 1 // number of times to reshoot glossy for reflection 
#define NUMMOTIONBLUR 5 // length of motion blur to move 
#define DRAWCYLINDER 1 // 0 means don't draw cylinder
#define ONLYPARTONE 0 // don't only do part 1, so  many shapes 
#define DOANIMATION 1 // do animation (render images for animation) 
#endif 
//----------------------------------------------------------------------------------------------------------

// Note: 
// Refraction is difficult cause need another helper recursive function
// that detects if you enter an object, go inside in a different path
// and exit with the new ray calculated based on the new object entered 
// basically the whatever refraction law, easy to understand, complicated to implement
// need to handle backfacing etc. 

//----------------------------------------------------------------------------------------------------------
// For offset Shadowing
static bool headLightPresent = false; // this is useless, for now 
//----------------------------------------------------------------------------------------------------------
// Constructors
//----------------------------------------------------------------------------------------------------------

// Initialize the root with no light source
Raytracer::Raytracer() : _lightSource(NULL) {
	// The root is initialize with a new node where everything in this node is NULL 
	_root = new SceneDagNode();
}

// Delete the root only 
Raytracer::~Raytracer() {
	delete _root;
}

//----------------------------------------------------------------------------------------------------------
// Adding Scene Objects to the Directed Acyclic Graph 
//----------------------------------------------------------------------------------------------------------
// This is a method from RayTracer but returns a node instead
SceneDagNode* Raytracer::addObject(SceneDagNode* parent, SceneObject* obj, Material* mat) 
{
	// Once Raytracer is created, root points to a created object, thus, parent points to that object 

	// Create a new node with SceneDagNode for the Directed Acyclic Graph  with parent object and 
    // makes a new sceneDagNode with sceneObject obj, and material mat.
    // Every node in the graph therefore links the object with the material 
    // as well as the matrix transformation associated with it. 
	// Note: Created with parent pointing to NULL 
	SceneDagNode* node = new SceneDagNode(obj, mat); 
	// Make this node point to the parent 
	node->parent = parent;  // Note: Originally _root points to NULL, so first object added has parent point to NULL 
	// Initialize the next to NULL (basically sibling), transformation applied to this node not applied to child 
	node->next = NULL;
	// Initialize the child to NULL, transformation applied to this node will be applied to child
	node->child = NULL;
	// Add the object to the parent's child list, this means
	// whatever transformation applied to the parent will also
	// be applied to the child.
	if (parent->child == NULL) 
	{
		parent->child = node; // no sibling, add yourself as the only child
	}
	// go to the child of your parent and add it to the siblings of that child 
	else {
		SceneDagNode* sibling = parent->child; // has sibling, add yourself to one of the linked list of siblings
		while (sibling->next != NULL) 
		{
			sibling = sibling->next;
		}
		sibling->next = node;
	}
	return node;
}

//----------------------------------------------------------------------------------------------------------
// Adding Light objects to the linked list  
//----------------------------------------------------------------------------------------------------------
LightListNode* Raytracer::addLightSource(LightSource* light) {
	LightListNode* tempHead = _lightSource; // head of linked list
	// Add the light source to the front of linkedList of light source
	_lightSource = new LightListNode(light, tempHead); 
	return _lightSource; // return head of the linked list
}

LightListNode* Raytracer::addHeadLightSource(double posX, double posY, double posZ, double colorX, double colorY, double colorZ) {
headLightPresent = true; 
    // 6 positions for positive and negative offsets in X, Y, Z directions.
	this->addLightSource(new PointLight(Point3D(posX + OFFSET, posY, posZ), Colour(colorX/NUMLIGHTS, colorY/NUMLIGHTS, colorZ/NUMLIGHTS)));
	this->addLightSource(new PointLight(Point3D(posX - OFFSET, posY, posZ), Colour(colorX/NUMLIGHTS, colorY/NUMLIGHTS, colorZ/NUMLIGHTS)));
	this->addLightSource(new PointLight(Point3D(posX, posY + OFFSET, posZ), Colour(colorX/NUMLIGHTS, colorY/NUMLIGHTS, colorZ/NUMLIGHTS)));
	this->addLightSource(new PointLight(Point3D(posX, posY - OFFSET, posZ), Colour(colorX/NUMLIGHTS, colorY/NUMLIGHTS, colorZ/NUMLIGHTS)));
	this->addLightSource(new PointLight(Point3D(posX, posY , posZ + OFFSET), Colour(colorX/NUMLIGHTS, colorY/NUMLIGHTS, colorZ/NUMLIGHTS)));
	this->addLightSource(new PointLight(Point3D(posX, posY , posZ - OFFSET), Colour(colorX/NUMLIGHTS, colorY/NUMLIGHTS, colorZ/NUMLIGHTS)));
	return _lightSource; // return the head of the linked list 
}

//----------------------------------------------------------------------------------------------------------
// Transformation on each individual scene object given the SceneDagNode that contains that scene object
//----------------------------------------------------------------------------------------------------------
void Raytracer::rotate(SceneDagNode* node, char axis, double angle) {
	// Create a matrix for transformation 
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;
	for (i = 0; i < 2; i++) {
        // Create the rotation angle for the given axis if it is given
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
        // If axis isn't given, create both inverse and forward rotations.
		if (i == 0) {
            // Right multiply for rotation
		    node->trans = node->trans * rotation; 	
			angle = -angle;
		} 
		else {
            // Left multiply for the inverse 
			node->invtrans = rotation * node->invtrans; 
		}	
	}
}

void Raytracer::translate(SceneDagNode* node, Vector3D trans) {
	// A temporary matrix to store the information for translation and it's inverse
	// Node => The SceneDagNode that stores the SceneObject and material properties
		// as well as all the transformation it has to be multiplied into 
	Matrix4x4 translation;
	// Update the node's transformation to multiply with the translation  (world to object) 
	// transformation = transformation*A
	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	node->trans = node->trans * translation; 	
	// Update the inverse transformation to multiply with translation (object to world)
	// inverseTransformation = Ainverse*inverseTransformation 
	translation[0][3] = -trans[0];	// put a negative, which basically means it is the inverse translation 
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	node->invtrans = translation * node->invtrans; 
}

void Raytracer::scale(SceneDagNode* node, Point3D origin, double factor[3]) {
	// Create a matrix called scale to hold A 
	Matrix4x4 scale;
	
	// Create A 
	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];

	// transformation = transformation*A
	node->trans = node->trans*scale; 	
		
	// Create Ainverse
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	node->invtrans = scale*node->invtrans; 
	// inverseTransformation = Ainverse*inverseTransformation 
}

//----------------------------------------------------------------------------------------------------------
// ShadeRay() Helper Functions, note: ShadeRay() itself is a helper function for Render()
//----------------------------------------------------------------------------------------------------------

// Helper function called by shadeRay() 
void Raytracer::traverseScene(SceneDagNode* node, Ray3D& ray) 
{
	SceneDagNode *childPtr;
	// Applies transformation of the current node to the global
	// transformation matrices.
	_modelToWorld = _modelToWorld*node->trans;
	_worldToModel = node->invtrans*_worldToModel; 
	if (node->obj)
	{
		// Perform intersection.
		if (node->obj->intersect(ray, _worldToModel, _modelToWorld)) 
		{
			ray.intersection.mat = node->mat;
		}
	}
	// Traverse the children.
	childPtr = node->child;
	while (childPtr != NULL) {
		traverseScene(childPtr, ray);
		childPtr = childPtr->next;
	}

	// Removes transformation of the current node from the global
	// transformation matrices.
	_worldToModel = node->trans*_worldToModel;
	_modelToWorld = _modelToWorld*node->invtrans;
}

// Helper function called by shadeRay() 
void Raytracer::computeShading(Ray3D& ray) 
{
	LightListNode* curLight = _lightSource;
	
	// Loop through all lights 
	for (;;) 
    {
		if (curLight == NULL) break;
		// Each lightSource provides its own shading function.

		// SHADOWS
		bool isInShadow = false;
		
		Point3D lightPos = curLight->light->get_position();
		if (ray.intersection.none == false) 
		{
			// Construct a ray from intersection point on ray to light source
			Point3D intersectPos = ray.intersection.point;
			double dx = lightPos[0] - intersectPos[0];
			double dy = lightPos[1] - intersectPos[1];
			double dz = lightPos[2] - intersectPos[2];
			double d = 0.01;

			// Set slight offset from intersection
			intersectPos[0] = dx > 0 ? intersectPos[0] + d : intersectPos[0] - d;
			intersectPos[1] = dy > 0 ? intersectPos[1] + d : intersectPos[1] - d;
			intersectPos[2] = dz > 0 ? intersectPos[2] + d : intersectPos[2] - d;

			Vector3D vec = Vector3D(dx, dy, dz);
			double epsilon = 0.001; 
			Ray3D shadowRay = Ray3D(intersectPos, lightPos - intersectPos);
			//Ray3D shadowRay = Ray3D(lightPos, intersectPos-lightPos);
			traverseScene(_root, shadowRay);
		
			// Check if light intersects at intersectPos on object, if not then it's in shadow 
			if (shadowRay.intersection.none == false) 
			{
				//double t_light = (lightPos[0] - shadowRay.origin[0]) / shadowRay.dir[0];
				Point3D shadowRayIntersect = shadowRay.intersection.point;
				if (shadowRay.intersection.t_value >= epsilon &&   shadowRay.intersection.t_value<= 1-epsilon)
				{
					isInShadow = true; 
				}
			}
/* // This is only a hack for making everything super bright and shiny 
// can fix by generating every ray as if it is the camera position. 
				// TODO: CHANGE TO RAY INTERSECTION POINT MAYBE 
				Vector3D l = Vector3D(intersectPos-lightPos);
//				Vector3D l = Vector3D(lightPos - intersectPos);
				l.normalize();
				// Reflection of l across n, vector 
				Vector3D r = Vector3D(-1 * l + 2 * ray.intersection.normal.dot(l) * ray.intersection.normal);
				// Uncomment below to check later if changes anything 
				// update the ray Direction to change from the view matrix to the world matrix 
				r = ray.viewToWorldRay3D * r;
				// normalize the ray direction 
				r.normalize();

				// Now create the actual ray given the origin and ray direction (from the point origin given above) 
				Ray3D reflectRay = Ray3D(ray.intersection.point, r);
				reflectRay.depth = ray.depth; // go to a reflectRay Depth of current ray depth 
				reflectRay.viewToWorldRay3D = ray.viewToWorldRay3D; 
				// Get the color for this pixel by calling the shade ray function 
				// shadeRay(ray) to generate pixel colour. 	
				Colour colTwo = shadeRay(reflectRay); 

				// Update colour of current ray (take into account lose of energy each time) 
				double numLightsTemp = 1; 
				if (headLightPresent) numLightsTemp =  NUMLIGHTS; 
				ray.col[0] += colTwo[0]/((RAYDEPTH - reflectRay.depth)*numLightsTemp); // need divide by remaining depth of recursion each time to lose energy 
				ray.col[1] += colTwo[1]/((RAYDEPTH - reflectRay.depth)*numLightsTemp); 
				ray.col[2] += colTwo[2]/((RAYDEPTH - reflectRay.depth)*numLightsTemp); 
				ray.col.clamp();
*/
		}
		if (NOSHADOW == 1) isInShadow = false; // always make it not in shadow 
		curLight->light->shade(ray, isInShadow, ONLYAMBIENT, DOTEXTUREMAPSPHERE);
		curLight = curLight->next;
	}
	//Done looping through all lights, now compute reflection 

// This is correct! You just need to only reflect if intersection point is on reflective objects 
// The r vector is computed correctly
	// Intersection point to point light source position
	Vector3D l = Vector3D(ray.origin - ray.intersection.point);
	l.normalize();
	Colour colTwo;
	colTwo[0] = 0;	colTwo[1] = 0;	colTwo[2] = 0;
	Ray3D reflectRay;

	for(int lmn = 0; lmn < NUMGLOSSY; lmn++)
	{
		// Reflection of l across n, vector 

		Vector3D r = Vector3D(-1 * l + 2 * ray.intersection.normal.dot(l) * ray.intersection.normal);
    		double rng1 = rand()/double(RAND_MAX);//generating random number
    		double rng2 = rand()/double(RAND_MAX);//generating random number
    		double rng3 = rand()/double(RAND_MAX);//generating random number
		if ((rand() % 10) > 4) rng1*= -1; // deduct instead 
		if ((rand() % 10) > 4) rng2*= -1; // deduct instead 
		if ((rand() % 10) > 4) rng3*= -1; // deduct instead 
		// only do glossy if NUMGLOSSY is more than 1, else do perfect 
		if (NUMGLOSSY > 1)
		{
			// Add a small random offset to r 
			r[0] += rng1; 	
			r[1] += rng2; 	
			r[2] += rng3; 	
		}
	// Uncomment below to check later if changes anything 
	// update the ray Direction to change from the view matrix to the world matrix 
		// normalize the ray direction 
		r.normalize();
		Point3D tempPoint = ray.intersection.point; 
		tempPoint[0] = 1000*r[0]; 
		tempPoint[1] = 1000*r[1]; 
		tempPoint[2] = 1000*r[2]; 

		// Now create the actual ray given the origin and ray direction (from the point origin given above) 
		 reflectRay = Ray3D(ray.intersection.point, r);
		//Ray3D reflectRay = Ray3D(tempPoint, -r); // note: This causes a recursive map into the object 
							// if intersection implemented correctly. 
		reflectRay.depth = ray.depth; // go to a reflectRay Depth of current ray depth 
		reflectRay.viewToWorldRay3D = ray.viewToWorldRay3D; // useless for now 
		// Get the color for this pixel by calling the shade ray function 
		// shadeRay(ray) to generate pixel colour. 	
		Colour tempCol = shadeRay(reflectRay); // get the colour 
		colTwo[0] += tempCol[0]; colTwo[1] += tempCol[1]; colTwo[2] += tempCol[2];  

		// Update colour of current ray (take into account lose of energy each time) 
		
	}

	ray.col[0] += (colTwo[0]/(double)NUMGLOSSY)/(double)(RAYDEPTH - reflectRay.depth); // need divide by remaining depth of recursion each time to lose energy 
	ray.col[1] += (colTwo[1]/(double)NUMGLOSSY)/(double)(RAYDEPTH - reflectRay.depth); 
	ray.col[2] += (colTwo[2]/(double)NUMGLOSSY)/(double)(RAYDEPTH - reflectRay.depth); 

	// pg 103/785 in pdf 
	// Option 2: Multiply with specular reflection to lose energy instead of the reflect ray 
/*
	if (reflectRay.intersection.none == false) 
	{
		Vector3D l2 = Vector3D(reflectRay.origin - reflectRay.intersection.point);
 //l2 = l; // TEMP
		l2.normalize();
		Vector3D r2 = Vector3D(-1 * l2 + 2 * reflectRay.intersection.normal.dot(l2) * reflectRay.intersection.normal);
 //r2 = r; // TEMP 
		r2.normalize();
			double alpha = reflectRay.intersection.mat->specular_exp;
			// Reflection of s across n
			// Intersection point to camera position
		Colour mat_s =  Colour(reflectRay.intersection.mat->specular);
//		Colour mat_d =  Colour(ray.intersection.mat->diffuse);
		ray.col[0] += colTwo[0]*mat_s[0]* pow(fmax(0, r2.dot(l2)), alpha); // need divide by remaining depth of recursion each time to lose energy 
		ray.col[1] += colTwo[1]*mat_s[1]*pow(fmax(0, r2.dot(l2)), alpha); 
		ray.col[2] += colTwo[2]*mat_s[2]* pow(fmax(0, r2.dot(l2)), alpha); 
	}
*/
	ray.col.clamp();

}

//----------------------------------------------------------------------------------------------------------
// Render() Helper Functions 
//----------------------------------------------------------------------------------------------------------

// This initializes the view matrix based on the position, direction and orientation of the eye/camera . 
Matrix4x4 Raytracer::initInvViewMatrix(Point3D eye, Vector3D view, Vector3D up) 
{
	// Create a 4by4 matrix 
	Matrix4x4 mat; 

	// Create a 3D Vector, w 
	Vector3D w;
	
	// Normalize the view/direction Vector given 
	view.normalize();
	
	// Question: What does this do? 
	// Update the orientation vector to be 
	// orientation = orientation - (orientation(dot)view)*view
	up = up - up.dot(view)*view; 
			// (make the view vector the length of the dot product)
			// make the up vector deduct the view vector 
	
	// normalize the orientation vector 
	up.normalize();
	
	// Set up the 3D vector created to be the normal of the view vector and the up vector 
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
// mat => resulting matrix
// 	wx	upx	-viewx	 eyex
//	wy	upy	-viewy	 eyey
// 	wz	upz	-viewz	 eyez
	return mat; 
}

// Helper function called by render()
// This function initializes the pixel buffer assuming _scrWidth & _scrheight has been given the right values 
void Raytracer::initPixelBuffer() 
{
	int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
	_rbuffer = new unsigned char[numbytes];
	_gbuffer = new unsigned char[numbytes];
	_bbuffer = new unsigned char[numbytes];
	for (int i = 0; i < _scrHeight; i++) 
	{
		for (int j = 0; j < _scrWidth; j++) 
		{
			_rbuffer[i*_scrWidth+j] = 0;
			_gbuffer[i*_scrWidth+j] = 0;
			_bbuffer[i*_scrWidth+j] = 0;
		}
	}
}

// Helper Function called by render() 
// This flushes the pixel buffer by writing it onto the bmp file, and then deleting all the buffers for thsi ray tracer 
void Raytracer::flushPixelBuffer(char *file_name) 
{
	// Note: Change this function if you want to write to an image with a different format or something else (VGA Screen, etc.) 
	bmp_write(file_name, _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer);
	delete _rbuffer;
	delete _gbuffer;
	delete _bbuffer;
}

// This function computes a color for a particular pixel, given a ray towards that pixel. 
Colour Raytracer::shadeRay(Ray3D& ray) 
{
	// Initialize the color to black 
	Colour col(0.0, 0.0, 0.0); 
	if (ray.depth == 0)
		return col; 
	ray.depth--; 
	// Traverse from the root 
	traverseScene(_root, ray); 
	
	// Only shade if ray hit something
	if (!ray.intersection.none) 
	{
		computeShading(ray); 
		col = ray.col;  
	}
	// You'll want to call shadeRay recursively (with a different ray, 
	// of course) here to implement reflection/refraction effects. (call shade ray in computeShading itself)  
	return col; 
}	

//----------------------------------------------------------------------------------------------------------
// Render()
//----------------------------------------------------------------------------------------------------------

// This function renders the entire scene once it is done being set up. 
void Raytracer::render(int width, int height, Point3D eye, Vector3D view, Vector3D up, double fov, char* fileName, SceneDagNode* node) 
{
	// Create a matrix called view to world 
	Matrix4x4 viewToWorld;
	// Set the screen width and the screen height of the raytracer to width and height given 
	_scrWidth = width;
	_scrHeight = height;
	
	// Calculate the factor 
	// QUESTION: WHAT IS THIS????? 
	double factor = (double(height)/2)/tan(fov*M_PI/360.0); // used to determine where the rays get projected 

	// Initialize the pixel buffer after setting _scrWidth and _scrHeight 
	initPixelBuffer();

	// Initialize the viewToWorld Matrix given the camera properties 
	viewToWorld = initInvViewMatrix(eye, view, up);
	
	int numMotion = 0; 
	int limit = 0; 
	if (DOMOTIONBLUR == 1) 
	{
		numMotion = NUMMOTIONBLUR; // execute normally if no motion blur 
		limit = 1; 
	}	
	// 20 times for motion blur 
	for(int lkj = numMotion; lkj >=limit ; lkj--)
	{
		// Construct a ray for each pixel.
		for (int i = 0; i < _scrHeight; i++) 
		{

			for (int j = 0; j < _scrWidth; j++) 
			{
				// Sets up ray origin and direction in view space, image plane is at z = -1.
				Point3D origin(0, 0, 0); // Create a point called origin 
				Point3D imagePlane;	// the image plane where it is some distance away from the camera 
	
							// It initializes the ray direction 
		// Original Using Center 
				imagePlane[0] = (-double(width)/2 + 0.5 + j)/factor; // x position of image plane (factor influences the final width of the image plane) 
				imagePlane[1] = (-double(height)/2 + 0.5 + i)/factor; // y position of image plane 
				Colour col; // initialize col 
				col[0] = 0; col[1] = 0;	col[2] = 0; 
				srand(time(NULL));  
			

				for(int ij =0; ij < NUMANTIALIASE; ij++)
				{
					// Anti-aliasing
	/*
	// ORIGINAL code given , looking from positive Z_axis 
					imagePlane[0] = (-double(width)/2 + ((double)rand()/(double)RAND_MAX) + j)/factor;
					imagePlane[1] = (-double(height)/2 + ((double)rand()/(double)RAND_MAX) + i)/factor;
					// The smaller this value, the closer the object gets 
					imagePlane[2] = -1; // z position of image plane 
	*/
	/*
	// Looking from negative X_axis 
	// With camera position: 
	//	Point3D eye3(-1, 0, 0);
	//	Vector3D view3(1, 0, 0);
					imagePlane[1] = (-double(width)/2 + ((double)rand()/(double)RAND_MAX) + j)/factor;
					imagePlane[2] = (-double(height)/2 + ((double)rand()/(double)RAND_MAX) + i)/factor;
					// The smaller this value, the closer the object gets 
					imagePlane[0] = -1; // x position of image plane 
	*/		
	
					// Note: These 3 lines are corrected, you assume it is in the view's coordinates, then you translate into actual world coordinates 
					imagePlane[0] = (-double(width)/2 + ((double)rand()/(double)RAND_MAX) + j)/factor;
					imagePlane[1] = (-double(height)/2 + ((double)rand()/(double)RAND_MAX) + i)/factor;
					// The smaller this value, the closer the object gets 
					imagePlane[2] = -1; // z position of image plane 

					// Construct the ray direction based on the points on the image plane computed 
					Vector3D rayDir = Vector3D(imagePlane[0], imagePlane[1], imagePlane[2]);
					// point Aimed is the position of the focal plane in specified direction 
					Point3D pointAimed; // for depth of field (depends on focal length) 
					Point3D start = eye; // initialize
				// If need to do depth of field 
				if (DODEPTHOFFIELD == 1) 
				{
					rayDir = viewToWorld * rayDir;
					// normalize the ray direction 
					rayDir.normalize();
					// start further instead of the usual eye position 
					pointAimed = eye + FOCALLENGTH*rayDir; 
					double radius = 1; // radius of the lens 

	    				double du = rand()/double(RAND_MAX+1);//generating random number
	    				double dv = rand()/double(RAND_MAX+1);
					// Need jitter position of start (important for depth of field) 
					// creating new camera position(or ray start using jittering)
					// Definition of u and v might be wrong,
					// may need to be direction orthogonal to light direction. 
					//double u = radius * cos(fov); 
					//double v = radius * sin(fov); 
					Vector3D u = up; 
					Vector3D v = up.cross(view); 
					start[0]= eye[0] -(radius/2)*u[0]-(radius/2)*v[0]+radius*(du)*u[0]+radius*(dv)*v[0]; 
					start[1]= eye[1] -(radius/2)*u[1]-(radius/2)*v[1]+radius*(du)*u[1]+radius*(dv)*v[1]; 
					start[2]= eye[2] -(radius/2)*u[2]-(radius/2)*v[2]+radius*(du)*u[2]+radius*(dv)*v[2]; 
					 //getting the new direction of ray
					 rayDir = pointAimed - start;
					 Ray3D ray = Ray3D(eye, rayDir);
					 // normalize the ray direction 
					 rayDir.normalize();
				
				}
			
				else
				{
					// update the ray Direction to change from the view matrix to the world matrix 
					rayDir = viewToWorld * rayDir;
					// normalize the ray direction 
					rayDir.normalize();
				}
					Ray3D ray = Ray3D(start, rayDir);
					// Now create the actual ray given the origin and ray direction (from the point origin given above) 
					// THIS IS DEFINITELY WRONG AS IT IS ONLY
					// INFLUENCED BY THE EYE DIRECTION AND NOT EYE POSITION 
					//Ray3D ray = Ray3D(origin, rayDir); // original code : THIS IS WRONG!! YOU NEED START FROM EYE, NOT ORIGIN 

					ray.depth = RAYDEPTH; // go to a ray depth of RAYDEPTH 
					ray.viewToWorldRay3D = viewToWorld; 

					// Get the color for this pixel by calling the shade ray function 
					// shadeRay(ray) to generate pixel colour. 	
					Colour colTwo = shadeRay(ray); 
					col[0] += colTwo[0]; 
					col[1] += colTwo[1]; 
					col[2] += colTwo[2]; 
				}

				col[0] /= (double) NUMANTIALIASE;  
				col[1] /= (double) NUMANTIALIASE; 
				col[2] /= (double) NUMANTIALIASE;  
				// Finally, color the pixel based on the color that was given 
				// Add less and less each sweep 
				_rbuffer[i*width+j] += int(col[0]*255)/(double)pow(2.0,lkj);
				_gbuffer[i*width+j] += int(col[1]*255)/(double)pow(2.0,lkj);
				_bbuffer[i*width+j] += int(col[2]*255)/(double)pow(2.0,lkj);;
			}
		}
		if (DOMOTIONBLUR == 1)
		{
			// move the sphere
			this->translate(node, Vector3D(0, 0.5 , 0));	
		}
	} 		// MOTION BLUR 
	if (DOMOTIONBLUR == 1)
	{
		// Reset the sphere to original position. 
		this->translate(node, Vector3D(0, -0.5*numMotion , 0));	
	}
	// Finally, flush the entire pixelbuffer to the file name given. 
	flushPixelBuffer(fileName);
}

//----------------------------------------------------------------------------------------------------------
// Main Function () 
//----------------------------------------------------------------------------------------------------------
int main(int argc, char* argv[])
{	

	// Build your scene and setup your camera here, by calling 
	// functions from Raytracer.  The code here sets up an example
	// scene and renders it from two different view points, DO NOT
	// change this if you're just implementing part one of the 
	// assignment.  
	//-----------------------------------------------------------------------------------------------------------
	// Step1 : Initialize the main raytracer, screen width and height and camera 
	//-----------------------------------------------------------------------------------------------------------
	Raytracer raytracer;
		//-------------------------------------------------
		// contains: int 		_scrWidth
		// 	     int 		_scrHeight
		// 	     LightListNode* 	_lightSource 	// a linked list of light sources
		//	     SceneDagNode*	_root		// The root of the linked list graph data structure 
		// 	     unsigned char*	_rbuffer;	
		// 	     unsigned char*	_gbuffer;	
		// 	     unsigned char*	_bbuffer;	
		//	     Matrix4x4		_modelToWorld
		//	     Matrix4x4		_worldToModel
		//-------------------------------------------------
	// Set default width and height of the pixel screen 
	// Update width and height if given in commandline 
	if (argc == 3) 
	{
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}

    //-------------------------------------------------------------------------------------
    if (DOTEXTUREMAPSPHERE == 1)
    {
        // For texture map, the buffers are defined globally  
        // Initialize the buffers 
        int numbytesTextureMap = width * height * sizeof(unsigned char);
        rbuffer = new unsigned char[numbytesTextureMap];
        gbuffer = new unsigned char[numbytesTextureMap];
        bbuffer = new unsigned char[numbytesTextureMap];
        for (int i = 0; i < height; i++) 
        {
            for (int j = 0; j < width; j++) 
            {
                rbuffer[i*width+j] = 0;
                gbuffer[i*width+j] = 0;
                bbuffer[i*width+j] = 0;
            }
        }
        int temp1 = width;
        int temp2 = height; 
        long unsigned int* widthTextureMap = new long unsigned int(width);
        long int* heightTexturemap = new long int(height); 
        // Read the texture map input image
        bmp_read("worldMap.bmp", widthTextureMap, heightTexturemap, &rbuffer, &gbuffer, &bbuffer); 
        std::cout <<" Finished reading: " << "worldMap.bmp" << std::endl; 
        unsigned char* _rbuffer;
        unsigned char* _gbuffer;    				
        unsigned char* _bbuffer;
        // test to see if data is correct
        _rbuffer = new unsigned char[numbytesTextureMap];
        _gbuffer = new unsigned char[numbytesTextureMap];
        _bbuffer = new unsigned char[numbytesTextureMap];
        for (int i = 0; i < height; i++) 
        {
            for (int j = 0; j < width; j++) 
            {
                _rbuffer[i*width+j] = (rbuffer[i*width+j]) ;
                _gbuffer[i*width+j] = (gbuffer[i*width+j]);
                _bbuffer[i*width+j] = (bbuffer[i*width+j]);
            }
        }
        // Write the texture map output buffer 
        // note: No algorithm done here, just to ensure reading and writing works
        bmp_write("output.bmp", width, height, _rbuffer, _gbuffer, _bbuffer);
        delete _rbuffer;
        delete _gbuffer;
        delete _bbuffer;
    }// end of if statement for texture map

    // Here, successfully read texture map into:
    // rbuffer
    // gbuffer
    // bbuffer 

    //-------------------------------------------------------------------------------------
    // Continue program as usual 
	// Camera parameters.
	Point3D eye(0, 0, 1);
	Vector3D view(0, 0, -1);
	Vector3D up(0, 1, 0);

	// Angle for field of view is 60 degrees
	double fov = 60;
	//-----------------------------------------------------------------------------------------------------------
	// Step2 : Add Light Sources and Scene Objects into the scene 
	//-----------------------------------------------------------------------------------------------------------
	// Create materials with AmbientColor, DiffuseColor, SpecularcColor, and SpecularExponent 
	// Note: Material is used for shading 
	// Material 1: Gold
	Material gold(Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648), Colour(0.628281, 0.555802, 0.366065), 51.2, DOTEXTUREMAPSPHERE); // 1 => sphere does texture map 
	//Material gold(Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648), Colour(0.628281, 0.555802, 0.366065), 1); // Test: Material with lower exponent 
	// Material 2: Jade
	int tempJade = 0; 
	if (DOTEXTUREMAPSPHERE == 1) tempJade = 2; 
	Material jade(Colour(0, 0, 0), Colour(0.54, 0.89, 0.63), Colour(0.316228, 0.316228, 0.316228), 12.8 , tempJade); // 0 => No texture map 
	//	Material jade(Colour(0, 0, 0), Colour(0.54, 0.89, 0.63), Colour(0.316228, 0.316228, 0.316228), 1);	
	Material mud(Colour(0.6, 0.0, 0.0), Colour(0.75164, 0.60648, 0.22648), Colour(0.628281, 0.555802, 0.366065), 32.2, 0);

	// Material for the cylinder 
	Material cylMat(Colour(0.0, 0.0, 0.6), Colour(0.75164, 0.60648, 0.22648), Colour(0.628281, 0.555802, 0.366065), 32.2, 0);

	// Material for the disk
	Material diskMat(Colour(0.5, 0.0, 0.5), Colour(0.75164, 0.60648, 0.22648), Colour(0.628281, 0.555802, 0.366065), 32.2, 0);

	// NOTE: Can only use HeadLight or Light Source but not both! 
	if (HEADLIGHTS==1)
	{
		// Note: Below takes way longer with reflection implemented 
	//	raytracer.addHeadLightSource(-5, 0, 0, 0.9, 0.9, 0.9); 
		raytracer.addHeadLightSource(5, 0, 0, 0.9, 0.9, 0.9);
	}

	else
	{
		// Defines a point light source with location in 3D as well as the color of the light. (close to white) 
		// Add the point light source to the ray tracer's list of lightSource 
	    // raytracer.addLightSource(new PointLight(Point3D(-5, 0, 0), Colour(0.9, 0.9, 0.9)));
		raytracer.addLightSource(new PointLight(Point3D(5, 0, 0), Colour(0.9, 0.9, 0.9)));
	}
	// raytracer.addLightSource(new PointLight(Point3D(5, 0, 0), Colour(0.4, 0.4, 0.4))); // Test: Darker light source 
	// Need associate created material with objects 
	// Adding actual objects to world  (known as SceneObject) 
	// Add a unit square into the scene with material mat
	// Note: It is pointer that points to objects created and returned by raytracer 
						// (SceneObject, Material) 
	SceneDagNode* sphere = raytracer.addObject(new UnitSphere(), &gold);
	SceneDagNode* plane = raytracer.addObject(new UnitSquare(), &jade);

	// now root has 2 childs, both the unitSphere and the unitSquare 
	// but really it is
	// root->child->Sphere
	// root->Sphere->next->plane
    
    // If don't do animation 
    if (DOANIMATION == 0) 
    {	
        //-----------------------------------------------------------------------------------------------------------
        // Step3 : Apply transformation to the objects to put them in place 
        //-----------------------------------------------------------------------------------------------------------
        // Apply some transformations to the unit square.

        // Create a factor for scaling (ScaleX, ScaleY, ScaleZ) 
        double factor1[3] = { 1.0, 2.0, 1.0 };
        double factor2[3] = { 6.0, 6.0, 6.0 };
        double factor3[3] = { 1.0, 1.0, 1.0 };
        double factor4[3] = { 0.25, 0.25, 0.25 };
        // NEW QUADRIC SURFACE 
        if (DRAWCYLINDER == 1)
        {
            SceneDagNode* cylinder = raytracer.addObject(new UnitCylinder(), &cylMat);
            raytracer.translate(cylinder, Vector3D(-2,-2,-5));
            raytracer.scale(cylinder, Point3D(0,0,0), factor3);

            // only draw disk if no glossy 		
            if (NUMGLOSSY <= 1)
            {
                SceneDagNode* disk = raytracer.addObject(new UnitDisk(), &diskMat);
                raytracer.translate(disk, Vector3D(-1,1,-2));
                raytracer.scale(disk, Point3D(0,0,0), factor4);
            }
        }
        // Do transformations on sphere 
        // Add the translation matrix to the sphere's list of transformation and inverse transformation 
        raytracer.translate(sphere, Vector3D(0, 0, -5));	
        // Similarly, add rotations and scaling 
        raytracer.rotate(sphere, 'x', -45); 
        raytracer.rotate(sphere, 'z', 45); 
        if (SCALESPHERE == 1 && DOMOTIONBLUR == 0)	raytracer.scale(sphere, Point3D(0, 0, 0), factor1);
        if (ONLYPARTONE != 1)
        {
            SceneDagNode* sphereTwo = raytracer.addObject(new UnitSphere(), &mud);
            // Add the translation matrix to the sphere's list of transformation and inverse transformation 
            raytracer.translate(sphereTwo, Vector3D(2, 2, -5));	
            // Similarly, add rotations and scaling 
            raytracer.rotate(sphereTwo, 'x', -30); 
            raytracer.rotate(sphereTwo, 'z', 25); 
            if (SCALESPHERE == 1)raytracer.scale(sphereTwo, Point3D(0, 0, 0), factor1);
        }
        if (DOMOTIONBLUR == 1) raytracer.scale(sphere, Point3D(0, 0, 0), factor3);
        // Similarly, do transformations on the plane 
        raytracer.translate(plane, Vector3D(0, 0, -7));	
        raytracer.rotate(plane, 'z', 45); 
        raytracer.scale(plane, Point3D(0, 0, 0), factor2);
    }
    // Do animation 
    else
    {
        double factorPlane[3] = {10, 10, 10};
        double factorSphere[3] = {1, 1, 1};
        double factorCyl[3] = { 0.5, 1, 0.5 };
        SceneDagNode* cylinder = raytracer.addObject(new UnitCylinder(), &cylMat);

        // Do transformations on the plane 
        /*raytracer.translate(plane, Vector3D(0, 0, -5));	
        raytracer.rotate(plane, 'x', -90);
        raytracer.scale(plane, Point3D(0, 0, 0), factorPlane);*/
        raytracer.translate(plane, Vector3D(0, 0, -7));	
        raytracer.rotate(plane, 'z', 45); 
        raytracer.scale(plane, Point3D(0, 0, 0), factorPlane);

        // Transform Sphere
        raytracer.translate(sphere, Vector3D(-2, 1.5, -5));	
        raytracer.rotate(sphere, 'x', -45); 
        raytracer.rotate(sphere, 'z', 45); 
            
        // Transform Cylinder		
        raytracer.translate(cylinder, Vector3D(-2, 0.5, -5));
        raytracer.scale(cylinder, Point3D(0,0,0), factorCyl);

        // Transform Sphere
        raytracer.translate(sphere, Vector3D(-3, 1.5, -7));	
        raytracer.rotate(sphere, 'x', -45); 
        raytracer.rotate(sphere, 'z', 45); 
            
        // Transform Cylinder		
        raytracer.translate(cylinder, Vector3D(-3, 0.5, -7));
        raytracer.scale(cylinder, Point3D(0,0,0), factorCyl);
    }
	//---------------------------------------------------------------------------------------------------------------
	// Step 4: Render the scene given the width, height, camera location, orientation, frame of view, and direction
	//---------------------------------------------------------------------------------------------------------------
	// Render the scene, feel free to make the image smaller for testing purposes.	
	// Note: Render creates and flushes the pixel buffer in 2 different functions 
    //       (it creates the pixel buffer in render() and flushes it in flushPixelBuffer() that is called by render() 
 
    // If don't do animation 
    if (DOANIMATION == 0)
    {
        raytracer.render(width, height, eye, view, up, fov, "view1.bmp", sphere);	// output the image with name view1.bmp
        
        std::cout << "View1.bmp done" << std::endl; 
        // Render it from a different point of view.
        Point3D eye2(4, 2, 1);
        Vector3D view2(-4, -2, -6);
        //Point3D eye2(2, 0, 1);
        //Vector3D view2(0, 0, -1);
        raytracer.render(width, height, eye2, view2, up, fov, "view2.bmp", sphere);	// output the image with name view2.bmp 

        std::cout << "View2.bmp done" << std::endl; 
        Point3D eye3(-1, 0, 0);
        Vector3D view3(1, 0, -1);
        raytracer.render(width, height, eye3, view3, up, fov, "view3.bmp", sphere);	// output the image with name
        std::cout << "View3.bmp done" << std::endl; 
        Point3D eye4(-2, 0, -6);
        Vector3D view4(1, 0, 0);
        raytracer.render(width, height, eye4, view4, up, fov, "view4.bmp", sphere);	// output the image with name
        std::cout << "View4.bmp done" << std::endl; 

    }
    // Do animation
    else
    {
            //Point3D eyeAnim(3, 1, 0);
            //Vector3D viewAnim(-2, 0, -1);
            //raytracer.render(width, height, eyeAnim, viewAnim, up, fov, "anim.bmp", sphere);	// output the image with name view1.bmp

            // Animation parameters
            const int NUM_FRAMES = 1;
            const int TRAVEL_DISTANCE = 10;
            const int NUM_LIGHTS = 2;
            const int LIGHT_HEIGHT = 10;
            const int CAM_SPEED = 0.1;
            char buffer[50]; // for keeping numbers of image names
            Point3D curEye;

        // FLYING THROUGH TEETH
        for(int i = 0; i < NUM_FRAMES; i++) {   
            // Move teeth together
            
            snprintf(buffer, sizeof(buffer), "%d%s", i, ".bmp");
            //raytracer.render(width, height, overViewEye, overViewView, overViewUp, overViewFOV, buffer, sphere);
            raytracer.render(width, height, curEye, view, up, fov, buffer, sphere);
            curEye[0] -= CAM_SPEED;
            curEye[1] -= CAM_SPEED;
            std::cout << buffer << " done" << std::endl; 
        }
        
            std::cout << "anim1.bmp done" << std::endl; 
    }
	return 0;
}

