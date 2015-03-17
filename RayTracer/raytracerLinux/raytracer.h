/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		This file contains the interface and 
		datastructures of the raytracer.  
		Simple traversal and addition code to 
		the datastructures are given to you.

***********************************************************/

// This file contains information about: 
//---------------------
// LightListNode (Linked List of lights) 
// SceneDagNode (Directed Acyclic Graph)
// RayTracer
//---------------------

#include "util.h"
#include "scene_object.h"
#include "light_source.h"

// Linked list containing light sources in the scene.
struct LightListNode {
	// The data
	LightSource* light;
	// The next data it points to 
	LightListNode* next;
	// Initialize the light to point to nothing 
	LightListNode() : light(NULL), next(NULL) {}
	// Initialize a new light to point to next
	LightListNode( LightSource* light, LightListNode* next = NULL ) : 
		light(light), next(next) {}
	// Destroy a light Node
	~LightListNode() { 
		if (!light) delete light; 
	}
};

// The scene graph, containing objects in the scene.
//  Basically has an object poining to all other objects 
// The scene graph, containing objects in the scene.
// Dag => Directed Acyclic Graph Data structure for all the nodes of objects in the scene (scene == world) 
struct SceneDagNode {
	// Pointer to geometry primitive, used for intersection.
	SceneObject* obj;
	// Pointer to material of the object, used in shading.
	Material* mat;

	// Each node maintains a transformation matrix, which maps the 
	// geometry from object space to world space and the inverse.
	// Matrix Maps Object->World 
	// Matrix Maps Object->Inverse
	Matrix4x4 trans;
	Matrix4x4 invtrans;
	
	// This is a graph data structure
	// next is whatever that the current node does not point down towards (sibling) 
	// parent is the node pointing to the current node 
	// child is whatever the current node points to 
	// Internal structure of the tree, you shouldn't have to worry about them.
	SceneDagNode* next;
	SceneDagNode* parent;
	SceneDagNode* child;

	//-----------------------------------------------------------------------------------------------
	// Constructors
	//-----------------------------------------------------------------------------------------------
	SceneDagNode() : 
		obj(NULL), mat(NULL), 
		next(NULL), parent(NULL), child(NULL) {
	}	
	// Initialize an Object with Material property 
	SceneDagNode( SceneObject* obj, Material* mat ) : 
		obj(obj), mat(mat), next(NULL), parent(NULL), child(NULL) {}
	//-----------------------------------------------------------------------------------------------
	// Destructor 
	~SceneDagNode() 
	{
		if (!obj) delete obj;
		if (!mat) delete mat;
	}
};

class Raytracer {
private:
//----------------------------------------------------------------------------------------------
// Private Members
//----------------------------------------------------------------------------------------------
	// Viewport is the sreen generated on your desktop 
	// Width and height of the viewport.
	// Width and height of the viewport.
	int _scrWidth;
	int _scrHeight;

	// Light list and scene graph.
	LightListNode *_lightSource;
	SceneDagNode *_root;

	// Pixel buffer.
	unsigned char* _rbuffer;
	unsigned char* _gbuffer;
	unsigned char* _bbuffer;
	// Each of these is dynamically allocated with a single array that has 
	// a size big enough for the entire screen in chars 

	// Maintain global transformation matrices similar to OpenGL's matrix
	// stack.  These are used during scene traversal. 
	Matrix4x4 _modelToWorld;
	Matrix4x4 _worldToModel;

//----------------------------------------------------------------------------------------------
// Private Methods
//----------------------------------------------------------------------------------------------
	// Allocates and initializes the pixel buffer for rendering, you
	// could add an interesting background to your scene by modifying 
	// this function.
	void initPixelBuffer();

	// Saves the pixel buffer to a file and deletes the buffer.
	void flushPixelBuffer(char *file_name);

	// Return the colour of the ray after intersection and shading, call 
	// this function recursively for reflection and refraction.  
	Colour shadeRay( Ray3D& ray ); 

	// Constructs a view to world transformation matrix based on the
	// camera parameters.
	Matrix4x4 initInvViewMatrix( Point3D eye, Vector3D view, Vector3D up );

	// Traversal code for the scene graph, the ray is transformed into 
	// the object space of each node where intersection is performed.
	void traverseScene( SceneDagNode* node, Ray3D& ray );

	// After intersection, calculate the colour of the ray by shading it
	// with all light sources in the scene.
	void computeShading( Ray3D& ray );
	

public:
	//--------------------------------------------------------------------------------------------------------------------
	// Constructor 
	//--------------------------------------------------------------------------------------------------------------------
	Raytracer();
	~Raytracer();
	//--------------------------------------------------------------------------------------------------------------------

	//--------------------------------------------------------------------------------------------------------------------
	// AddObject()
	//--------------------------------------------------------------------------------------------------------------------
	// Add an object into the scene with a specific parent node, 
	// don't worry about this unless you want to do hierarchical 
	// modeling.  You could create nodes with NULL obj and mat, 
	// in which case they just represent transformations.  
	// So then it gets added to a different object as part of the linked list and not always the root
	SceneDagNode* addObject( SceneDagNode* parent, SceneObject* obj, Material* mat );

	// Add an object into the scene, with material mat.  
	// The function returns a Node to the object node you just added
 	// use the Node to apply transformations to the object at is in the linked list 
	SceneDagNode* addObject( SceneObject* obj, Material* mat ) 
	{
		return addObject(_root, obj, mat); // call the add object from parent above 
	}
	//--------------------------------------------------------------------------------------------------------------------

	// Add a light source to the ray tracer, it appends it to the linked list of lights 
	LightListNode* addLightSource( LightSource* light );


	LightListNode* addHeadLightSource(double posX, double posY, double posZ, double colorX, double colorY, double colorZ );

	// Transformation functions are implemented by right-multiplying 
	// the transformation matrix to the node's transformation matrix.
	
	// Apply rotation about axis 'x', 'y', 'z' angle degrees to node.
	// Done by multiplying transformation matrix to given node's transofrmation matrix 
	void rotate( SceneDagNode* node, char axis, double angle );

	// Apply translation in the direction of trans to node.
	void translate( SceneDagNode* node, Vector3D trans );

	// Apply scaling about a fixed point origin.
	// Therefore, to scale object itself, the node has to be already be transformed to the origin (object space) before using this function
	// Apply scaling about a fixed point origin.
	void scale( SceneDagNode* node, Point3D origin, double factor[3]);

	// Renders an image fileName with width and height and a camera positioned at eye, with view vector view, up vector up, and field of view fov.
	// Note: View Vector is direction camera is looking at (the tracing from the camera)
	// Note: Up Vector is when you rotate the camera, and image gets rotated although you are tracing towards the same point) 
	// Field Of View => angle at which you can see 
	void render( int width, int height, Point3D eye, Vector3D view, Vector3D up, double fov, char* fileName, SceneDagNode* node );
};
