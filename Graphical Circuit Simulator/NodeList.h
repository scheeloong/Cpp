#ifndef NODELIST_H_INCLUDED
#define NODELIST_H_INCLUDED

#include <string>
#include <iomanip>
#include <iostream>
#include <cmath> // for abs()
#include <sstream> // for concatenation of strings with ints and doubles
using namespace std; // Note: If forgot ; from here, can put it at top of Resistor.h and CodeBlocks will still compile
#include "Node.h"

// for solve() & draw()
#include "ResistorList.h"
#include "Resistor.h"
#include "easygl.h" // For WIN32 Graphics

#define MIN_ITERATION_CHANGE 0.0001 // for solve()

// for draw()
// These parameters can be changed for better graphics
#define NODE_WIDTH 10
#define RESISTOR_RADIUS NODE_WIDTH/3
#define DISTANCEBETWEENRESISTORS RESISTOR_RADIUS*4
#define DISTANCEBETWEENNODES NODE_WIDTH*2

extern easygl window;
class NodeList
{
private:
    Node* head;

    // Methods
    Node* findNode(int nodeID, Node* &prev); // used by deleteNode, findInsertNode, printNode
    Node* findInsertNode (int nodeID_); // used by insertResistor
    bool deleteNode (int nodeID_); // used by deleteResistor and deleteAllResistors
    double calcNodeVoltage(int nodeID_); // used by solve

public:
    // Constructors
    NodeList();

    // Destructors
    ~NodeList();

    // Methods
    bool hasResistor(string label_);
    void insertResistor(string label, double resistance,  int endpoints[2]);
    bool modifyResistor(string label_, double resistance_);
    bool deleteResistor(string label_);
    void deleteAllResistors();
    double getResistance(string label_);
    void setVoltage(int nodeID_, double voltage_);
    void unsetVoltage(int nodeID_);
    bool solvable();
    void solve();
    void solveDraw();

    // Prints
    bool printResistor(string label_);
    void printAllNodes();
    bool printNode(int nodeID);
    void printNodeVoltages(); // used by solve

    // Draw
    void set_draw_coords (float & xleft, float & ybottom, float & xright, float & ytop);
    void draw(float ytop_);

};
#endif // NODELIST_H_INCLUDED
