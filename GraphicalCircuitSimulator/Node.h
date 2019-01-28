#ifndef NODE_H
#define NODE_H

#include <string>
#include <iostream>
#include <iomanip>
#include "ResistorList.h" // each node contains a Resistor List

using namespace std;

class Node
{
private:
   int nodeID; // node ID of this node
   double voltage; // if setV is false, this is arbitrary based on last solve command
                   // if setV is true, this is a fixed constant
   bool setV;
   int numRes; // number of resistors of this node
   ResistorList resList; // each node has only 1 resistor list
   Node* next;
   double posX; // for draw

// Note: (From Discussion Board)
//   Case1: Maintaining any setV on voltages
//      Nodes are only deleted from "deleteR name" if no resistors and setV is false
//      Nodes are not deleted if no resistors and setV is true if
//   Case2: "deleteR all"
//      All nodes are deleted regardless of setV

public:
   // Constructor and Destructor
   Node();
   Node(int nodeid_);
   ~Node();

   // Methods
   void addResistor (string label_, double resistance_, int endpoints[2]);
   bool hasResistor(string label_);
   bool modifyResistor(string label_, double resistance_);
   bool deleteResistor (string label_);
   void deleteAllResistors();
   double getResistance(string label_);

   bool modifyResistorPosY(string label_, double posY_);
   bool getResistorPosY(string label_, double &posY_);
   // Get
   double getVoltage() const;
   bool getSetV() const;
   int getNumRes() const;
   Node* getNext() const; // return pointer to next node
   int getNodeID() const;
   const ResistorList* getResistorList() const;
   const double getPosX() const;

   // Set
   void setTemporaryVoltage(double voltage_); // change voltage values for unset voltages
   void setVoltage(double voltage_);// changes setV here , if 0, setV is false, if not 0, setV is true
   void setNumRes (int numRes_);
   void setNext(Node* p); // set the nextNode to point to something else
   void setNodeID(int nodeid_);
   void setPosX(double posX_);
   void unsetVoltage();


   // Print
   bool printResistor(string label_);
   void print ();
};
#endif
