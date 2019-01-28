#ifndef RESISTORLIST_H_INCLUDED
#define RESISTORLIST_H_INCLUDED

#include <string>
#include <iostream>
#include <iomanip>

#include "Resistor.h"

using namespace std;

class Resistor; // let the compiler now this class is defined

class ResistorList
{
private:
    Resistor* head;
    // Methods
    Resistor* findResistor(string resname, Resistor * &prev);

public:
    // Constructor
    ResistorList ();

    // Destructor
   ~ResistorList ();

    // Methods
    void insertResistor(string label, double resistance,  int endpoints[2]);
    bool modifyResistor(string label_, double resistance_);
    bool deleteResistor (string label_);
    bool hasResistor(string label_);
    double getResistance(string label_);
    bool printResistor(string label_);
    void printResistors();
    void deleteAllResistors();

    bool modifyResistorPosY(string label_, double posY_);
    bool getResistorPosY(string label_, double &posY_);
    // Get
    Resistor* getHead() const;
};
#endif // RESISTORLIST_H_INCLUDED
