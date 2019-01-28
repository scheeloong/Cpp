#ifndef RESISTOR_H_INCLUDED
#define RESISTOR_H_INCLUDED

#include <iostream>
#include <string>
#include <iomanip>

using namespace std;

class Resistor
{
private:
   double resistance; // resistance (in Ohms)
   string label; // name of this resistor
   int endpoints[2]; // IDs of nodes it attaches to
   Resistor* next;
   double posY; // for draw()

public:
   // Constructor
   Resistor(string label_, double resistance_, int endpoints_[2], Resistor* next);
   // Default Constructor
   Resistor();
   // Destructor (delete all pointers to next
   ~Resistor();

   // Method
   bool hasNext();

   // Get
   string getLabel() const;
   double getResistance() const;
   const int* getEndPoints() const;
   Resistor* getNext() const; // returns pointer to next node
   const double getPosY() const;

   // Set
   void setLabel(string label_);
   void setResistance(double resistance_);
   void setNodes(int endpoints_[2]);
   void setNext(Resistor* p);
   void setPosY(double posY_);


   // Print
   // you *may* create either of the below to print your resistor
   // Dear TA, I used print() because I am reusing my code from last lab (LAB 3).
   void print ();
   friend ostream& operator<<(ostream&,const Resistor&);

   // Overload = operator
   Resistor& operator = (const Resistor right_side);
};
#endif
