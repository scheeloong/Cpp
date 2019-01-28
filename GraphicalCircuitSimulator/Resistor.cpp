#include "Resistor.h"

//-------------------------------------------------
// Constructors & Destructors
//-------------------------------------------------
// Constructor for class Resistor
// Parameters:
//      label_ is used to identify the resistor
//      resistance_ is the resisance of this resistor
//      endpoints_ holds the nodeIDs to which this resistor connects
//      next is a pointer to the next resistor added into the resistor list
Resistor::Resistor(string label_, double resistance_, int endpoints_[2], Resistor* next)
{
    this->resistance = resistance_;
    this->label = label_;
    this->endpoints[0] = endpoints_[0];
    this->endpoints[1] = endpoints_[1];
    this->next = next; // Note: Did not create a new object here, simply points to an object given
    this->posY = -1; // Default vaue
}

// Note: Default constructor does not work anymore
// Thus, re-initialize the default constructor that does nothing.
Resistor::Resistor()
{
    // do nothing
}

// Note: Destructor of ResistorList destroys the whole linked list
// Destructor for class Resistor
Resistor::~Resistor()
{
    // This destructor only destroys this resistor node itself
    // Do nothing
}

//-------------------------------------------------
// Overload Assignment (=) operator
Resistor& Resistor::operator = (const Resistor right_side)
{
    this->resistance = right_side.resistance;
    this->label = right_side.label;
    this->endpoints[0] = right_side.endpoints[0];
    this->endpoints[1] = right_side.endpoints[1];
    // shallow copy, need to handle making a new copy in ResistorList's assignment operator
    this->next = right_side.next;
    return *this;
}

//-------------------------------------------------
//----------
// Methods
//----------

// This method returns if the given Resistor is pointing to another Resistor
// Parameters: None
// Return Value: True if it is, False otherwise
bool Resistor::hasNext()
{
    if (this->next == NULL)
        return false;
    return true;
}

//-------------------------------------------------
//-------
// GET
//-------

// This method gets the Resistor's object's member data name.
string Resistor::getLabel() const
{
    return this->label;
}

// This method gets the Resistor's object's member data resistance.
double Resistor::getResistance() const
{
    return this->resistance;
}

// return pointer to next node
Resistor* Resistor::getNext() const
{
    return this->next;
}

// return pointer to endpoints
const int* Resistor::getEndPoints() const
{
    return this->endpoints;
}

const double Resistor::getPosY() const
{
    return this->posY;
}



//-------------------------------------------------
//-------
// SET
//-------

// This method sets the Resistor's object's member data resistance.
void Resistor::setLabel(string label_)
{
    this->label = label_;
}

// This method sets the Resistor's object's member data resistance.
void Resistor::setResistance(double resistance_)
{
    this->resistance = resistance_;
}
// This method sets the Resistor's object's member data endpoints_.
void Resistor::setNodes(int endpoints_[2])
{
    this->endpoints[0] = endpoints_[0];
    this->endpoints[1] = endpoints_[1];
}

void Resistor::setNext(Resistor* p) // set the nextNode to point to something else
{
    this->next = p;
}
// Print method for class Resistor
void Resistor::print()
{
    cout << setw(20) << left << this->label << " " <<
    setw(8) << right << showpoint << setprecision(2) << this->resistance << " Ohms "
    << this->endpoints[0] << " -> " << this->endpoints[1] << endl;
}

void Resistor::setPosY(double posY_)
{
    this->posY = posY_;
}
