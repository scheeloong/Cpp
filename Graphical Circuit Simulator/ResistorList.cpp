#include "ResistorList.h"

// Constructors
ResistorList::ResistorList()
{
    this->head = NULL;
}

// Destructor
ResistorList::~ResistorList ()
{
    Resistor * del;
    while (this->head != NULL)
    {
        del = (this->head);
        this->head = this->head->getNext(); // point to next element
        delete del;
    }
}

//------------
// Methods
//------------
// This method finds a resistor in the resistor list by resname and
// returns a pointer to it if it exists.
// For special cases (first node, not found), return value of prev is insignificant
// Parameters: string resname, Resistor * &prev
// Return Value: pointer to resistor if it exists, NULL otherwise
// Return by Reference: prev is a pointer to the resistor before the found resistor
Resistor* ResistorList::findResistor(string resname, Resistor * &prev)
{
    if (head == NULL) return NULL;
    Resistor *found = head;
    prev = found;
    // Find the resistor in the ResistorList Linked List
    while(found != NULL)
    {
        // If resistor found
        if(found->getLabel() == resname)
        {
            return found;
        }
        // Continue looking for resistor
        prev = found;
        found = found->getNext();
    }
    // Resistor not found
    return NULL;
}

// This method inserts a resistor object at the end of the resistor list
// Parameters: *Resistor res
// Return Value: None
void ResistorList::insertResistor(string label, double resistance,  int endpoints[2])
{
    // Note: Resistor needs to be printed in order added to node,
    //       append new resistor to the end of list for easy print()
    // Make a new resistor object using Resistor's constructor
    Resistor* res = new Resistor(label, resistance, endpoints, NULL);
    // If inserting first resistor
    if (this->head == NULL)
    {
        this->head = res;
        return;
    }
    // Head already has a resistor element
    Resistor* last;
    last = this->head;
    // make last point to last resistor in list
    while (last->getNext() != NULL)
    {
        last = last->getNext();
    }
    // here, last->next points to NULL,
    // which means we found the last resistor in list
    // add new resistor to end of list
    last->setNext(res);
    return;
}

// This method modifies a resistor in Resistor List
// given a resistor name and its resistance
// Parameters: double resistance, string label
// Return value: True if modified, false otherwise
bool ResistorList::modifyResistor(string label_, double resistance_)
{
    Resistor* temp = NULL;
    Resistor* res = this->findResistor(label_, temp);
    if (res==NULL)  return false; // no resistor of name found
    // Resistor was found
    res->setResistance(resistance_);
    return true;
}

bool ResistorList::modifyResistorPosY(string label_, double posY_)
{
    Resistor* temp = NULL;
    Resistor* res = this->findResistor(label_, temp);
    if (res==NULL)  return false; // no resistor of name found
    // Resistor was found
    res->setPosY(posY_);
    return true;
}

bool ResistorList::getResistorPosY(string label_, double &posY_)
{
    Resistor* temp = NULL;
    Resistor* res = this->findResistor(label_, temp);
    if (res==NULL)  return false; // no resistor of name found
    // Resistor was found
    posY_ = res->getPosY();
    return true;
}


// This method deletes a resistor in Resistor List given a resistor name
// Parameters: string label
// Return Value: true if Resistor is deleted, false otherwise
bool ResistorList::deleteResistor (string label_)
{
    Resistor * prev = NULL;
    Resistor* res = this->findResistor(label_, prev);
    // If resistor not found, return false
    if (res==NULL)  return false;
    // Resistor found,
    // re-connect link list without this resistor
    // Cases:
    // i) If deleting first resistor node in linkedlist
    if (prev == head && res == head)
    {   // Note: Must always have this first case. If do not reconnect head to NULL when
        // deleting whatever head points to, will cause crash during run time.
        head = res->getNext(); // res->getNext() might be NULL if only element in list
    }
    // ii) if not deleting first node
    else
    {
        prev->setNext(res->getNext());
    }
    delete res;
    return true;
}

// This method deletes all resistors currently in the resistor list
// Parameters: None
// Return Value: None
void ResistorList::deleteAllResistors()
{
    Resistor* res = head;
    Resistor* del = res;
    // If no resistor of name found, return false
    while (res!=NULL)
    {
        del = res;
        res = res->getNext();
        delete del;
    }
    // Point head back to NULL
    head = NULL;
    return;
}

// This method returns true if this resistor list contains a given resistor label and false otherwise
// Parameters: string label_
// Return Value: True if head points to a resistor object, false if head points to NULL
bool ResistorList::hasResistor(string label_)
{
    Resistor * res;
    Resistor *temp;
    res = this->findResistor(label_, temp);
    if(res == NULL) return false;
    else return true;
}

// This method returns resistance of resistor label_
// It will only be called if resistorList has the resistor
// and that is checked in Main.cpp
// Parameters: string label_
// Return Value: Value of resistance of label_
double ResistorList::getResistance(string label_)
{
    Resistor *res;
    Resistor *temp;
    res = this->findResistor(label_, temp);
    return res->getResistance();
}

// This method gets the ResistorList member data head
Resistor* ResistorList::getHead() const
{
    return this->head;
}



// This method prints resistor label_
// Parameters: string label_
// Return Value: true if printed, false if resistor does not exist
bool ResistorList::printResistor(string label_)
{
    Resistor *temp;
    Resistor * res = this->findResistor(label_,temp);
    if (res== NULL) return false;
    cout <<"Print:"<<endl;
    res->print();
    return true;
}

// This method prints all the resistors connected to this Resistor List
// Parameters: None
// Return Value: None
void ResistorList::printResistors()
{
    Resistor * res = head;
    while (res!= NULL)
    {
        cout << "  "; // indent 2 spaces for printing resistors from printNode
        res->print();
        res = res->getNext();
    }
}
