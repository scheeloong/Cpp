#include "Node.h"

//-------------------------------------------------
// Constructors & Destructors
//-------------------------------------------------
// Default Constructor
Node::Node()
{
    nodeID = 0;
    setV = false;
    voltage = 0; // this is arbitrary based on last solve command as long setV is false
    numRes = 0;
    next = NULL;
    posX = -1; // default value
}

Node::Node(int nodeid_)
{
    nodeID = nodeid_;
    setV = false;
    voltage = 0;
    numRes = 0;
    next = NULL;
}

// Destructor for class Node
Node::~Node()
{
    // This Node destroys itself
}

//-------------------------------------------------
//----------
// Methods
//----------

// This method adds a resistor to the current Node given the resistor's label, resistance and its endpoints
// Parameters: string label_, double resistance_, int endpoints[2]
// Return Value: None
void Node::addResistor (string label_, double resistance_, int endpoints[2])
{
    this->resList.insertResistor(label_, resistance_, endpoints);
    this->numRes++;
}

// This method returns true if this node contains a given resistor label and false otherwise
// Parameters: string label_
// Return Value: True if head points to a resistor object, false if head points to NULL
bool Node::hasResistor(string label_)
{
    if (resList.hasResistor(label_))
        return true;
    else return false;
}

// This method modifies a resistor in Resistor List
// given a resistor name and its resistance
// Parameters: double resistance, string label
// Return value: True if modified, false otherwise
bool Node::modifyResistor(string label_, double resistance_)
{
    if (resList.modifyResistor(label_,resistance_))
        return true;
    else return false;
}


bool Node::modifyResistorPosY(string label_, double posY_)
{
    if (resList.modifyResistorPosY(label_,posY_))
        return true;
    else return false;
}

bool Node::getResistorPosY(string label_, double &posY_)
{
    if (resList.getResistorPosY(label_,posY_))
        return true;
    else return false;
}



// This method deletes a resistor in Resistor List given a resistor name
// Parameters: string label
// Return Value: true if Resistor is deleted, false otherwise
bool Node::deleteResistor (string label_)
{
    if (resList.deleteResistor(label_))
    {
        this->numRes--;
        return true;
    }
    else
    {
        return false;
    }
}

// This method deletes all resistors currently in the resistor list
// Parameters: None
// Return Value: None
void Node::deleteAllResistors()
{
    resList.deleteAllResistors();
    this->numRes = 0;
}

// This method returns resistance of resistor label_
// Parameters: string label_
// Return Value: Value of resistance of label_
double Node::getResistance(string label_)
{
    return resList.getResistance(label_);
}

// This method prints resistor with label label_
// Parameters: string label_
// Return Value: true if printed, false if resistor does not exist
bool Node::printResistor(string label_)
{
    return resList.printResistor(label_);
}

// This method prints out the connections of Resistor objects connected to this Node.
// Parameter: nodeIndex is the position of this node in the node array.
void Node::print()
{
    if (this->numRes != 0)
    {
    cout << "Connections at node " << this->nodeID <<": "
    << this->numRes << " resistor(s)" <<endl;
    resList.printResistors(); // print all the resistors, indentation is done in printResistors function
    }
}

//-------------------------------------------------
//-------
// GET
//-------

// This method gets the Node's object's member data nodeID.
int Node::getNodeID() const
{
    return this->nodeID;
}


// This method gets the Node's object's member data numRes.
int Node::getNumRes() const
{
    return this->numRes;
}

// This method gets the Node's object's member data voltage.
double Node::getVoltage() const
{
    return this->voltage;
}

// This method gets the Node's object's member data setV.
bool Node::getSetV() const
{
    return this->setV;
}

// This method gets the Node's object's member data next.
Node* Node::getNext() const
{
    return this->next;
}

// This method gets the Node's object's member data resList.
const ResistorList* Node::getResistorList() const
{
    return &(this->resList);
}

const double Node::getPosX() const
{
    return this->posX;
}


//-------------------------------------------------
//-------
// SET
//-------

// This method sets the Resistor's object's member data nodeID.
void Node::setNodeID(int nodeid_)
{
    this->nodeID = nodeid_;
}

// This method sets the Resistor's object's member data numRes.
void Node::setNumRes (int numRes_)
{
    this->numRes = numRes_;
}

// This method sets the Resistor's object's member data voltage
// and updates its setV status
void Node::setVoltage(double voltage_)
{
    this->voltage = voltage_;
    this->setV = true;
}

void Node::unsetVoltage()
{
    this->voltage = 0;
    this->setV = false;
}

// This method temporary sets the Resistor's object's member data voltage
// but keeps its setV status as false
void Node::setTemporaryVoltage(double voltage_)
{
    this->voltage = voltage_;
    this->setV = false;
    return;
}

// This method sets the Resistor's object's member data next.
void Node::setNext(Node* p)
{
    this->next = p;
}


void Node::setPosX(double posX_)
{
    this->posX = posX_;
}

