#include "NodeList.h"

//-------------------------------------------------
// Constructors & Destructors
//-------------------------------------------------

// Constructor
NodeList::NodeList() { head = NULL;}

// Destructor
NodeList::~NodeList()
{
    Node * del;
    while (this->head != NULL)
    {
        del = (this->head);
        this->head = this->head->getNext(); // point to next element
        delete del;
    }
}

//-------------------------------------------------
//----------
// Methods
//----------

// This method accepts a node ID and return a pointer to the corresponding Node, or NULL if it does not exist
// Parameters: int nodeid, Node* & prev
// Return Value: pointer to node if it exists, NULL it does not.
//               If it exists, return by reference the value of pointer before the Node
Node* NodeList::findNode(int nodeID, Node* &prev)
{
    if (head == NULL) return NULL;
    Node* found = head;
    prev = found;
    while (found!= NULL)
    {
        if (found->getNodeID() == nodeID)
            return found;
        prev = found;
        found = found->getNext();
    }
    // No more Nodes which means Node is not found
    return found; // NULL
}

// This method accepts a node ID and return a pointer to the corresponding Node if it exists, or creates a new one
// Parameters: int nodeID
// Return Value: pointer to node if it exists, or a pointer to the new created Node
Node* NodeList::findInsertNode (int nodeID)
{
    Node *temp;
    Node *newnode;
    newnode = this->findNode(nodeID, temp);
    // If given nodeID not found, insert Node
    if(newnode == NULL)
    {
        newnode = new Node(nodeID);
        if (this->head == NULL)
        {
            this->head = newnode;
            return newnode;
        }
        // head already has node, insert according to value of nodeID
        Node* position;
        position = this->head;
        Node* prev = position;
        // Look for proper position
        while (((position->getNext()) != NULL ))
        {
            if (position->getNodeID() > newnode->getNodeID())
            {
                // Found proper position, insert it there
                // If proper position is first node,
                if (prev == this->head && position == this->head)
                {
                    newnode->setNext(position);
                    head = newnode;
                    return newnode;
                }
                // Proper position is not first node
                newnode->setNext(position);
                prev->setNext(newnode);
                return newnode;
            }
            prev = position;
            position = position->getNext(); // make position point to proper position in list
        }
        // Proper position is at the last Node,
        // which might be the first node
        // Check if last node is bigger or smaller than new node
        // If last node is bigger than new node, insert it before
        if (position->getNodeID() > newnode->getNodeID())
        {
            // Note: position's NodeID will never be the same as a newnode's NodeID

            // Found proper position before last element, insert it
            // If first element is last
            if (prev == this->head && position == this->head)
            {
                newnode->setNext(position);
                this->head = newnode; // make p the first in list
            }
            // If first element is not last
            else
            {
                newnode->setNext(prev->getNext());
                prev->setNext(newnode);
            }
            return newnode;
        }
        // If last node is smaller than new node, insert it after
        else
        {
            newnode->setNext(position->getNext()); // will point to NULL
            position->setNext(newnode);
            return newnode;
        }
    }
    else
    {
        // node existed, return node
        return newnode;
    }
}

// This method determines if a given nodeID already exists in the NodeList and delete it if it does
// Parameters: int nodeID_
// Return Value: true if node exist and is deleted, false otherwise
bool NodeList::deleteNode (int nodeID_)
{
    Node * prevNode = NULL;
    Node * delNode = this->findNode(nodeID_, prevNode);
    // If given nodeID is not found, there is nothing to delete
    if(delNode == NULL)
        return false;
    // If deleting first node
    if (prevNode == head && delNode == head)
    {
        head = delNode->getNext();
    }
    // nodeID is found, delete it after rearranging list
    // If not deleting first Node
    else
    {
        prevNode->setNext(delNode->getNext());
    }
    delete delNode;
    return true;
}

// This method returns value of the newly calculated NodeVoltage for a Node in this NodeList
// Parameters: nodeID_
// Return Value: new calculated voltage value for given nodeID_
double NodeList::calcNodeVoltage(int nodeID_)
{
    // p is a Node pointer pointing to whatever the Node's head is pointing to
    Node* p = this->findInsertNode(nodeID_);
    Node* q;
    const ResistorList* resList;
    Resistor* res;
    const int* nodeid;
    resList = p->getResistorList();
    res = resList->getHead();
    // Initialize Values for Kirchoff's Current Equation
    double finalVoltage = 0; // V = IR = I/G
    double finalConductance = 0;  // G = 1/R
    double finalResistance = 0;   // R = 1/G
    double finalCurrent = 0;      // I'
    double resistanceTemp = 0;
    double voltageTemp = 0;
    int otherNodeID = 0;

    // Calculate total conductance and total current
    while(res!=NULL)
    {
        // Step 1: Add conductance for each resistor
        resistanceTemp = res->getResistance();
        finalConductance += (1/resistanceTemp) ;

        // Step 2: Add current for each resistor
        // Get other NodeID that this resistor connects to
        nodeid = res->getEndPoints();
        if (nodeid[0] == nodeID_)
        {
            otherNodeID = nodeid[1];
        }
        // nodeid[1] == nodeID_
        else
        {
            otherNodeID = nodeid[0];
        }
        // Get otherNodeID's Node and its voltage
        q = this->findInsertNode(otherNodeID);
        voltageTemp = q->getVoltage();
        // Add to current finalCurrent
        finalCurrent += (voltageTemp/resistanceTemp);
        // Step 3: Iterate & Repeat till no more resistors
        res = res->getNext();
    }
    // R = 1/G
    finalResistance = 1 / finalConductance;
    // V = IR
    finalVoltage = finalCurrent * finalResistance ;
    return finalVoltage;
}

// This method determines if a resistor with a given label exists in any of the Nodes
// Parameters: string label_
// Return Value: true if the resistor exist, false otherwise
bool NodeList::hasResistor(string label_)
{
    if (head == NULL) return false; // no nodes in node list
    Node* traverse = head;

    while(traverse != NULL)
    {
        if(traverse->hasResistor(label_))
            return true;
        traverse = traverse->getNext(); // check next Node
    }
    return false;
}

// This method inserts a resistor object at the end of the resistor list of each node it connects to.
// Parameters: string label, double resistance, int endpoints[2])
// Return Value: None
void NodeList::insertResistor(string label, double resistance,  int endpoints[2])
{
    // Get pointers to the node IDs
    Node *firstNode =  this->findInsertNode(endpoints[0]); // returns node pointer to this nodeid
    Node *secondNode = this->findInsertNode(endpoints[1]); // returns node pointer to this nodeid
    firstNode->addResistor(label, resistance, endpoints);
    secondNode->addResistor(label, resistance, endpoints);
    return;
}

// This method change the resistance of a resistor by name or returns false otherwise
// Parameters: string label_, double resistance_
// Return Value: true if the resistor with given name's resistance has been changed, false otherwise
bool NodeList::modifyResistor(string label_, double resistance_)
{
    Node* traverse = head;
    while(traverse != NULL)
    {
        if(traverse->modifyResistor(label_,resistance_))
        {
            // Here, Changed resistor on its first node
            // need change resistor on its second node
            traverse = traverse->getNext();
            // loop through other nodes until modified the same resistor at other places
            while(!(traverse->modifyResistor(label_, resistance_)))
            {
                traverse = traverse->getNext();
            }
            return true;
        }
        traverse = traverse->getNext(); // check next Node
    }
    return false;
}

// This method deletes a Resistor by name or returns false otherwise
// Parameters: string label_
// Return Value: true if resistor is deleted, false otherwise
bool NodeList::deleteResistor(string label_)
{
    Node* traverse = head;
    while(traverse != NULL)
    {
        if(traverse->deleteResistor(label_))
        {
            Node* del = traverse;
            // Here, Deleted resistor on its first node
            // need delete resistor on its second node
            traverse = traverse->getNext();
            // If first node now has no resistors and no set voltage
            if ((del->getNumRes()== 0) && (del->getSetV() == false))
            {
                // remove first node from linked list
                this->deleteNode(del->getNodeID());
            }
            // loop through other nodes until modified the same resistor at other places
            while(!(traverse->deleteResistor(label_)))
            {
                traverse = traverse->getNext();
            }
            del = traverse;
            // If second node now has no resistors and no set voltage
            if ((del->getNumRes()== 0) && (del->getSetV() == false))
            {
                // remove second node from linked list
                this->deleteNode(del->getNodeID());
            }
            return true;
        }
        traverse = traverse->getNext(); // check next Node
    }
    return false;
}

// This method deletes all resistors in the list
// Parameters: None
// Return Value: None
void NodeList::deleteAllResistors()
{
    Node* traverse = head;
    Node* del;
    while (traverse != NULL)
    {
        traverse->deleteAllResistors();
        del = traverse;
        traverse = traverse->getNext();
        this->deleteNode(del->getNodeID()); // remove this node from list
    }
    return;
}

// This method returns resistance of resistor label_
// Parameters: string label_
// Return Value: Value of resistance of label_
double NodeList::getResistance(string label_)
{
    Node* traverse = head;
    double resistance = 0;
    while (traverse!= NULL)
    {
        if (traverse->hasResistor(label_))
        {
            resistance = traverse->getResistance(label_);
            break;
        }
        traverse = traverse->getNext();
    }
    return resistance;
}


// This method sets the voltage of a given nodeID_
// It creates the Node if the nodeID_given doesn't exist
// Parameters:int nodeID_ double voltage_
// Return Value: None
void NodeList::setVoltage(int nodeID_, double voltage_)
{
    Node* node = this->findInsertNode(nodeID_);
    node->setVoltage(voltage_);
    return;
}

// This method unsets the voltage of a given nodeID_
// It destroys the Node if it does not have any resistors attached to it
// Parameters:int nodeID_
// Return Value: None
void NodeList::unsetVoltage(int nodeID_)
{
    Node* node = this->findInsertNode(nodeID_);
    node->unsetVoltage();

    // if this node has no resistors attached to it, delete it
    if (node->getNumRes() == 0)
    {
        this->deleteNode(nodeID_);
    }
    return;
}

// This method returns true if at least one node with resistors attached has its voltage set
// and false otherwise
// Parameters: None
// Return Value: true if solvable, false if no nodes with resistors have its voltage set.
bool NodeList::solvable()
{
    // If no nodes, return false
    if (head == NULL) return false;
    Node* p = head;
    while (p!= NULL)
    {
        // If this node has resistors attached,
        if (p->getNumRes() > 0)
        {
            // If the node with resistors attached has its voltage set
            if (p->getSetV() == true)
            {
                return true;
            }
        }
        p = p->getNext();
    }
    // Here p is NULL and no node with resistors attached
    // and has its voltage set can be found
    return false;
}

// This method solves the voltages of all the nodes on this nodeList
// and prints out each's (node with resistor)'s voltage value
// Note: This method is only called once it is determined that it is solvable in Main
// Parameters: None
// Return Value: None
void NodeList::solve()
{
    // intialize change to a number larger than MIN_ITERATION_CHANGE
    double change = MIN_ITERATION_CHANGE + 1;
    // To calculate how much a voltage of a node has change
    double oldVoltage = 0;
    double newVoltage = 0;
    double newChange = 0;
    Node *p = head;
    // initialize all node voltages of unset voltages to 0
    while (p!= NULL)
    {
        // If current node has no set voltages, initialize its voltage to 0
        if (p->getSetV() == false)
        {
            p->setTemporaryVoltage(0);
        }
        // Iterate
        p = p->getNext();
    }
    // Finish initializing node voltages of unsetVoltages to 0
    p = head;
    // Apply numerical calculation to solve for node Voltages
    while (change > MIN_ITERATION_CHANGE)
    {
        // Initialize change to less than MIN_ITERATION_CHANGE
        change = MIN_ITERATION_CHANGE - 1;
        while (p!= NULL)
        {
            // If current node has no set voltages, calculates its voltage change
            if (p->getSetV() == false)
            {
                // Get old voltage
                oldVoltage = p->getVoltage();
                newVoltage = this->calcNodeVoltage(p->getNodeID());
                p->setTemporaryVoltage(newVoltage);
                newChange = abs((newVoltage - oldVoltage));
                if (newChange > change)
                {
                    // Update change to the largest changes in voltages
                    change = newChange;
                }
            }
            // Iterate
            p = p->getNext();
        }
        p = head;
    }
    cout << "Solve:" << endl;
    this->printNodeVoltages();
    return;
}

// Solve function for draw
void NodeList::solveDraw()
{
    // intialize change to a number larger than MIN_ITERATION_CHANGE
    double change = MIN_ITERATION_CHANGE + 1;
    // To calculate how much a voltage of a node has change
    double oldVoltage = 0;
    double newVoltage = 0;
    double newChange = 0;
    Node *p = head;
    // initialize all node voltages of unset voltages to 0
    while (p!= NULL)
    {
        // If current node has no set voltages, initialize its voltage to 0
        if (p->getSetV() == false)
        {
            p->setTemporaryVoltage(0);
        }
        // Iterate
        p = p->getNext();
    }
    // Finish initializing node voltages of unsetVoltages to 0
    p = head;
    // Apply numerical calculation to solve for node Voltages
    while (change > MIN_ITERATION_CHANGE)
    {
        // Initialize change to less than MIN_ITERATION_CHANGE
        change = MIN_ITERATION_CHANGE - 1;
        while (p!= NULL)
        {
            // If current node has no set voltages, calculates its voltage change
            if (p->getSetV() == false)
            {
                // Get old voltage
                oldVoltage = p->getVoltage();
                newVoltage = this->calcNodeVoltage(p->getNodeID());
                p->setTemporaryVoltage(newVoltage);
                newChange = abs((newVoltage - oldVoltage));
                if (newChange > change)
                {
                    // Update change to the largest changes in voltages
                    change = newChange;
                }
            }
            // Iterate
            p = p->getNext();
        }
        p = head;
    }
    return;
}

//-------------------------------------------------
//----------
// Prints
//----------

// This method print resistor of label label_ if it exists and false otherwise
// Parameters: string label_
// Return Value: true if resistor exists and printed, false otherwise
bool NodeList::printResistor(string label_)
{
    // Note: Same resistor will be connected to 2 different nodes,
    // just need print from first node
    Node* traverse = head;
    while (traverse != NULL)
    {
        // Just need reach first node that the resistor connects to, then return true
        if (traverse->printResistor(label_))
            return true;
        traverse = traverse->getNext();
    }
    return false;
}

// This method prints all Nodes connected to this NodeList
// Parameters: None
// Return Value: None
void NodeList::printAllNodes()
{
    Node* traverse = head;
    // If there is a node, it means there must be at least a resistor in network
    // Always cout "Print: \n" even if there is no resistor
    // or nodes in network as shown in lab handout
    cout <<"Print:"<<endl;
    while (traverse != NULL)
    {
        traverse->print();
        traverse = traverse->getNext();
    }
    return;
}

// This method prints the Node with node nodeID
// It returns true if nodeID exists and printed or false otherwise
// Parameters: int nodeID
// Return Value: true if node exist and printed, false otherwise
bool NodeList::printNode(int nodeID)
{
    Node* temp;
    Node* p = findNode(nodeID, temp);
    // From discussion board, "Print:" gets outputted even if the node doesnt exist
    cout <<"Print:"<< endl;
    if (p == NULL) // Node does not exist
    {
        return false;
    }
    else
    {
        p->print();
        return true;
    }
}

// This method prints the voltage values of all (nodes with resistors attached)
// Parameters: None
// Return Value: None
void NodeList::printNodeVoltages()
{
    Node* p = head;
    while (p!=NULL)
    {
        // if this node has resistors
        if ((p->getNumRes()) > 0)
        {
            // Print out its node voltage value
            cout <<"  Node " << p->getNodeID() <<": " << fixed << setprecision(2) << showpoint << p->getVoltage() <<" V"<<endl;
        }
        // Iterate
        p = p->getNext();
    }
}

//-------------------------------------------------
//----------
// Draw
//----------

// Note: This code assumes every resistor inserted into network have a unique name
// Note: This code only draws nodes with resistors attached.
void NodeList::set_draw_coords (float & xleft, float & ybottom, float & xright, float & ytop)
{
    // Iterate over Nodes and ressistors, setting drawing coordinate member variables
    // Return minimum and maximum coordinates will use, so can be passed to
    // easy gl package via set_world_coordinates()
    // Initialize just in case no nodes in list
    xleft = 0;
    ybottom = 0;
//-------------------------------
//  NO NODES (DEFAULT VALUE)
//-------------------------------
    xright = 1000;
    ytop = 1000;
    if (head == NULL) return ; // no nodes to draw
// Nodes exist, calculate draw coordinates to be set
    Node* p = head;
    const ResistorList* resList;
    Resistor* res;
    resList = p->getResistorList();
    res = resList->getHead();
//-------------------------------
// INITIALIZE TO -1
//-------------------------------
// Step 1:Initialize all resistors and nodes to -1 coordinates.
    while (p!= NULL)
    {
        // Set the Node's X position to -1
        p->setPosX(-1);
        // If this node has resistors attached,
        if (p->getNumRes() > 0)
        {
            // Get its resistor list
            resList = p->getResistorList();
            res = resList->getHead();
            // Initialize all resistor's y position to -1 position
            while (res!= NULL)
            {
                res->setPosY(-1);
                res = res->getNext();
            }
        }
        p = p->getNext();
    }
//---------------------------------------------------------------------------------------
//  Compute X for every Node while computing Y for every Resistor in the Node
//---------------------------------------------------------------------------------------
// Step 2: Compute X for every Node
    p = head;
    resList = p->getResistorList();
    res = resList->getHead();
    double coordinateX = 10; // initial value for first Node
    double coordinateY = 10 + RESISTOR_RADIUS; // initial value for first Resistor in first Node
    string resName;
    Node * traverse;
    while (p!= NULL)
    {
        // If this node has resistors attached,
        if (p->getNumRes() > 0)
        {
            // Set the Node's X position to its current position
            p->setPosX(coordinateX);
            // Get its resistor list
            resList = p->getResistorList();
            res = resList->getHead();
            // Step 3: Compute Y for every Resistor
            while (res!= NULL)
            {
                // If this resistor hasn't been given a Y-Coordinate yet
                if (res->getPosY() == -1)
                {
                    // Set first resistor to coordinateY
                    res->setPosY(coordinateY);
                    // Need set second resistor of the same resistor to coordinateY
                    resName = res->getLabel();
                    traverse = p->getNext();
                    // loop through other nodes until modified the same resistor at other places
                    while(!(traverse->modifyResistorPosY(resName,coordinateY)))
                    {
                        traverse = traverse->getNext();
                    }
                    // updated both resistors coordinate Y
                    coordinateY += DISTANCEBETWEENRESISTORS;
                }
                res = res->getNext();
            }
            coordinateX += DISTANCEBETWEENNODES; // increment coordinateX
        }
        p = p->getNext();
    }
    // Here, we have finished updating all Nodes and Resistor Coordinates
    xright = coordinateX;
    ytop = coordinateY;
    return;
}

// Uses global variable window
void NodeList::draw(float ytop_)
{
    // Initialize window
    window.gl_clearscreen();
    window.gl_setfontsize (10);
    window.gl_setlinestyle (SOLID);
    window.gl_setlinewidth (1);
    // Initialize variables
    Node* p;
    Node* traverse;
    Resistor * res;
    const ResistorList * resList;
    string NodeText;
    string ResText;
    string resName;
    stringstream s;
    // If this is solvable, it means there are voltage set,
    if (this->solvable())
    {
        // solve the voltages without printing them out
        this->solveDraw();
    }
    p = head;
    // Step 1: Draw each Node with Resistors as GREY rectangles
    window.gl_setcolor(LIGHTGREY);
    while (p!= NULL)
    {
        if (p->getNumRes() > 0)
        {
            resList = p->getResistorList();
            res = resList->getHead();
            window.gl_fillrect(p->getPosX(), 10, p->getPosX() + NODE_WIDTH, ytop_); // (x1, y1, x2, y2)
            // Step 2: Set up the text and draw text for current Node
            window.gl_setcolor(BLACK);
            s << "Node " << p->getNodeID() ;
            NodeText = s.str();
            // Reinitialize for next
            s.clear(); s.str(string());
            window.gl_drawtext(p->getPosX() + (NODE_WIDTH/2), 5, NodeText, 250);
            s << fixed << setprecision(2) << showpoint << p->getVoltage() << " V";
            NodeText = s.str();
            // If this node has set voltage, write it in RED
            if (p->getSetV()) window.gl_setcolor(RED);
            window.gl_drawtext(p->getPosX() + (NODE_WIDTH/2), 2.5, NodeText, 250);
            s.clear(); s.str(string());
                // (xcenter, ycenter, text, bound
            // Step 3: Draw each Resistor as BLACK Circles
            window.gl_setcolor(BLACK);
            while (res != NULL)
            {
                window.gl_drawarc(p->getPosX() + (NODE_WIDTH/2) , res->getPosY(), RESISTOR_RADIUS, 0, 360);
                                   // (xcenter, ycenter, radius, startangle, angleextent)
                res = res->getNext();
            }
        }
        // Iterate through Nodes
        p = p->getNext();
        window.gl_setcolor(LIGHTGREY);
    }
    // No more nodes to draw
    // Step4: Draw RED lines between each resistor
    window.gl_setcolor(RED);
    p = head;
    double x1, x2, y1, y2; // for drawing Resistor Lines
    y2 = 0; // initialize y2 to 0
    while (p != NULL)
    {
        if (p->getNumRes() > 0)
        {
            resList = p->getResistorList();
            res = resList->getHead();
            // Reinitialize for next
                // (xcenter, ycenter, text, bound
            while (res != NULL)
            {
                // Get first resistor's position
                x1 = p->getPosX()+ (NODE_WIDTH/2);
                y1 = res->getPosY();
                // Get second resistor's position
                // Go to next Node,
                traverse = p->getNext();
                resName = res->getLabel();

                //Check if Traverse is not NULL
                if (traverse != NULL)
                {
                    // loop through other nodes until get the same resistor at other places,
                    // if the resistor doesn't exist, it means you looped through it before,
                    // break out of while loop
                    // If traverse has no resistor, this will return false as well
                    while(!(traverse->getResistorPosY(resName,y2)))
                    {
                        traverse = traverse->getNext();
                        if (traverse == NULL) break;
                    }
                }
                // If never looped through this resistor before, y2 wont be 0
                if (y2 != 0)
                {
                    x2 = traverse->getPosX() + (NODE_WIDTH/2);
                    // Draw a red line between the resistors
                    window.gl_setcolor(RED);
                    window.gl_drawline(x1, y1, x2, y2);
                    // Step5: Draw out resistor name and value below the resistor
                    window.gl_setcolor(BLACK);
                    s << resName << " (" << fixed << setprecision(2) << showpoint << res->getResistance() << " Ohm)";
                    ResText = s.str();
                    // Reinitialize for next
                    s.clear(); s.str(string());
                    window.gl_drawtext((x1+x2)/2, y1 - 2, ResText, 250);
                    // Reinitialize y2
                    y2 = 0;
                }
                res = res->getNext();
            }
        }
        // Iterate through Nodes
        p = p->getNext();
        window.gl_setcolor(RED);
    }
    return;
}
