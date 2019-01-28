//-----------------------------------------------------------------------------------
// Main Function
//-----------------------------------------------------------------------------------
// Note: (From Discussion Board)
// 1. Won't be tested with 0 resistance since that will break the provided numerical solution to solve
// 2. Solve error gets outputted as long as : NO (nodes with resistors) have setV = true;

// Code Design:
// 1. This program deletes a Node whenever it is ( not set and has no resistors attached to it)
// 2. This program deletes all Nodes regardless of setV whenever all deleteR all command is run
//    as deleteR command is treated as a "Reset Network" command in this code

#include "Rparser.h" // file tha deals with parsing of commands
#include "NodeList.h"
#include "Node.h"
#include "ResistorList.h"
#include "Resistor.h"
#include "easygl.h" // For WIN32 Graphics
using namespace std;

// Algorithm:
//  1. Parser handles the parsing of commands and if command is valid,
        // it gives cases for main function to handle.
//  2. Main function executes each cases differently.
//  3. Data storage is done using linkedlists NodeList->Node->ResistorList->Resistor

// These additional errors are dealed with in Main and not Rparser.cpp
const string er8("no nodes have their voltage set");
const string er9(" not found");// needs to concatenate with arguments

// switch(action)
// -1; => Do nothing cause invalid command input: implemented as Default
// 0; => Do nothing cause empty string "" given. Also implemented as Default
// 1; => solve
// 2; => setV
// 3; => unsetV
// 4; => insert Resistor
// 5; => modify Resistor
// 6; => print Resistor
// 7; => print node
// 8; => print all nodes
// 9; => delete resistor
// 10; => delete all resistors
// 11; => draw

// Main Function's Prototypes
void parse_draw (NodeList& nodeList);

// Global Variables
// Create a global easygl object called window
easygl window ("Resistor display", WHITE);
NodeList nodeList; // Create a global NodeList called nodeList
float xleft, ybottom, xright, ytop;

int main()
{
    int resCount = 0;
    // Declare NodeList to store list of nodes,
    // each which has a resistor list,
    // which stores list of resistors
    int action = 0; // return value by Rparser
    string resName; // used for printing, modifying and inserting resistor name.
    int nodeid[2];  // used for printing nodes and inserting resistor into nodes
    double oldResistance, resistance; // used for getting resistance values of resistor in insertR and modifyR
    double voltage; // for setV, unsetV, solve


    // if EOF (Windows: ctrl + Z, UNIX: ctrl + D) read, break out of while loop
    while (!cin.eof())
    {
        // Prompt for next command and parse it.
        cout << "> ";
        action = parser(resName, resistance, nodeid, voltage);
        switch (action)
        {
            // case -1 and case 0 is included in default
            //-----------------------------------
            case 1: // solve
            //-----------------------------------
                if (nodeList.solvable())
                {
                    nodeList.solve();
                }
                else
                {
                    cout << errormsg+er8 << endl;
                }
                break;

            //-----------------------------------
            case 2: // setV nodeid voltage
            //-----------------------------------
                nodeList.setVoltage(nodeid[0], voltage);
                cout << "Set: node " << nodeid[0] << " to " << fixed << setprecision(2)
                     << showpoint << voltage << " Volts" << endl;
                break;
            //-----------------------------------
            case 3: // unsetV nodeid
            //-----------------------------------
                nodeList.unsetVoltage(nodeid[0]);
                cout << "Unset: the solver will determine the voltage of node " << nodeid[0] << endl;
                break;

            //----------------------------------------------------
            case 4: // insertR name resistance nodeid1 nodeid2
            //----------------------------------------------------
            // Note: Resistor name does not have to be unique as it is not specifed in lab handout
            // Insert resistor into both nodes it connects to into both nodes in nodelist
                nodeList.insertResistor(resName, resistance, nodeid);
                cout << "Inserted: resistor " << resName << " " << fixed << setprecision(2) << showpoint <<
                resistance << " Ohms " << nodeid[0] << " -> " << nodeid[1] << endl;
                resCount++;
                break;

            //-----------------------------------
            case 5: // modifyR name resistance
            //-----------------------------------
                // Note: It does not handle the case for resistor with same names
                // if that is the case, it will modify the very first 2 resistors connecting to the
                // 2 nodes with the lowest nodeID as it breaks once it finds the first 2 resistors.
                // Check to see if the given resistor name exist
                if (!(nodeList.hasResistor(resName)))
                {
                    cout<< errormsg+"resistor "<< resName << er9 << endl;
                }
                else
                {   // Resistor exists,
                    // Get old resistance first
                    oldResistance = nodeList.getResistance(resName);
                    // Condition below should always evaluate to true
                    if(nodeList.modifyResistor(resName, resistance))
                    {
                        cout << "Modified: resistor " << resName << " from " << oldResistance
                        << " Ohms to " << fixed << setprecision(2) << showpoint << resistance << " Ohms" << endl;
                    }
                }
                break;
            //-----------------------
            case 6: // printR name
            //-----------------------
                if(!(nodeList.printResistor(resName)))
                {
                    // If resistor name given does not exist, print error
                    cout<< errormsg+"resistor "<< resName << er9 << endl;
                }
                break;

            //----------------------------
            case 7: // printNode nodeid
            //----------------------------
            // Note: If the nodeID has no resistor, it does not print.
            // Note: If the nodeID does not exists, it does not print because it has no resistor.
            // Note: Checking if given nodeid is already in node and if it has resistors
            //       is done in NodeList::printNode(); method
                // print all resistors connected to this nodeid => nodeid[0]
                nodeList.printNode(nodeid[0]);
                break;
            //----------------------------
            case 8: // printNode all
            //----------------------------
                    // Print all nodes
                    nodeList.printAllNodes();
                break;

            //----------------------------
            case 9: // deleteR name
            //----------------------------
            // If deleted resistor successfully
                if ((nodeList.deleteResistor(resName)))
                {
                    // Resistor exists, deleted resistor from both nodes it connects t
                    cout << "Deleted: resistor " << resName << endl;
                    resCount --;
                }
            // If resistor not found, print error message
                else
                {
                    cout<< errormsg+"resistor "<< resName << er9 << endl;
                }
                break;

            //----------------------------
            case 10: // deleteR all
            //----------------------------
                // Delete all resistors and
                // Set number of resistors in network to 0
                resCount = 0;
                nodeList.deleteAllResistors();
                // by reinitializing all nodes
                cout << "Deleted: all resistors" << endl;
                break;

            //----------------------------
            case 11: // draw
            //----------------------------
            // Create the window with name My Window and background colour white.
                parse_draw(nodeList);
                break;

            default:
                // do nothing if case 0 or -1
                break;
        }
    }
    // After end of file detected, remove all memory allocated
    nodeList.deleteAllResistors();
    return 0;
}

void parse_draw (NodeList& nodeList)
{
    // Get the maximum and minimum coordinates needed for current network
    nodeList.set_draw_coords(xleft, ybottom, xright, ytop);
    // Set the new world coordinates
    window.set_world_coordinates(xleft, ybottom, xright,ytop);
    // Draw all resistors and nodes
    nodeList.draw(ytop);
    cout << "Draw: control passed to graphics window" << endl;
    // Note: This function will not return until Proceed button is pushed.
    // Used to zoom in and out of pictures
    window.gl_event_loop();
    cout <<"Draw: complete; responding to commands again" << endl ;
}

// Redrawing the screen when graphic packages calls callback
// This function is called whenever the mouse moves over the screen or when it is clicked or anything with
// zooms, clicks or mousemovement
void easygl::drawscreen (void)
{
    window.gl_clearscreen();
    // world coordinates was already set from before for window object
    nodeList.draw(ytop); // use nodeList's draw
}
