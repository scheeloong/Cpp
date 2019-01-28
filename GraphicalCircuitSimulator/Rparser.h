#ifndef RPARSER_H_INCLUDED
#define RPARSER_H_INCLUDED
#include <iostream> //cout, cin
#include <sstream> // string stream
#include <string> // strings
#include <iomanip> // setw, setprec

using namespace std;

//-----------------------------------------------------------------------------------
// Function Prototypes
//-----------------------------------------------------------------------------------
// Note: Definitions are given above the implementation of the functions below
int parser(string &resName, double &resistance, int *nodeid, double &voltage);
bool checkName(string name);
bool checkGotArguments(stringstream& lineStream);
bool checkNoArguments(stringstream& lineStream);
bool checkType(stringstream& lineStream);
bool checkTypeID(stringstream& lineStream);
bool checkNegative(double resistance);
bool checkNodeIDequals(int nodeid1, int nodeid2);
void removeWhiteSpaces(stringstream& lineStream);
bool checkWhiteSpace(stringstream& lineStream);
//-----------------------------------------------------------------------------------
// Error Messages
//-----------------------------------------------------------------------------------
const string errormsg("Error: ");
const string er1("invalid command");
const string er2("invalid argument");
const string er3("negative resistance");
const string er4("resistor name cannot be the keyword \"all\"");
const string er5("both terminals of resistor connect to node "); // needs to concatenate with arguments
const string er6("too many arguments");
const string er7("too few arguments");
// Extra errors handled in Main.cpp
#endif // RPARSER_H_INCLUDED
