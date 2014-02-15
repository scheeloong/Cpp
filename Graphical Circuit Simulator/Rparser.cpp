// Soon Chee Loong, 999295793, cheeloong.soon@mail.utoronto.ca
// ECE244Lab4 Fall 2013
// Objective: Creating and Manipulating Linked List in C++ , simulate hardware circuitry
// Note: '\n' = (Pressing Enter), eof = (CTRL + D)
// Note: This program disregards any extra whitespace character trailing all arguments

//-----------------------------------------------------------------------------------
// Header Files
//-----------------------------------------------------------------------------------
#include "Rparser.h"

//-----------------------------------------------------------------------------------
// Main Parser Function
//-----------------------------------------------------------------------------------
// This function takes accepts a line of input from the cin and parse it.
// It returns a different integer value depending on which command input was given
// Parameters: None
// Return Value:
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


int parser(string &resName, double &resistance, int *nodeid, double &voltage)
{
    // Declare strings to store commands and arguments.
    string line, command, name;
    // Reads an entire line from input stream cin into string variable line.
    // Includes trailing white spaces. It stops reading at '\n'
    getline(cin, line);

    //Put the line in as a string in a stringstream for parsing
    // Make a new stringstream for each line so flags etc. are in a known state
    stringstream lineStream(line);
    if(!lineStream.eof())
{
    // Get first string in line
    lineStream >> command;
    //-------------------------------------------------------------------------------
    // Command: solve
    //-------------------------------------------------------------------------------
    if (command == "solve")
    {
        //Check for no extra arguments
        if (!checkGotArguments(lineStream))
        {
            return 1;
        }
    }

    //-------------------------------------------------------------------------------
    // Command: setV nodeid voltage
    //-------------------------------------------------------------------------------
    else if (command == "setV")
    {
        // Check if no more arguments
        if (!checkNoArguments(lineStream))
        {
            // Read first argument
            lineStream >> name;
            stringstream temp; // To check if it's the right type for nodeid
            temp << name;
            temp >> nodeid[0];
            // If it is an integer, execute command
            if(!temp.fail())
            {
                // check if temp has anything left
                // if not end of file for temp
                if(!temp.eof())
                {
                    cout << errormsg+er2 << endl;
                }
                else
                {
                    // check for no arguments
                    if (!checkNoArguments(lineStream))
                        {
                            lineStream >> voltage;
                            // check if voltage is double and for too few arguments
                            if (!checkType(lineStream) && !checkWhiteSpace(lineStream) && !checkGotArguments(lineStream))
                            {
                                //All conditions satisfied for parsing, further checkings done in Main
                                lineStream.ignore(); // delete the rest of the stream in lineStream
                                return 2;
                            }

                        }
                }
            }
            // first argument is not an integer
            else
            {
                cout << errormsg+er2 <<endl;
            }

        }
    }

    //-------------------------------------------------------------------------------
    // Command: unsetV nodeid
    //-------------------------------------------------------------------------------
    else if (command == "unsetV")
    {
        // Check if no more arguments
        if (!checkNoArguments(lineStream))
        {
            // Read first argument
            lineStream >> name;
            stringstream temp; // To check if it's the right type for nodeid
            temp << name;
            temp >> nodeid[0];
            // If it is an integer, execute command
            if(!temp.fail())
            {
                // check if temp has anything left
                // if not end of file for temp
                if(!temp.eof())
                {
                    cout << errormsg+er2 << endl;
                }
                else
                {
                    // Check for too many arguments
                    if (!checkGotArguments(lineStream))
                        {
                            //All conditions satisfied, print output
                            lineStream.ignore(); // delete the rest of the stream in lineStream
                            return 3;
                        }
                }
            }
            // First argument is not an integer
            else
            {
                cout << errormsg+er2 <<endl;
            }
        }
    }
    //-------------------------------------------------------------------------------
    // Command: insertR name resistance nodeid1 nodeid2
    //-------------------------------------------------------------------------------
    else if (command == "insertR")
    {
        // Read first argument
        lineStream >> name;
        // check if name is not all and check for too few arguments
        if (!checkName(name) && !checkNoArguments(lineStream))
        {   // Note: Removed all white space at checkNoArguments
            // Read second argument
            lineStream >> resistance;

            // Cases for error handling here onwards:
            // i) 1alphabet 2                   -> invalid argument
            // ii) 1 2alphabet                  -> invalid argument
            // iii) 1alphabet 2alphabet         -> invalid argument
            // iv) 1 2 alphabet                 -> too many argument
            // v) 1 2alphabet alphabet          -> invalid argument
            // vi) 1alphabet 2alphabet alphabet -> invalid argument
            // viii) 1alphabet 2 alphabet       -> invalid argument
            // All satisfied

            // check if resistance is double and check if resistance is negative and check for too few arguments
            if (!checkType(lineStream) && !checkNegative(resistance) && !checkNoArguments(lineStream))
            {
                // Read third argument
                lineStream >> nodeid[0];

                // check if nodeid1 is int and check for too few arguments
                if (!checkTypeID(lineStream) && !checkNoArguments(lineStream))
                {
                    // Note: If given "1alphabet" as value, "1" will be recorded above,
                    //       and "alphabet" will be recorded below as error
                    // Check if value given was a double,
                    // Read fourth argument
                    // Note: If given "2alphabet" as value, "2" will be recorded below,
                    //       but "alphabet" will be recorded as extra arguments.
                    lineStream >> nodeid[1];
                    // check if nodeid2 is int and
                    // check if final argument was an integer & not an integer followed by characters e.g. 20alphabet and
                    //  check if nodeid1 == nodeid2 and check if too many arguments
                    if (!checkTypeID(lineStream) && !checkWhiteSpace(lineStream) && (!checkNodeIDequals(nodeid[0], nodeid[1])) && !checkGotArguments(lineStream))
                    {
                        //All conditions satisfied, print output
                        resName = name;
                        lineStream.ignore(); // delete the rest of the stream in lineStream
                        return 4;
                    }
                }
            }
        }
    }

    //-------------------------------------------------------------------------------
    // Command: modifyR name resistance
    //-------------------------------------------------------------------------------
    else if (command == "modifyR")
    {
        // Read first argument
        lineStream >> name;
        // check if name is not all and check for too few arguments
        if (!checkName(name) && !checkNoArguments(lineStream))
        {
            // Read second argument
            lineStream >> resistance;
            // check if resistance is double and check if resistance is negative and check for too few arguments
            if (!checkType(lineStream) && !checkWhiteSpace(lineStream) && !checkNegative(resistance) && !checkGotArguments(lineStream))
            {
                //All conditions satisfied, print output
                lineStream.ignore(); // delete the rest of the stream in lineStream
                resName = name;
                return 5;
            }
        }
    }
    //-------------------------------------------------------------------------------
    // Command: printR name || printR all
    //-------------------------------------------------------------------------------
    else if (command == "printR")
    {
        // Read first argument
        if (!checkNoArguments(lineStream))
        {
            lineStream >> name;
            // check if name is not all and check for too few arguments
            if (!checkGotArguments(lineStream))
            {
                // name given is a specific name
                //All conditions satisfied, print output
                lineStream.ignore(); // delete the rest of the stream in lineStream
                resName = name; // assign resName to name
                return 6;
            }
        }
    }

    //-------------------------------------------------------------------------------
    // Command: printNode nodeid || printNode all
    //-------------------------------------------------------------------------------
    else if (command == "printNode")
    {
        if (!checkNoArguments(lineStream))
        {
            // Read first argument
            lineStream >> name;
            stringstream temp; // To check if it's the right type for nodeid
            temp << name;
            temp >> nodeid[0];
            // If it is an integer execute command
            if(!temp.fail())
            {
            // check if temp has anything left
            // if not end of file
                if(!temp.eof())
                {
                    cout << errormsg+er2 << endl;
                }
                else
                {
                    // check for too many arguments
                    if (!checkGotArguments(lineStream))
                        {
                            //All conditions satisfied, print output
                            lineStream.ignore(); // delete the rest of the stream in lineStream
                            return 7;
                        }
                }
            }
            // If it is the argument "all", and not too many arguments execute command
            else if(name == "all")
            {
                if (!checkGotArguments(lineStream))
                {
                    // All conditions satisfied, print output
                    lineStream.ignore(); // delete the rest of the stream in lineStream
                    return 8;
                }
                // Here, there are still arguments, error outputted in function above,
                // do nothing
            }
            else
            {
                cout << errormsg+er2 <<endl;
            }
        }
    }
    //-------------------------------------------------------------------------------
    // Command: deleteR name || deleteR all
    //-------------------------------------------------------------------------------
    else if (command == "deleteR")
    {
        // Check for too few arguments
        if (!checkNoArguments(lineStream))
        {
            // Read first argument
            lineStream >> name;
            // check if name is not all and check for too many arguments
            if (!checkGotArguments(lineStream))
            {
                if (name == "all")
                {
                    //All conditions satisfied, return case 10
                    lineStream.ignore(); // delete the rest of the stream in lineStream
                    return 10;
                }
                else
                {
                    //All conditions satisfied
                    resName = name;
                    lineStream.ignore(); // delete the rest of the stream in lineStream
                    return 9;
                }
            }
        }
    }

    //-------------------------------------------------------------------------------
    // Command: draw
    //-------------------------------------------------------------------------------
    else if (command == "draw")
    {
        //Check for no extra arguments
        if (!checkGotArguments(lineStream))
        {
            return 11;
        }
    }
    //-------------------------------------------------------------------------------
    // Command: Invalid command or no command given
    //-------------------------------------------------------------------------------
    else if (command == "")
    {
	// do nothing
        lineStream.ignore(); // delete the rest of the stream in lineStream
        return 0;
    }

    else
    {
        cout << errormsg+er1 <<endl;
    }
}
    //-------------------------------------------------------------------------------
    // Clean Up and return for next line to be parsed.
    //-------------------------------------------------------------------------------
    // If you are here, error occured in one of the if statements
    lineStream.ignore(); // delete the rest of the stream in lineStream
    return -1;
}

//-----------------------------------------------------------------------------------
// Helper Functions (for Parser Function)
//-----------------------------------------------------------------------------------
// This function checks that the string name is not "all"
// It prints the error message if it is.
// Parameters: string name
// Return Value: false if successful, true if error occurs
bool checkName(string name)
{
    if (name == "all")
    {
        cout << errormsg+er4 << endl;
        return true;
    }
    return false;
}

// This function checks that the stringstream lineStream has reached the end of file
// It prints the error of too many arguments if it hasn't
// Parameters: stringstream& lineStream
// Return Value: false if successful, true if error occurs
bool checkGotArguments(stringstream& lineStream)
{
    removeWhiteSpaces(lineStream);
    if (!lineStream.eof())
    {
        cout << errormsg+er6 << endl;
        return true;
    }
    return false;
}

// This function checks that the stringstream lineStream has not reached the end of file
// It prints the error of too few arguments if it has.
// Parameters: stringstream& lineStream
// Return Value: false if successful, true if error occurs
bool checkNoArguments(stringstream& lineStream)
{
    removeWhiteSpaces(lineStream);
    if (lineStream.eof())
    {
        cout << errormsg+er7 << endl;
        return true;
    }
    return false;
}

// This function checks that the type given from stringstream lineStream is the right type.
// It prints the error of wrong type if it isn't
// Parameters: stringstream& lineStream
// Return Value: false if successful, true if error occurs
bool checkType(stringstream& lineStream)
{
    if (lineStream.fail())
    {
        cout << errormsg+er2 <<endl;
        return true;
    }
    return false;
}

// This function checks that the type given from stringstream lineStream is the integer type.
// It prints the error of wrong type if it isn't
// Parameters: stringstream& lineStream
// Return Value: false if successful, true if error occurs
bool checkTypeID(stringstream& lineStream)
{
    if (checkType(lineStream))
    {
        return true;
    }
    if (lineStream.peek() == '.')
    {
        cout << errormsg+er2 <<endl;
        return true;
    }
    return false;
}

// This function checks that the resistance value given is not negative
// It prints the error of negative resistance if it is.
// Parameters: double resistance
// Return Value: false if successful, true if error occurs
bool checkNegative(double resistance)
{
    if (resistance < 0.0)
    {
        cout << errormsg+er3 << endl;
        return true;
    }
    return false;
}

// This function checks that nodeid1 is not equal to nodeid2
// It prints the error of resistance coneccted to the same node if it is.
// Parameters: int nodeid1, int nodeid2
// Return Value: false if successful, true if error occurs
bool checkNodeIDequals(int nodeid1, int nodeid2)
{
    if (nodeid1 == nodeid2)
    {
        cout << errormsg+er5 << nodeid1 << endl;
        return true;
    }
    return false;
}

// This function removes all trailing whitespaces on a stream
// Parameters: stringstream& lineStream
// Return Value: void
void removeWhiteSpaces(stringstream& lineStream)
{
    while (lineStream.peek() == ' ')
    {
	lineStream.ignore(1,' ');
    }
    return;
}

//-------------------------------------------------------------------------------------------------------
// This function checks that the next part is a white space or end of file at the end of last argument
// If not, an invalid argument must have been given
// It prints the error of invalid argument if it is not a white space and not end of file
// Parameters: stringstream& lineStream
// Return Value: false if successful, true if error occurs
bool checkWhiteSpace(stringstream& lineStream)
{
    // it means a value followed by alphabets were given
   if (lineStream.peek() != ' ' && !lineStream.eof())
    {
        cout << errormsg+er2<<endl;
        return true;
    }
    return false;
}
// Note: Example where checkWhiteSpace is used is if given input is
// "insertR R1 100 1 2abcd" => "abcd" will be recorded as too many arguments if checkWhite space is not used
// need check white space to know that "abcd" is part of "2" and thus it is an invalid argument error
