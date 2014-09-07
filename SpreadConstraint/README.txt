SpreadConstraint.cpp contains my entire Spread Constraint code and it's propagation algorithm. The main(void) function is
essentially the propagation function.  You can compile it and run it on its own. To implemented this into SCIP, 
all you need to do is copy paste the main() function as the CONS_PROP() function in SCIP, 
read all the domain variables from SCIP, create a new problem using my code, solve it in my code, and then copy paste the solved values
into the domain variables from SCIP. 

SpreadConstraint.m contains the code to calculate the second derivative of the complicated equation, 
as well as graph the dmax(q) to show that it is not continuous, whereas the paper's 
incorrect claim states that dmax(q) is continuous. 

SpreadConstraint.tex contains the LaTeX code to generate the report that I have written 