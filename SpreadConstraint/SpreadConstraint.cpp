/*
Hi, I did not have time to think about these and I realize this today on my last day on research.
So, when testing, if you encounter errors, the errors might be at the below cases.
Important ReadMe:
1. For the general mean value, I pick the q that resulted in the largest
DMAX(q) function. Note: I did not account for the cases where q results in shifting to a previous interval,
and therefore need to calculate q recursively for every q at each intervals end points.
The change can be done easily, but   am not sure if the current code is correct,
or I am suppose to calculate the largest q from DMAX(q) function recursively to select q.

2. There is an error in the code. I prune each X individually before I updated them.
That shouldn't be right as it would end up pruning more than allowed based on given stdDevUB.
The right method should be break out of the current iteration as soon as 1 X changes bounds,
and to re-iterate from the beginning. This happens in my code for both the specific and general mean values case.
You can find it by searching for TODO
All the TODO comments are things that needs to be done that I haven't completed.

3. There is a part of the code under TODO: that I did not account for m being equal to 0. Currently,
it does not cause any problems from my test cases, but if you start seeing infinite, nan or negative
values in the bounds, this should be the problem.

4. Need fix the inf and nan values problem from finding qLast intelligently. Then,
everything from DMax should be correct.

5. Once done 4., use the same solution to update DMin as well. Shouldn't be too much work since
now instead of just comparing the end points, you compare both endpoints and the single localmax for DMIN.
*/

// This file contains my implementation of Spread Constraint
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <algorithm>
#include <math.h>

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdarg.h>  // to accept variable arguments
using namespace std;

// Note: The paper does not account for the special case where
// i) m = 0, if m = 0, lots of equations divides by m, division by 0 is undefined
// ii) m = n, if m = n, some equatons divides by (n-m), division by 0 is undefined

// m == n is used to calculate q0 for opt(q), which is needed to prune U using S.
// However, can be proven if m == n, the q0 won't be in that interval, so no worries.

// m == 0  is used to calculate DMax(q)
// However, can note that if m = 0, we don't need to calculate v as no X values will be assigned to v
// so we re-calculate everything by ignoring v, in other words, set v = 0.
// this results in

// if (m == 0)
//{
//    a = 1.0;
//    v = 0;
//    opt doesnt add the ((q-ES)^2/m) term
//}

// else (normal case)
// a = 1.0 + 1.0/m
// v = (q-ES)/m
// opt adds the ((q-ES)^2/m) term

// Ignore all
class Interval
{
    friend class domainVar; // allows domainVar to easily access all public and private members of this class
public: // to make life easier, just call these public,
        // so that main function that contains propagation algorithm can access these easily.
    double Ilb; // interval lower bound
    double Iub; // interval upper bound
    double Vlb; // q's lower bound in this interval to locate v
    double Vub; // q's upper bound in this interval to locate v
    double ES; // extreme sum for all X in this interval, i.e. min + max
    double C; // extreme^2 sum for all X in this interval, i.e. (min)^2 + (max)^2
    int *M; // array of indexes of all X in the M(Iq)
    int *R; // array of indexes of all X in the R(Iq)
    int *L; // array of indexes of all X in the L(Iq)
    int m; // cardinality of M(Iq)
    int r; // cardinality of R(Iq)
    int l; // cardinality of l(Iq)

    // The range of nVar that can be in this interval
     // Note: nVarUb may be lower than nVarLb
    // It may also not contain the minimum value if it contains the global minimum since nVar is a convex function
    double nVarUb;
    double nVarLb;

    double q0; // the value of q0 if this interval was the global min
                // -1 => m == n for this interval

public:
    Interval()
    {
        Ilb = 0;
        Iub = 0;
    }

    // Constructer, initializes with the upper and lower bound of this class
    Interval(int lb, int ub)
    {
        this->Ilb = lb;
        this->Iub = ub;
    }
};

// This class represents a domainVariable
// It includes additional private data for:
// i) Spread Constraint
class domainVar
{
private:
// General for Every Constraint
    int n; // number of X variables
    int **X; // a 2D array of X's and its upper and lower bound
    char name = 'a';

// mergesort for create interval, returns number of values currently in array
// note: This merge sort is unique as it removes duplicate elements while sorting
template <class X>
int mergesort(X a[], int n)
{
    if (n==1)
    {
        return 1;
    }
    int q, p;
    q = n/2;
    p = n/2 + 1;
    X b[q];
    X c[p];
    int i = 0;
    for (i = 0; i < q; i++)
    {
        b[i] = a[i];
    }
    int k = 0;
    for (int j = q; j < n; j++)
    {
        c[k] = a[j];
        k++;
    }
    q = mergesort(b, i);
    p = mergesort(c, k);
    int r, s, t;
    t = 0; r = 0; s = 0;
    while( (r!=q) && (s!= p))
    {
        if (b[r] < c[s])
        {
            a[t] = b[r];
            r++;
        }
        else if (b[r] > c[s])
        {
            a[t] = c[s];
            s++;
        }
        else // here it means both of them are equal
        {
            a[t] = b[r];
            r++;
            s++; // also incrememnt s to skip it
        }
        t++;

    }
    if ( r==q)
    {
        while(s!=p)
        {
            a[t] = c[s];
            s++;
            t++;
        }
    }
    else
    {
        while(r != q)
        {
            a[t] = b[r];
            r++;
            t++;
        }
    }
    return t;
}

public:
// Specially For Spread Constraint
// Let these variables be public so main function can access these easily.
    int numI; // number of contiguous intervals
    Interval **I; // array of contiguous interval bounds (note: This is a 1D array, but simple hack works with 2 used as 1)
    bool needDelete; // to determine if this is the first time interval I is being created, or needs to be clean up
                    // before creating again

    // constructor
    domainVar()
    {
        n = 10;
        X = new int* [n];
        int i = 0;
        for (i = 0; i < n; i++)
        {
            X[i] = new int [2];
            X[i][0] = 0;
            X[i][1] = 0;
        }
        needDelete = false;
    }

    domainVar(int n_args, char symbol)
    {
        this->n = n_args;
        name = symbol;
        X = new int* [n_args];
        int i = 0;
        for (i = 0; i < n_args; i++)
        {
            X[i] = new int [2];
            X[i][0] = 0;
            X[i][1] = 0;
        }
        needDelete = false;
    }

    void updateDomainVar(int index, int lb, int ub)
    {
        this->setDomVarLb(index, lb);
        this->setDomVarUb(index, ub);
    }

    domainVar(char symbol, int n_args, ...)
    {
        name = symbol;
        va_list head;
        va_start(head, n_args); // initialize head to be every variable after n_args
        n = n_args;
        X = new int* [n];

        int lb, ub;
        for(int i = 0; i < n; i++)
        {
            X[i] = new int [2];
            lb = va_arg(head, int);
            ub = va_arg(head, int);
            X[i][0] = lb;
            X[i][1] = ub;
        }
        va_end(head);
        needDelete = false;
    }

    // Public Methods
    void printDomVar()
    {
        int i = 0;
        for(i = 0; i < n; i++)
        {
            cout << "Bounds of "<< name << "[ " << i << "] : " <<  X[i][0] <<" , " << X[i][1] <<endl;
        }
    }

    void setDomVarLb(int i, int value)
    {
        X[i][0] = value;
    }

    int getDomVarLb(int i)
    {
        return X[i][0];
    }

    void setDomVarUb(int i, int value)
    {
        X[i][1] = value;
    }

    int getDomVarUb(int i)
    {
        return X[i][1];
    }

    void setDomN(int value)
    {
        n = value;
    }

    int getDomN()
    {
        return n;
    }
    double maxStdDevThiago()
    {
        // Maximum std. deviation possible from all values of X
        double maxScalc = 0;
        double sumXi = 0;
        double sumXj = 0;

        double Sub = 0; // Sum of upper bounds
        double Slb = 0; // Sum of lower bounds
        int i = 0;
        for (i = 0; i < this->getDomN(); i++ )
        {
            Slb += this->getDomVarLb(i);
            Sub += this->getDomVarUb(i);
        }

        double minUcalc = Slb/this->getDomN();
        double maxUcalc = Sub/this->getDomN();

   double min_mean_;		//
   double max_mean_;		//
   double distance_x_max;	//
   double distance_x_min;	//
    double X[this->getDomN()];
   // Assign each variable, x[k] to its extrema, such that the std dev is maximized
   // Assigning it to the extrema will maximize the std. dev as explained in the paper
   // "pg 13 of 15  Simplification and extension of Spread Constraint by Pierre et al"
   for (int k=0;k < this->n;k++)
   {
	  // assume x = xmax => min_mean increases
      min_mean_ = minUcalc + (this->getDomVarUb(k) - this->getDomVarLb(k))/this->n ;

// As explained in the paper, find out which values to use before assigning values of x in order to get max std. deviation.
		// lower bound of x = xmax
        if ( fabs(this->getDomVarUb(k) - min_mean_) <= fabs(this->getDomVarUb(k) - maxUcalc) )
        {
            distance_x_max = fabs(this->getDomVarUb(k) - min_mean_);
        }
        else
        {
            distance_x_max = fabs(this->getDomVarUb(k) - maxUcalc);
        }

      //assume x = xmin => max_mean decreases
      max_mean_ = maxUcalc - (this->getDomVarUb(k) - this->getDomVarLb(k))/this->n;
	  // upper bound of x = xmin
	  if (fabs(this->getDomVarLb(k) - max_mean_) >= fabs(this->getDomVarLb(k) - minUcalc))
      {
          distance_x_min = fabs(this->getDomVarLb(k) - max_mean_);
      }
      else
      {
          distance_x_min =  fabs(this->getDomVarLb(k) - minUcalc);
      }

      if (distance_x_max > distance_x_min)
      {
			X[k] = this->getDomVarUb(k); // assign to upper bound extreme
			minUcalc = min_mean_; // update the minimum mean
      }

      else // LINE HAHA THIS SHOULDN'T BE CORRECT
      {
         X[k] = this->getDomVarLb(k); // assign to lower bound extreme
		 maxUcalc = max_mean_; // update the maximum mean
      }
   }
   double std_dev = 0;
   // Here, we're done assigning X to its appropriate extreme values,
    // calculate the maximum std. deviatio that can occur
   for(int i = 0; i < this->n; ++i )
      std_dev += pow(X[i] - 3, 2.0);  // numerator of variance calculation
   std_dev = std_dev/this->n; // variance
   std_dev = sqrt(std_dev); // std_deviation
   return std_dev;
    }

    // This algorithm from Wen Yang is faster but not as tight as the general algorithm given on paper
    double maxStdDevSpecific(double mean)
    {
        double Xi[this->getDomN()];
        double maxScalc = 0;
        double diffOne = 0;
        double diffTwo = 0;
        double sumXi = 0;
        // Assign all Xi
        for (int i = 0; i < this->getDomN(); i++)
        {
            diffOne = fabs(this->getDomVarLb(i) - mean);
            diffTwo = fabs(this->getDomVarUb(i) - mean);
      //      cout <<" diffONe is : " << diffOne<< " and DiffTWO is : " << diffTwo << endl;
            if (diffOne >= diffTwo)
            {
                Xi[i] = this->getDomVarLb(i);
            }
            else
            {
                Xi[i] = this->getDomVarUb(i);
            }
     //       cout << " X[ " << i << " ] is " << Xi[i] << endl;
            sumXi += fabs(Xi[i] - mean);
        }
        maxScalc = sumXi * sumXi;
        maxScalc = maxScalc/this->getDomN();
        maxScalc = sqrt(maxScalc);
        return maxScalc;
    }
    // Calculate and returns max Std Dev from current values of X
    double maxStdDev()
    {
        // Maximum std. deviation possible from all values of X
        double maxScalc = 0;
        double sumXi = 0;
        double sumXj = 0;

        double Sub = 0; // Sum of upper bounds
        double Slb = 0; // Sum of lower bounds
        int i = 0;
        for (i = 0; i < this->getDomN(); i++ )
        {
            Slb += this->getDomVarLb(i);
            Sub += this->getDomVarUb(i);
        }

        double minUcalc = Slb/this->getDomN();
        double maxUcalc = Sub/this->getDomN();

        // To keep track of the terms to substitute in when calculating maxScalc
        double Xi[this->getDomN()];
        double Xj[this->getDomN()];
        double minUcalc_ , maxUcalc_; // _ => The updated values for changing the mean and max to other extreme
        double distXmax[2]; // Lower and Upper bound for distance when X = Xmax
        double distXmin[2]; // Lower and Upper bound for distance when X = Xmin

        for (i = 0; i < this->getDomN(); i++)
        {
            minUcalc_ = minUcalc + ((this->getDomVarUb(i) - this->getDomVarLb(i))/(this->getDomN()));
            maxUcalc_ = maxUcalc - ((this->getDomVarUb(i) - this->getDomVarLb(i))/(this->getDomN()));
            distXmin[0] = std::min(fabs(this->getDomVarLb(i) - minUcalc), fabs(this->getDomVarLb(i) - maxUcalc_));
            distXmin[1] = std::max(fabs(this->getDomVarLb(i) - minUcalc), fabs(this->getDomVarLb(i) - maxUcalc_));
            distXmax[0] = std::min(fabs(this->getDomVarUb(i) - maxUcalc), fabs(this->getDomVarUb(i) - minUcalc_));
            distXmax[1] = std::max(fabs(this->getDomVarUb(i) - maxUcalc), fabs(this->getDomVarUb(i) - minUcalc_));

            if (distXmin[0] >= distXmax[1])
            {
                Xi[i] = this->getDomVarLb(i);
                Xj[i] = this->getDomVarLb(i);
            }

            else if (distXmax[0] >= distXmin[1])
            {
                Xi[i] = this->getDomVarUb(i);
                Xj[i] = this->getDomVarUb(i);
            }

            else if (this->getDomVarUb(i) < minUcalc )
            {
                Xi[i] = this->getDomVarLb(i);
                Xj[i] = this->getDomVarLb(i);
            }

            else if (this->getDomVarLb(i) > maxUcalc)
            {
                Xi[i] = this->getDomVarUb(i);
                Xj[i] = this->getDomVarUb(i);
            }

            else // you assign them to get maximum deviation possible if nothing can be inferred
            {
               Xi[i] = this->getDomVarUb(i);
               Xj[i] = this->getDomVarLb(i);
            }
            sumXi += Xi[i] * Xi[i];
            sumXj += Xj[i];

        }
        sumXj = (sumXj * sumXj)/ this->getDomN();
        maxScalc = sumXi - sumXj; // n*variance
        maxScalc /= this->getDomN(); // variance
        maxScalc = sqrt(maxScalc); // std.deviation
        return maxScalc;
    }

    // This function creates the interval and it's values based on the current values of X
    void createInterval()
    {
        if(this->needDelete)
        {
//       cout << "Deleting previous intervals" << endl;
            for (int i = 0; i < this->numI; i++)
            {
                delete [] this->I[i]->M;
                delete [] this->I[i]->L;
                delete [] this->I[i]->R;
                delete this->I[i];
            }
            delete [] this->I;
        }
        int B[2*n]; // create an array to hold all current values of X
        for (int i = 0; i < n; i++)
        {
            B[i] = this->getDomVarLb(i);
            B[i+n] = this->getDomVarUb(i);
        }
        this->numI = this->mergesort(B, 2*n) - 1; // if B has 5 elements, this means you only have 4 intervals
        // Since now have all the unique domain variables, can create variable I
        this->I = new Interval*[this->numI];

        for (int i = 0; i < this->numI; i++)
        {
            // Set the upper and lower bound
            this->I[i] = new Interval(B[i], B[i+1]);
            // Now set values for M(Iq), R(Iq), L(Iq) as well as their cardinalities: m, r, l
            // While calculating ES
            this->I[i]->M = new int[this->getDomN()];
            this->I[i]->L = new int[this->getDomN()];
            this->I[i]->R = new int[this->getDomN()];
            this->I[i]->m = 0;
            this->I[i]->l = 0;
            this->I[i]->r = 0;
            this->I[i]->ES = 0;
            this->I[i]->C = 0;

            for (int j = 0; j <  this->getDomN(); j++)
            {
                if ( this->getDomVarUb(j) <= this->I[i]->Ilb )
                {
                    this->I[i]->L[this->I[i]->l] = j; // store the index of x which belongs to L(Iq)
                    this->I[i]->l++;
                    this->I[i]->ES += this->getDomVarUb(j);
                    this->I[i]->C += (this->getDomVarUb(j) * this->getDomVarUb(j));
                }
                else if (this->getDomVarLb(j) >= this->I[i]->Iub)
                {
                    this->I[i]->R[this->I[i]->r] = j; // store the index of x which belongs to R(Iq)
                    this->I[i]->r++;
                    this->I[i]->ES += this->getDomVarLb(j);
                    this->I[i]->C += (this->getDomVarLb(j) * this->getDomVarLb(j));
                }
                else
                {
                    this->I[i]->M[this->I[i]->m] = j; // store the index of x which belongs to M(Iq)
                    this->I[i]->m++;
                }
            }
            // Calculate the Vmin and Vmax for this interval
            // Vmin = ES + Imin*m
            this->I[i]->Vlb = this->I[i]->ES + (this->I[i]->Ilb * this->I[i]->m);
            this->I[i]->Vub = this->I[i]->ES + (this->I[i]->Iub * this->I[i]->m);

            // Calculate the range of n*variance that can be in this interval
            if (this->I[i]->m == 0)
            {
                // If m is 0, there means there are no Middle intervals
                this->I[i]->nVarLb = this->I[i]->C - ((this->I[i]->Vlb * this->I[i]->Vlb)/this->getDomN());
                this->I[i]->nVarUb = this->I[i]->C - ((this->I[i]->Vub * this->I[i]->Vub)/this->getDomN());
            }
            else
            {
                this->I[i]->nVarLb = this->I[i]->C + (((this->I[i]->Vlb - this->I[i]->ES)*(this->I[i]->Vlb - this->I[i]->ES))/this->I[i]->m) - ((this->I[i]->Vlb * this->I[i]->Vlb)/this->getDomN());
                this->I[i]->nVarUb = this->I[i]->C + (((this->I[i]->Vub - this->I[i]->ES)*(this->I[i]->Vub - this->I[i]->ES))/this->I[i]->m) - ((this->I[i]->Vub * this->I[i]->Vub)/this->getDomN());
            }

            // If special case where 2nd derivative is 0
            if ( this->I[i]->m == this->getDomN())
            {
                // do nothing, cause we're guaranteed the global min does not belong to this interval
                this->I[i]->q0 = -1; // m == n
            }
            else
            {
                // calculate the global min q if it was this interval
                this->I[i]->q0 = ((this->I[i]->ES * this->n)/(this->n  - this->I[i]->m)) ;
            }
            // Print only the first time
            if (this->needDelete == false)
            {
                // Print this to check interval initialize properly
                cout << "Interval: " << i << " Ilb: " << this->I[i]->Ilb << " Iub: " <<
                this->I[i]->Iub <<  " ES: " << this->I[i]->ES << " C: " << this->I[i]->C
                << " Vlb: " <<  this->I[i]->Vlb << " Vub: " << this->I[i]->Vub << endl
                << " nVarLb: " << this->I[i]->nVarLb << " nVarUb: " <<
                this->I[i]->nVarUb << " q0: " << this->I[i]->q0 << endl;
                cout << "m: " << this->I[i]->m << " l: " << this->I[i]->l << " r: " << this->I[i]->r <<endl;
            }
        }
        this->needDelete = true;
        return;
    }

    // Returns N*Variance value from the OPT graph given a value of q
    // can only be called after createInterval() function is used at least once
    // Note: It returns 0 for every q value outside the given domains of the intervals
    double OPT(double q)
    {
        double nVariance = 0;
        // If q is not within the intervals, return 0
        if((q < this->I[0]->Vlb) || (q > this->I[this->numI -1]->Vub))
        {
            return 0;
        }
        // Figure out which interval q is in
        int qIndex = 0;
        for (int i = 0; i < this->numI; i++)
        {
            if((q >= this->I[i]->Vlb) && (q <= this->I[i]->Vub))
            {
                qIndex = i;
                break;
            }
        }
        if (((double)this->I[qIndex]->m) == 0)
        {
            nVariance = (this->I[qIndex]->C - (pow((double) q, 2.0)/((double)this->getDomN()))) ;
        }
        else
        {
            nVariance = (this->I[qIndex]->C + ((pow((double)(q - this->I[qIndex]->ES), 2.0))/((double)this->I[qIndex]->m)) - (pow((double) q, 2.0)/((double)this->getDomN()))) ;
        }
        return nVariance;
    }

    // Note: optMax is from the SpreadConsS's upper bound
    // get the variance from S and multiply by N
    // bound = 0 => lower bound of interval
    // bound = 1 => upper bound of interval
    // Need bound cause for discontinuous graph,
    // bound can be upper or lower bound
    double DMAX(double q, double XMIN, double optMAX, int bound)
    {
        double dmax = 0;
        // If q is not within the intervals, return 0
        if((q < this->I[0]->Vlb) || (q > this->I[this->numI -1]->Vub))
        {
            return 0;
        }
        // Figure out which interval q is in
        int qIndex = 0;
        // If lower bound, start from lowest interval and look for it
        if (bound ==  0)
        {
            for (int i = 0; i < this->numI; i++)
            {
                if((q >= this->I[i]->Vlb) && (q <= this->I[i]->Vub))
                {
                    qIndex = i;
                    break;
                }
            }
        }
        // If upper bound, start from highest interval and look for it
        else // if (bound == 1)
        {
            for (int i = numI - 1; i >= 0; i--)
            {
                if((q >= this->I[i]->Vlb) && (q <= this->I[i]->Vub))
                {
                    qIndex = i;
                    break;
                }
            }
        }
        double v = 0;
        double m = 1.0*this->I[qIndex]->m;
        double a = 0;

        if (m == 0)
        {
            a = 1.0; // since the 1\m term will not exist if v is not involved
            v = 0; // v value is not needed since no X carries the value v
        }
        else
        {
            v = ((q - this->I[qIndex]->ES)/((double)m));
            a = 1.0 + (1.0/((double)m));
        }
        double n = this->getDomN();
        double opt = this->OPT(q); // takes into account if m is 0
        double b = XMIN - v;
        double c = opt - n*(pow((double) optMAX, 2.0));
        dmax = ((-b + (sqrt((b*b) - (a*c))))/((double)a));
        return dmax;
    }

    double DMIN(double q, double XMAX, double optMAX, int bound)
    {
        double dmin = 0;
        // If q is not within the intervals, return 0
        if((q < this->I[0]->Vlb) || (q > this->I[this->numI -1]->Vub))
        {
            return 0;
        }
        // Figure out which interval q is in
        int qIndex = 0;
        // If lower bound, start from lowest interval and look for it
        if (bound ==  0)
        {
            for (int i = 0; i < this->numI; i++)
            {
                if((q >= this->I[i]->Vlb) && (q <= this->I[i]->Vub))
                {
                    qIndex = i;
                    break;
                }
            }
        }
        // If upper bound, start from highest interval and look for it
        else // if (bound == 1)
        {
            for (int i = numI - 1; i >= 0; i--)
            {
                if((q >= this->I[i]->Vlb) && (q <= this->I[i]->Vub))
                {
                    qIndex = i;
                    break;
                }
            }
        }
        double m = 1.0*this->I[qIndex]->m;
        double v = 0;
        double a = 0;
        if (m == 0)
        {
            a  = 1.0;
            v  = 0;
        }
        else
        {
            v = ((q - this->I[qIndex]->ES)/((double)m));
            a = 1.0 + (1.0/((double)m));
        }
        double n = this->getDomN();
        double opt = this->OPT(q); // handles case where m == 0

        double b = v - XMAX; // note: This is different from DMAX
        double c = opt - n*(pow((double) optMAX, 2.0));
        dmin = ((-b + (sqrt((b*b) - (a*c))))/((double)a));
        return dmin;
    }

    double findDMax(double qSelected,int indexQSelected, int xIndex, int stdDevUb)
    {
   //     cout << " Parameters are: " << endl
   //     << "qSelected: " << qSelected << " indexQSelected: "
   //      << indexQSelected << " xIndex: " << xIndex << " stdDevUb: " << stdDevUb << endl;
    //    cout << " Calculated nVariance from stdDevUB: " << (pow((double) stdDevUb, 2.0)*this->n) <<endl;
        // Initialize d11 to maximum shift
        double d11 =  qSelected - this->I[indexQSelected]->Vlb;
        double dMax = 0;
        double a = 0;
        double b = 0;
        double c = 0;
        double v = 0;
        double nVarCurrent = 0;

        if (this->I[indexQSelected]->m == 0 )
        {
            v = 0;
            a = 1.0;
        }
        else
        {
            v = (qSelected - this->I[indexQSelected]->ES)/((double) this->I[indexQSelected]->m);
            nVarCurrent += (pow((double)(v - qSelected/((double)this->n)), 2.0) * this->I[indexQSelected]->m);
            a = 1.0+ 1.0/this->I[indexQSelected]->m;

        }

        int index = 0;
        for (int i = 0; i < this->I[indexQSelected]->l; i++)
        {
            index = this->I[indexQSelected]->L[i];
            nVarCurrent += pow((double)(this->getDomVarUb(index) - qSelected/((double)this->n)), 2.0);
        }
        for (int i = 0; i < this->I[indexQSelected]->r; i++)
        {
            index = this->I[indexQSelected]->R[i];
            nVarCurrent += pow((double)(this->getDomVarLb(index) - qSelected/((double)(this->n))),2.0) ;
        }

        b = this->getDomVarLb(xIndex) - v;
        c = nVarCurrent - (pow((double) stdDevUb, 2.0)*this->n);

    //    cout << " A: " << a << " B: " << b << " C: " << c <<" nVarCurrent " <<nVarCurrent << endl ;

        dMax = (-b + sqrt(pow (b, 2.0) - a*c))/a;
//cout << "dMax: " << dMax << " d11: " << d11 << endl;
        if (dMax < d11)
        {
            return dMax;
        }
        else
        {
            if((indexQSelected == 0) || (d11 == 0))
                return d11;
            else
            {
                // Update X's minimum value
                this->setDomVarLb(xIndex, getDomVarLb(xIndex) + floor(d11));
                this->createInterval(); // create the new intervals
                // Determine which indexQ is in
                for (int i = 0 ; i < this->numI ; i++)
                {
                    if ((qSelected >= this->I[i]->Vlb) && (qSelected <= this->I[i]->Vub))
                    {
                        indexQSelected = i;
                        break;
                    }
                }
                return (d11 + this->findDMax(qSelected, indexQSelected, xIndex, stdDevUb));
            }
        }
        // SHOULD NEVER COME HERE
        cout << "ERROR: WENT TO return dMAX" << endl;
        return dMax;
    }

    double findDMin(double qSelected,int indexQSelected, int xIndex, int stdDevUb)
    {
   //     cout << " Parameters are: " << endl
  //      << "qSelected: " << qSelected << " indexQSelected: "
    //     << indexQSelected << " xIndex: " << xIndex << " stdDevUb: " << stdDevUb << endl;
   //     cout << " Calculated nVariance from stdDevUB: " << (pow((double) stdDevUb, 2.0)*this->n) <<endl;

        // Initialize d11 to maximum shift
        double d11 =   this->I[indexQSelected]->Vub - qSelected;
        double dMin = 0;
        double a = 0;
        double b = 0;
        double c = 0;
        double v = 0;
        double nVarCurrent = 0;
        if(this->I[indexQSelected]->m == 0)
        {
            v = 0;
            a = 1.0;
        }
        else
        {
            v = (qSelected - this->I[indexQSelected]->ES)/((double) this->I[indexQSelected]->m);
            nVarCurrent += (pow((double)(v - qSelected/((double)this->n)), 2.0) * this->I[indexQSelected]->m);
            a = 1.0+ 1.0/this->I[indexQSelected]->m;

        }

        int index = 0;
        for (int i = 0; i < this->I[indexQSelected]->l; i++)
        {
            index = this->I[indexQSelected]->L[i];
            nVarCurrent += pow((double)(this->getDomVarUb(index) - qSelected/((double)this->n)), 2.0);
        }
        for (int i = 0; i < this->I[indexQSelected]->r; i++)
        {
            index = this->I[indexQSelected]->R[i];
            nVarCurrent += pow((double)(this->getDomVarLb(index) - qSelected/((double)(this->n))),2.0) ;
        }

        b = v - this->getDomVarUb(xIndex) ; // note: This is different from DMax
        c = nVarCurrent - (pow((double) stdDevUb, 2.0)*this->n);

    //    cout << " A: " << a << " B: " << b << " C: " << c <<" nVarCurrent " <<nVarCurrent << endl ;

        dMin = (-b + sqrt(pow (b, 2.0) - a*c))/a;
        cout << "dMin: " << dMin << " d11: " << d11 << endl;
        if (dMin < d11)
        {
            return dMin;
        }
        else
        {
            // return if last interval
            if((indexQSelected == (this->numI-1)) || (d11 == 0))
                return d11;
            else
            {
                // Update X's maximum value
                this->setDomVarUb(xIndex, getDomVarLb(xIndex) - floor(d11));
                this->createInterval(); // create the new intervals
                // Determine which indexQ is in
                for (int i = 0 ; i < this->numI ; i++)
                {
                    if ((qSelected >= this->I[i]->Vlb) && (qSelected <= this->I[i]->Vub))
                    {
                        indexQSelected = i;
                        break;
                    }
                }
                return (d11 + this->findDMin(qSelected, indexQSelected, xIndex, stdDevUb));
            }
        }
        // SHOULD NEVER COME HERE
        cout << "ERROR: WENT TO return dMin" << endl;
        return dMin;
    }

};

int main(void)
{
 double minSTDDEV = 99999;
 double tempMean = 0;
 double haha = 0;
double tempX1 = 0;
double tempX2 = 0;
double tempX3 = 0;
    // Temporary to brute force to calculate optimal answer
    for(int i = 1; i < 4; i++)
    {
        for(int j = 2; j<7; j++)
        {
            for(int k = 3; k<9; k++)
            {
                tempMean = (i+j+k)/3.0;
                if(tempMean == 5)
                {
                    haha = ((i - tempMean)*(i - tempMean) +  (j - tempMean)*(j - tempMean) +  (k - tempMean)*(k - tempMean));
                    haha /= 3.0;
                    haha = sqrt(haha);
                    if (haha < minSTDDEV)
                    {
                        minSTDDEV = haha;
                        tempX1 = i;
                        tempX2 = j;
                        tempX3 = k;
                    }
                }
            }
        }
    }
    cout << " Optimal value for X1: " << tempX1 << " X2: " << tempX2 << " X3: " << tempX3 << endl;

    // TODO Reverse Engineer: Brute Force find optimal solution for X1, X2, X3 without U and S
    // Find the minimum S able to be obtained,
    // and check if propagation prunes the optimal solution.

    cout << "Hello world!" << endl;
    int i = 0;

// Currently, 12th August 2014,
// my propagation algorithm appears to be correct (from brute force calculations)
// and prunes tighter than CP Optimizer for both general and specific mean case


// Assuming read N, lb and ub values from SCIP
    int N = 3; int lb[3] ={1, 2, 3}; int ub[3] = {3,6,9};
    domainVar SpreadConsX(N, 'X'); // create space for 3 domain variables
    for (int i = 0; i < N; i++)
    {
        SpreadConsX.updateDomainVar(i, lb[i], ub[i]);
    }
    SpreadConsX.printDomVar();
    domainVar SpreadConsU('U', 1, 3, 3);
    SpreadConsU.printDomVar();
    domainVar SpreadConsS('S', 1, 1, 10);
    SpreadConsS.printDomVar();
    // Create the intervals
    SpreadConsX.createInterval();

/*
    domainVar SpreadConsX('X', 3, 10, 20, 12, 20, 3, 30);
    SpreadConsX.printDomVar();
    domainVar SpreadConsU('U', 1, 1, 9);
    SpreadConsU.printDomVar();
    domainVar SpreadConsS('S', 1, 1, 30);
    SpreadConsS.printDomVar();
    // Create the intervals
    SpreadConsX.createInterval();
*/

/*
    domainVar SpreadConsX('X', 3, 1, 3, 2, 6, 3, 9);
    SpreadConsX.printDomVar();
    domainVar SpreadConsU('U', 1, 2, 2);
    SpreadConsU.printDomVar();
    domainVar SpreadConsS('S', 1, 1, 10);
    SpreadConsS.printDomVar();
    // Create the intervals
    SpreadConsX.createInterval();
*/
/*
// Tested and appears to be correct
// For debugging OPT
    for (int yp = 0; yp < 19; yp++)
    {
        double hehe = SpreadConsX.OPT(yp);
        cout << "q: " << yp << " OPT: " << hehe << endl;
    }
*/
/*
// Tested and appears to be correct
// For debugging DMAX Graph
    for (int hp = 0; hp < 21; hp++)
    {
        // check lower bound
        double huhuLB = SpreadConsX.DMAX(hp, 1.0, 50, 0);
        cout << "q: " << hp << " OPTlb: " << huhuLB << endl;
        double huhuUB = SpreadConsX.DMAX(hp, 1.0, 50, 1);
        cout << "q: " << hp << " OPTub: " << huhuUB << endl;
    }
*/

/*
    domainVar SpreadConsX('X', 3, 1, 3, 2, 6, 3, 9);
    SpreadConsX.printDomVar();
    domainVar SpreadConsU('U', 1, 5,5);
    SpreadConsU.printDomVar();
    domainVar SpreadConsS('S', 1, 1, 100);
    SpreadConsS.printDomVar();
    // Create the intervals
    SpreadConsX.createInterval();
*/

    // Note: Found out need to prune both ways for M(Iq)
    // Note: Found out can only change bounds after pruning all X, so that don't result in
    // QIndex not existing segmentation fault.
    // Also found out that adding the second way to prune M(Iq) can result in
    // S being (0,0) as it prunes further and when that happens,
    // I cannot run any X propagation anymore, so break right away
    // Note: Found out that forgot to handle special case where m = 0,
    // if m = 0, cannot divide by m!! thus, that's why previous code got lazy and just added 0.000001
    // but the right thing is to handle everything where m = 0 => need to recalculate a simpler equation
    // which was done and now it works.

/*
    domainVar SpreadConsX('X', 3, 1, 3, 2, 6, 3, 9);
    SpreadConsX.printDomVar();
    domainVar SpreadConsU('U', 1, 6, 6);
    SpreadConsU.printDomVar();
    domainVar SpreadConsS('S', 1, 1, 4);
    SpreadConsS.printDomVar();
    // Creation of intervals is where R, M, L is determined for each of the original intervals
    // Note: Creation of intervals only depends on X and not S or U
    // Create the intervals
    SpreadConsX.createInterval();
*/

    // Propagation Algorithm
    bool changed = true;
    // While something is pruned , keep running propagation
    while (changed)
    {
        // Initialize to false
        changed = false;

        //------------------------------------------------------------
        // Prune the Std Deviation Maximum value using X
        double maxScalc  = 0;
        if ( SpreadConsU.getDomVarLb(0) == SpreadConsU.getDomVarUb(0))
        {
            // Specifc case, faster but not as tight
            maxScalc = SpreadConsX.maxStdDevSpecific( (double) SpreadConsU.getDomVarUb(0));
            // General case, prunes tighter
            maxScalc = SpreadConsX.maxStdDev();
        }
        // General Case
        else
        {
            maxScalc = SpreadConsX.maxStdDev();
           //  maxScalc = SpreadConsX.maxStdDevThiago();// Note: Has been confirmed to be wrong
        }
        if (maxScalc < SpreadConsS.getDomVarUb(0))
        {
           // changed = true; // Note: Dont need to re-run since this is the first to get propagated
            SpreadConsS.setDomVarUb(0,floor(maxScalc));
            if (SpreadConsS.getDomVarUb(0) < SpreadConsS.getDomVarLb(0))
            {
                SpreadConsS.setDomVarLb(0,SpreadConsS.getDomVarUb(0));
            }
        }
        cout << "Propagate S using X" <<endl;
        SpreadConsS.printDomVar();
        //------------------------------------------------------------
        // Prune the Mean bounds using X
        double Sub = 0; // Sum of upper bounds
        double Slb = 0; // Sum of lower bounds

        for (i = 0; i < SpreadConsX.getDomN(); i++ )
        {
            Slb += SpreadConsX.getDomVarLb(i);
            Sub += SpreadConsX.getDomVarUb(i);
        }
        double minUcalc = Slb/SpreadConsX.getDomN();
        double maxUcalc = Sub/SpreadConsX.getDomN();
        if (minUcalc > SpreadConsU.getDomVarLb(0))
        {
            if (minUcalc > SpreadConsU.getDomVarUb(0))
            {
                cout << "No feasible solution as minUCalc from X is > Umax"  << endl;
                return 0;
            }
            SpreadConsU.setDomVarLb(0, ceil(minUcalc));
            changed = true;
        }
        if (maxUcalc < SpreadConsU.getDomVarUb(0))
        {
            if (maxUcalc < SpreadConsU.getDomVarLb(0))
            {
                cout << "No feasible solution as maxUCalc from X is < Umin"  << endl;
                return 0;
            }
            SpreadConsU.setDomVarUb(0, floor(maxUcalc));
            changed = true;
        }
        cout << "Propagate U using X" <<endl;
        SpreadConsU.printDomVar();
        //------------------------------------------------------------
        // Prune the Mean Bounds using S

        // Calculate the max N*Variance based on given std. deviation S
        // Note: For SCIP , due to Bound Consistency, can only prune using max variance and not minimum.
        double nVarUbCalc = SpreadConsX.getDomN() * SpreadConsS.getDomVarUb(0) * SpreadConsS.getDomVarUb(0);
//   cout << "nVarUbCalc " << nVarUbCalc << endl;

        // Check for consistency of solution
        double lowestPossibleNVar = SpreadConsX.I[0]->nVarLb; // Initialize to first lower bound
        for (int i = 0; i < SpreadConsX.numI; i++)
        {
            // If this interval contains the global minimum
            if ((SpreadConsX.I[i]->q0 <= SpreadConsX.I[i]->Vub) && (SpreadConsX.I[i]->q0 >= SpreadConsX.I[i]->Vlb))
            {
                lowestPossibleNVar = SpreadConsX.I[i]->C + (((SpreadConsX.I[i]->q0 - SpreadConsX.I[i]->ES)*(SpreadConsX.I[i]->q0 - SpreadConsX.I[i]->ES))/SpreadConsX.I[i]->m) - ((SpreadConsX.I[i]->q0 * SpreadConsX.I[i]->q0)/SpreadConsX.getDomN());
                break; // break out of for loop cause guaranteed to be minimal
            }
            // If it does not contain the global minimum, one of the end points must be reaching the global minimum,
            // thus, set it to one of the end points.
            else
            {
                if(SpreadConsX.I[i]->nVarLb < lowestPossibleNVar)
                {
                    lowestPossibleNVar = SpreadConsX.I[i]->nVarLb;
                }
                if(SpreadConsX.I[i]->nVarUb < lowestPossibleNVar)
                {
                    lowestPossibleNVar = SpreadConsX.I[i]->nVarUb;
                }
            }
        }
        cout << "nVarUBCalc is " << nVarUbCalc << " whereas lowestPossibleNvar is " << lowestPossibleNVar << endl;
        if (lowestPossibleNVar > nVarUbCalc)
        {
        // Note: Can't do this as it turns out the algorithm below doesn't work
        // on the specific example of [1,3], [2,6]  , [3,9] variables
        // It will end up pruning to [3,3] for all 3 variables
        // which prunes S to [0,0], which makes algorithm doesn't work.
        // However, still need above code to solve for nVarUbCalc
        // which is used below
        //    cout << "NO FEASIBLE SOLUTION " << endl
        //       << "as maximum n*Variance from given S is lower " << endl <<"than minimum possible n*Variance from given X values " << endl;
        //            return -1; // NO FEASIBLE SOLUTION
        }

        // Solve for q to get nVarUbCalc for each loop
        double qtemp1 = 0;
        double qtemp2 = 0;
        double qMeanUb = SpreadConsU.getDomVarUb(0) * SpreadConsX.getDomN();
        double qMeanLb = SpreadConsU.getDomVarLb(0) * SpreadConsX.getDomN();
        // TODO:
        // NOTE: I am not sure if the m = 0 case causes any problems here.
        // Cause if m = 0, the graph can no longer be proven to be convex as the 2nd derivative
        // will be < 0 and is concave instead.

        double a = 0;
        double b = 0; // note: This is b' which is b/2
        double c = 0;
        double inRoot = 0;
        for (int i = 0; i < SpreadConsX.numI; i++)
        {
            // Note: Always write 1.0 instead of 1 when working with doubles
            // if not it will divide as integers and result in an integer answer
            a =((1.0/SpreadConsX.I[i]->m) - (1.0/SpreadConsX.getDomN()));
            b = ((-1.0) * (SpreadConsX.I[i]->ES / SpreadConsX.I[i]->m));
            c = (((SpreadConsX.I[i]->ES * SpreadConsX.I[i]->ES)/((double)SpreadConsX.I[i]->m))+ SpreadConsX.I[i]->C - nVarUbCalc);
            inRoot = b*b - a*c;
            qtemp1 =  ((((-1.0) * b) + sqrt(inRoot))/a);
            qtemp2 = ((((-1.0) * b) - sqrt(inRoot))/a);
// cout << "a: "<< a << " b: " << b << " c: " << c << " inRoot: " << inRoot <<endl;
// cout << "qtemp: " << qtemp1 << "  " << qtemp2 << endl;

            // if the calculated qtemp1 is within the bounds of this interval
            if (qtemp1 >= SpreadConsX.I[i]->Vlb && qtemp1 <= SpreadConsX.I[i]->Vub )
            {
                if (ceil(qtemp1/SpreadConsX.getDomN()) < SpreadConsU.getDomVarUb(0))
                {
                    qMeanUb = qtemp1; // update qMeanUb
                    SpreadConsU.setDomVarUb(0, ceil(qtemp1/SpreadConsX.getDomN()));
                    changed = true;
                    cout << "Success: Pruned U's upper bound to " <<ceil(qtemp1/SpreadConsX.getDomN()) << endl;
                }
            }
            if (qtemp2 >= SpreadConsX.I[i]->Vlb && qtemp2 <= SpreadConsX.I[i]->Vub )
            {
                if (floor(qtemp2/SpreadConsX.getDomN()) > SpreadConsU.getDomVarLb(0))
                {
                    qMeanLb = qtemp2; // update qMeanLb
                    SpreadConsU.setDomVarLb(0, floor(qtemp2/SpreadConsX.getDomN()));
                    changed = true;
                    cout << "Success: Pruned U's lower bound to " <<floor(qtemp2/SpreadConsX.getDomN()) << endl;
                }
            }
        }
        //------------------------------------------------------------

        // Note: This helps in dealing with weird arithmetic error and calculations
        // However, it may result in a less tightly pruned algorithm.
        // It shouldn't matter much in terms of pruning though, cause this case is rare and is often
        // close to ideality.
        if (SpreadConsS.getDomVarUb(0) == 0)
        {
            break; // get out from propagation algorithm for X  if upper bound is 0
        }

        // Prune the X Bounds using S (U is used to determine qSelected)
        // Note: U helps by selecting qSelected
        double qSelected  = 0;
        int indexQSelected = 0;
        int originalIndexQSelected = 0; // to restore original indexQSelected at end of iteration
        int xIndex = 0;
        double DMax = 0;
        int originalLowerBound = 0;
        int maxShift = 0;
        double v = 0;
        double DMin = 0;
        int originalUpperBound = 0;

        // Need these 2 to handle deadling with changing intervals
        // and qSelected
        int XLB[SpreadConsX.getDomN()]; // get all the updated lower bounds
        int XUB[SpreadConsX.getDomN()]; // get all the updated upper bounds

        // For Given Mean Value, when lower bound is equal to upper bound
 //       if (SpreadConsU.getDomVarLb(0) == SpreadConsU.getDomVarUb(0))
// TODO: Uncomment above and comment below.
 if ( 2 == 1) // for debugging else case
//if (1 == 1)
        {
            // note: It is called qSelected here NOT q0, as q0 refers to something else to propagate U below
            qSelected = (double) SpreadConsU.getDomVarLb(0) * SpreadConsX.getDomN();
            //   qSelected = 10; // For Debugging
            cout << "Pruning X using GIVEN U with qSelected: " << qSelected << endl;
            // Determine which interval qSelected is in
            indexQSelected = 0;
            for (int i = 0; i < SpreadConsX.numI; i++)
            {
                if ((qSelected >= SpreadConsX.I[i]->Vlb) && (qSelected <= SpreadConsX.I[i]->Vub))
                {
                    indexQSelected = i;
                    break;
                }
            }
            originalIndexQSelected = indexQSelected;

            xIndex = 0;
            DMax = 0;
            originalLowerBound = 0;
            maxShift = 0;
            // Prune the X variables on R(I) upper bound
            for (int i = 0; i < SpreadConsX.I[originalIndexQSelected]->r; i++)
            {
                xIndex = SpreadConsX.I[originalIndexQSelected]->R[i];
                // Save the originalLowerBound of Current X as it will be changed
                originalLowerBound = SpreadConsX.getDomVarLb(xIndex);
                DMax =  SpreadConsX.findDMax(qSelected,originalIndexQSelected, xIndex, SpreadConsS.getDomVarUb(0));
                cout << "Dmax R(Iq) for X[" << xIndex <<"] is: " << DMax << endl;
                SpreadConsX.setDomVarLb(xIndex, originalLowerBound); // reset original lower bound
                maxShift = ceil(SpreadConsX.getDomVarLb(xIndex) + DMax);
                if (maxShift < SpreadConsX.getDomVarUb(xIndex))
                {
                        // SpreadConsX.setDomVarUb(xIndex, maxShift);
                    XUB[xIndex] = maxShift;
                    changed = true;
                    cout << "Success: Pruned X[" << xIndex << "]'s upper bound to " <<maxShift << endl;
                }
                else
                {
                    XUB[xIndex] = SpreadConsX.getDomVarUb(xIndex);
                }
                // If change upper bound, will not change lower bound,
                XLB[xIndex] = SpreadConsX.getDomVarLb(xIndex);
                // Recreate the intervals to original
                SpreadConsX.createInterval();
            }
            v = 0;
            // Prune the X variables on M(I) upper bound
            // Note: If you did not have orignalIndexQSelected and just indexQSelected,
            // the bottom code will change the value of indexQSelected

            for (int i = 0; i < SpreadConsX.I[originalIndexQSelected]->m ; i++)
            {
                xIndex = SpreadConsX.I[originalIndexQSelected]->M[i];
                cout<< " m is " << xIndex << endl;
                // Save the originalLowerBound of Current X as it will change
                originalLowerBound = SpreadConsX.getDomVarLb(xIndex);
                // Since you can shift all the min to at least up to V, shift it
                // First , get the indexQSelected without shifting
                for (int j = 0; j < SpreadConsX.numI; j++)
                {
                    if ((qSelected >= SpreadConsX.I[j]->Vlb) && (qSelected <= SpreadConsX.I[j]->Vub))
                    {
                        indexQSelected = j; // get new indexQSelected
                        break;
                    }
                }
                // Calculate for v
                v = (qSelected - SpreadConsX.I[indexQSelected]->ES)/((double) SpreadConsX.I[indexQSelected]->m);
                DMax = v - SpreadConsX.getDomVarLb(xIndex); // can shift by at least v - current lower bound
                // Shift current xIndex to v
                SpreadConsX.setDomVarLb(xIndex, floor(v));
                SpreadConsX.createInterval(); // create the new intervals
                // Find the new intervals for indexQSelected
                for (int j = 0; j < SpreadConsX.numI; j++)
                {
                    if ((qSelected >= SpreadConsX.I[j]->Vlb) && (qSelected <= SpreadConsX.I[j]->Vub))
                    {
                        indexQSelected = j;
                        break;
                    }
                }
                DMax +=  SpreadConsX.findDMax(qSelected,indexQSelected, xIndex, SpreadConsS.getDomVarUb(0));
                cout << "Dmax M(Iq) for X[" << xIndex <<"] is: " << DMax << endl;
                SpreadConsX.setDomVarLb(xIndex, originalLowerBound); // reset original lower bound
                maxShift = ceil(SpreadConsX.getDomVarLb(xIndex) + DMax);

                if (maxShift < SpreadConsX.getDomVarUb(xIndex))
                {
                    //SpreadConsX.setDomVarUb(xIndex, maxShift);
                    XUB[xIndex] = maxShift;
                    changed = true;
                    cout << "Success: Pruned X[" << xIndex << "]'s upper bound to " <<maxShift << endl;
                }
                else
                {
                    XUB[xIndex] = SpreadConsX.getDomVarUb(xIndex);
                }
                // If X is in M(Iq), will need to prune both ways
                // Reset the  intervals to original
                SpreadConsX.createInterval();
            }

            DMin = 0;
            originalUpperBound = 0;
            // Prune the X variables on L(I) lower bound
            for (int i = 0; i < SpreadConsX.I[originalIndexQSelected]->l; i++)
            {
                xIndex = SpreadConsX.I[originalIndexQSelected]->L[i];
                originalUpperBound = SpreadConsX.getDomVarUb(xIndex);
                DMin =  SpreadConsX.findDMin(qSelected,originalIndexQSelected, xIndex, SpreadConsS.getDomVarUb(0));
                cout << "Dmin L(Iq) for X[" << xIndex <<"] is: " << DMax << endl;
                SpreadConsX.setDomVarUb(xIndex, originalUpperBound); // reset original lower bound
                maxShift = floor(SpreadConsX.getDomVarUb(xIndex) - DMin);
                if (maxShift > SpreadConsX.getDomVarLb(xIndex))
                {
                    XLB[xIndex] = maxShift;
                    //SpreadConsX.setDomVarLb(xIndex, maxShift);
                    changed = true;
                    cout << "Success: Pruned X[" << xIndex << "]'s lower bound to " <<maxShift << endl;
                }
                else
                {
                    XLB[xIndex] =  SpreadConsX.getDomVarLb(xIndex);
                }
                XUB[xIndex] = SpreadConsX.getDomVarUb(xIndex);
                // Reset to original intervals
                SpreadConsX.createInterval();
            }// end of L algorithm

            v = 0;
            DMin = 0;
            originalUpperBound = 0;
            // Prune the X variables on M(I) lower bound
            for (int i = 0; i < SpreadConsX.I[originalIndexQSelected]->m ; i++)
            {
                xIndex = SpreadConsX.I[originalIndexQSelected]->M[i];
                cout<< " m is " << xIndex << endl;
                // Save the originalLowerBound of Current X as it will change
                originalUpperBound = SpreadConsX.getDomVarUb(xIndex);
                // Since you can shift all the min to at least up to V, shift it
                // First , get the indexQSelected without shifting
                // Calculate for v
                v = (qSelected - SpreadConsX.I[originalIndexQSelected]->ES)/((double) SpreadConsX.I[originalIndexQSelected]->m);
                DMin = SpreadConsX.getDomVarUb(xIndex) - v; // can shift by at least current upper bound - v to v
                // Shift current xIndex to v
                SpreadConsX.setDomVarUb(xIndex, ceil(v));
                SpreadConsX.createInterval(); // create the new intervals
                // Find the new intervals for indexQSelected
                for (int j = 0; j < SpreadConsX.numI; j++)
                {
                    if ((qSelected >= SpreadConsX.I[j]->Vlb) && (qSelected <= SpreadConsX.I[j]->Vub))
                    {
                        indexQSelected = j;
                        break;
                    }
                }
                DMin +=  SpreadConsX.findDMin(qSelected,indexQSelected, xIndex, SpreadConsS.getDomVarUb(0));
                cout << "Dmin M(Iq) for X[" << xIndex <<"] is: " << DMin << endl;
                SpreadConsX.setDomVarUb(xIndex, originalUpperBound); // reset original upper bound
                maxShift = floor(SpreadConsX.getDomVarUb(xIndex) - DMin);

                if (maxShift  > SpreadConsX.getDomVarLb(xIndex))
                {
                    //SpreadConsX.setDomVarUb(xIndex, maxShift);
                    XLB[xIndex] = maxShift;
                    changed = true;
                    cout << "Success: Pruned X[" << xIndex << "]'s lower bound to " <<maxShift << endl;
                }
                else
                {
                    XLB[xIndex] = SpreadConsX.getDomVarLb(xIndex);
                }
                // If X is in M(Iq), will need to prune both ways
                // Reset the  intervals to original
                SpreadConsX.createInterval();
            }

            // Note: Important TODO:
            // Doing this below that you only update values after checking all X
            // is incorrect as you may over-prune but it makes the code a lot faster than
            // simply re-calculating everthing as soon as you find one single X that gets to be pruned.

            // Update all bound values
            for(i = 0; i < SpreadConsX.getDomN(); i++)
            {
                SpreadConsX.setDomVarLb(i, XLB[i]);
                SpreadConsX.setDomVarUb(i, XUB[i]);
            }
            // Recreate the new intervals
            SpreadConsX.createInterval();
            // Re-print latest bounds
            cout<<"Latest Bounds:" << endl;
            SpreadConsX.printDomVar();
            SpreadConsU.printDomVar();
            SpreadConsS.printDomVar();
        } // End of given mean algorithm

        // For general mean values, when lower bound not equal to upper bound of mean
        // Note: The paper is wrong as the DMax(q) is shown to be convex and not derivable from Matlab Code
        // I have confirmed this error with the author himself via email.
        else
        {
            cout << "Pruning X using GENERAL U" << endl;
            qSelected  = SpreadConsU.getDomVarLb(0); // Initialize to the lower bound
            indexQSelected = 0;
            originalIndexQSelected = 0; // to restore original indexQSelected at end of iteration
            xIndex = 0;
            DMax = 0;
            originalLowerBound = 0;
            maxShift = 0;
            v = 0;
            DMin = 0;
            originalUpperBound = 0;
            double XMin = 0; // current XMin used to calculate DMax for L(Iq), R(Iq)
            double XMax = 0; // current XMax used to calculate DMin for L(Iq), R(Iq)
            double QMax = 0; // Qmax from current U's upper bound
            double QMin = 0; // Qmin from current U's lower bound
            double stdDevUb = SpreadConsS.getDomVarUb(0);
            double optMax = 1.0*stdDevUb*stdDevUb*SpreadConsX.getDomN();
            QMax = SpreadConsU.getDomVarUb(0) * SpreadConsX.getDomN();
            QMin = SpreadConsU.getDomVarLb(0) * SpreadConsX.getDomN();
            stdDevUb = SpreadConsS.getDomVarUb(0);
            double dMaxTemp = 0; // to get value of dMax from DMAX() function
            double dMaxFinal = 0; // initialize to minimum possible value
            double dMinTemp = 0;  // to get value of dMin from DMIN() function
            double dMinFinal = 0; // initialize to minimum possible value
            double qArgumentLB = 0;
            double qArgumentUB = 0;
            bool pruneUseR = false;
            bool pruneUseL = false;
            bool pruneUseM = false;
            bool stopPruning = false;
            // Loop through each variable and calculate qMax to be used fpr R, L, M or none
            for(int s = 0; s < SpreadConsX.getDomN(); s++)
            {
                // Note: Each X can either belong to (R or M) & (L or M) but not both
                // Stop pruning is used to not bother looking for X in M if already prune for R for DMax case
                pruneUseR = false;
                pruneUseM = false;
                stopPruning = false;
                pruneUseL = false;

                XMin = SpreadConsX.getDomVarLb(s);
                XMax = SpreadConsX.getDomVarUb(s);
                xIndex = s;
                bool localMax = false; // false means not yet detected local max yet


                //--------------------------------------------------------------------------------------
                // Prune Upper Bound of X using R(Iq) and M(Iq)
                //--------------------------------------------------------------------------------------
                // Loop through every interval and keep track of which interval has the largest dMaxTemp
                for(int t = 0; t < SpreadConsX.numI; t++)
                {
                    qArgumentLB = SpreadConsX.I[t]->Vlb;
                    qArgumentUB = SpreadConsX.I[t]->Vub;
                    // Case 1: Refer to left end point (lower bound)
                    // Note: If you sketch the graph, you can have 6 different scenarios,
                    // each scenario below covers 2 of the 6 by these inequality tests
                    // Scenario 1:  (qArgumentLB >= QMin) && (qArgumentLB <= QMax)
                                    // set qArgumentLB to qArgumentLB
                    // Scenario 2: ( qArgumentLB<= QMin) && (qArgumentUB >= QMin)
                                    // set qArgumentLB to QMin
                    // Scenario 3:  (qArgumentUB < QMin) || (qArgumentLB > QMax)
                                    // continue; cause nothing can be done (including for Case 2)
                    // Scenario 1
                    if((qArgumentLB >= QMin) && (qArgumentLB <= QMax))
                    {
                        // Do nothing cause already assigned correctly
                    }
                    // Scenario 2
                    else if((qArgumentLB<= QMin) && (qArgumentUB >= QMin))
                    {
                        qArgumentLB = QMin;
                    }
                    // Scenario 3
                    else if ((qArgumentUB < QMin) || (qArgumentLB > QMax))
                    {
                        continue; // go to next iteration of loop
                    }
                    else
                    {
                        cout << " SHOULDN'T BE IN HERE" << endl;
                        continue;
                    }
                    dMaxTemp = SpreadConsX.DMAX(qArgumentLB,XMin,optMax,0);
                    cout << "dMaxTemp: " << dMaxTemp << " dMaxFinal: " << dMaxFinal << endl;
                    if (dMaxTemp > dMaxFinal)
                    {
                        dMaxFinal = dMaxTemp;
                        qSelected = qArgumentLB;
                    }
                    // Only need execute this if localMax not found yet,
                    if (localMax == false)
                    {
                        // Check right end point
                        // Case 2: Refer to right end point (upper bound)
                        // Scenario 1:  (qArgumentUB >= QMin) && (qArgumentUB <= QMax)
                                        // set qArgumentUB to qArgumentUB
                        // Scenario 2: (qArgumentUB >= QMax) && (qArgumentLB <= QMax)
                                        // set qArgumentLB to QMin
                        // Scenario 3:  (qArgumentUB < QMin) || (qArgumentLB > QMax)
                        // Scenario 1
                        if ((qArgumentUB >= QMin) && (qArgumentUB <= QMax))
                        {
                            // do nothing, qArgumentUB is already the right value
                        }
                        // Scenario 2
                        else if((qArgumentUB >= QMax) && (qArgumentLB <= QMax))
                        {
                            qArgumentUB = QMax;
                        }
                        // Scenario 3
                        else if ((qArgumentUB < QMin) || (qArgumentLB > QMax))
                        {
                            cout << " Should never be in here since above Case 1 should have continue to next iteration" << endl;
                            continue;
                        }
                        else
                        {
                            cout << " Definitely should never be in here" << endl;
                            continue;
                        }
                        dMaxTemp = SpreadConsX.DMAX(qArgumentUB,XMin,optMax,1);
                        if (dMaxTemp > dMaxFinal)
                        {
                            dMaxFinal = dMaxTemp;
                            qSelected = qArgumentUB;
                        }
                        // Check if contain local max
                        // Variables to calculate qSelected or qLast
                        double qLast = 0; // The optimal value of q from the last interval
                        double qLastVer1 = 0, qLastVer2 = 0, qLastVer3 = 0; // for debugging
                        double qLastSelected = 0; // from qLast to determine if bigger or smaller
                        double a = 0;
                        double b = 0;
                        double c = 0;
                        double d = 0;

                        double m = 0;  double n = 0;
                        double inRoot = 0;

                        // Calcualte for the optimal qLast
                        m = SpreadConsX.I[t]->m;
                        n = SpreadConsX.getDomN();
  cout << " CHECKING WORK" << endl;
  cout << "M: " << m << " N: " << n << endl;

                        a = ((-1.0)/m + ((1.0 + 1.0/m)*(1.0/n)));
                        b = ((-2.0)/m) * (XMin + SpreadConsX.I[t]->ES);
                        c = pow((double)XMin, 2.0) + ((2.0*(XMin)*SpreadConsX.I[t]->ES)/(double)m)
                        - (1.0+1.0/m)*(SpreadConsX.I[t]->C - (pow((double) stdDevUb, 2.0)*SpreadConsX.getDomN()))
                         - (pow( (double)SpreadConsX.I[t]->ES, 2.0)/m);
                        d = 1.0/m;
  cout << "a: " << a << " b: " << b << " c: " << c << " d: " << d <<  endl;

                        inRoot = ((-1.0)*b*b + 4.0*a*c)*((-1.0)*d*d + a);
                        qLastVer1 = ((d* sqrt(inRoot) + a*b - b*d*d)/(2.0*((-1.0)*a*a + a*d*d)));
                        if (qLastVer1 < 0)
                        {
                            qLastVer1 = (-1.0)*qLastVer1; // Version 1: My theory with Matlab Solving the algebraic eqn
                        }                           // Note: If  use version 1, need account for infinity answers
            cout << "Version 1 qLast: " << qLastVer1 << endl;
            // Note: Although will get infinity for q, answer still remains correct.


// Using Purely Matlab Solution
// % Version 2: Solution 1
// % (ES*n^2 - Xjmin*n^2 + n*(-(m - n)*(C + C*m - C*n - OptMax^2*m + OptMax^2*n - Xjmin^2*m + Xjmin^2*n + ES^2 - OptMax^2 - 2*ES*Xjmin))^(1/2) - ES*m*n + Xjmin*m*n)/(m^2 - 2*m*n + m + n^2 - n)

// % Version 3: Solution 2
// % -(Xjmin*n^2 - ES*n^2 + n*(-(m - n)*(C + C*m - C*n - OptMax^2*m + OptMax^2*n - Xjmin^2*m + Xjmin^2*n + ES^2 - OptMax^2 - 2*ES*Xjmin))^(1/2) + ES*m*n - Xjmin*m*n)/(m^2 - 2*m*n + m + n^2 - n)
                        double ES, C, OptMax, Xjmin;
                        Xjmin = (double) 1.0*XMin;
                        OptMax = (double) 1.0*stdDevUb*stdDevUb*n;
                        C =  (double) SpreadConsX.I[t]->C;
                        ES = (double) SpreadConsX.I[t]->ES;

  cout << "Xjmin: " << Xjmin << " OptMax: " << OptMax << " C: " << C << " ES: " << ES <<  endl;

double inRootTEMP1 = -(m - n)*(C + C*m - C*n - OptMax*OptMax*m + OptMax*OptMax*n - Xjmin*Xjmin*m + Xjmin*Xjmin*n + ES*ES - OptMax*OptMax - 2.0*ES*Xjmin);
double inRootTEMP2 = -(m - n)*(C + C*m - C*n - OptMax*OptMax*m + OptMax*OptMax*n - Xjmin*Xjmin*m + Xjmin*Xjmin*n + ES*ES - OptMax*OptMax - 2.0*ES*Xjmin);
cout << "inRootTemp1: " <<inRootTEMP1 << "inRootTemp2: " << inRootTEMP2 << endl;
                        // Solution 1 from purely Matlab
                        qLastVer2 = ((ES*n*n - Xjmin*n*n + n*sqrt(-(m - n)*(C + C*m - C*n - OptMax*OptMax*m + OptMax*OptMax*n - Xjmin*Xjmin*m + Xjmin*Xjmin*n + ES*ES - OptMax*OptMax - 2.0*ES*Xjmin)) - ES*m*n + Xjmin*m*n)/((double)(1.0*(m*m - 2*m*n + m + n*n - n))));
            cout << "Version 2 qLast: " << qLastVer2 << endl;

                        // Solution 2
                        qLastVer3 = (-1.0*(Xjmin*n*n - ES*n*n + n*sqrt(-(m - n)*(C + C*m - C*n - OptMax*OptMax*m + OptMax*OptMax*n - Xjmin*Xjmin*m + Xjmin*Xjmin*n + ES*ES - OptMax*OptMax - 2.0*ES*Xjmin)) + ES*m*n - Xjmin*m*n)/((double)((1.0)*(m*m - 2*m*n + m + n*n - n))));
            cout << "Version 3 qLast: " << qLastVer3 << endl;

                        if(qLastVer2 > 0)
                        {
                            qLast = qLastVer2;
                        }
                        else if (qLastVer3 > 0)
                        {
                            qLast = qLastVer3;
                        }
                        else
                        {
                            qLast = qLastVer1;
                        }
            cout << "Chosen qLast: " << qLast << endl;

// FIXME: qLastVer1 has inf sometimes when a = 0, qLastVer2 & qLastVer3 has inf and nan sometimes

                        if ((qLast >= SpreadConsX.I[t]->Vlb) && (qLast <= SpreadConsX.I[t]->Vub))
                        {
                            // qLast is in this interval, assign it to proper X Values
                            // compare it with  the sides of all other intervals to ensure
                            // qSelected is indeed assigned to the maximum
                            localMax  = true; // localMax was found, don't have to do this again.

                            // TODO: Check if the DMax with current qLast is bigger than current
                            // and update qSelected if it is
                            // TODO: Check the equation and see how it will be divided by 0,
                            // then, account for all special cases.

                            // Case 1: Mean  Bounds (Q0 and Q1) covers qlast
                            if ((qLast >= QMax)&& (qLast <= QMin))
                            {
                                //qLast = qLast;
                            }
                            else if (qLast <= QMin)
                            {
                                qLast = QMin;
                            }
                            else // if (qLast >= QMax);
                            {
                                qLast = QMax;
                            }

                            dMaxTemp = SpreadConsX.DMAX(qLast,XMin,optMax,0);
                            cout << "dMaxTemp: " << dMaxTemp << " dMaxFinal: " << dMaxFinal << endl;
                            if (dMaxTemp > dMaxFinal)
                            {
                                dMaxFinal = dMaxTemp;
                                qSelected = qLast;
                            }
                        }
                        else
                        {
                            // if qLast does not belong to this interval,
                            continue;
                        }
                    } // end of if localMax found
                    // else, skip this interval and move on to next interval
                    // NOTE: A backup fix would be to loop through all q from U and pick the highest DMax.
                } // End of current interval
                // Here, have finish assigning dMaxFinal and qSelected for R(Iq) and M(Iq)
                // To prune upper bound of X
                // Need calculate dMax from qSelected to prune R(Iq) and M(Iq)
                // First, get QIndex
                cout << "Using DMAX to prune upper bound of current X" << endl;
                cout << "xIndex is " << xIndex << endl;
                cout << "qSelected is " << qSelected << endl;
                // Since you can shift all the min to at least up to V, shift it
                // First , get the indexQSelected without shifting
                for (int j = 0; j < SpreadConsX.numI; j++)
                {
                    if ((qSelected >= SpreadConsX.I[j]->Vlb) && (qSelected <= SpreadConsX.I[j]->Vub))
                    {
                        indexQSelected = j; // get new indexQSelected
                        break;
                    }
                }
                cout << "indexQSelected is " << indexQSelected << endl;
                originalIndexQSelected = indexQSelected;
                for (int i = 0; i < SpreadConsX.I[originalIndexQSelected]->r; i++)
                {
                    // If this X is part of R
                    if (xIndex == SpreadConsX.I[originalIndexQSelected]->R[i])
                    {
                        pruneUseR = true;
                        stopPruning = true; // don't bother pruning M for this X
                    }
                }
                // Step 1: Prune this X using R(Iq) if X is in R(Iq)
                if (pruneUseR)
                {
                    cout << "Pruning upperbound using R(Iq)" << endl;
                    // Save the originalLowerBound of Current X as it will change
                    originalIndexQSelected = indexQSelected;
                    DMax = 0;
                    originalLowerBound = 0;
                    maxShift = 0;
                    // Prune current X on R(I) upper bound
                    // Save the originalLowerBound of Current X as it will be changed
                    originalLowerBound = SpreadConsX.getDomVarLb(xIndex);
                    DMax =  SpreadConsX.findDMax(qSelected,originalIndexQSelected, xIndex, SpreadConsS.getDomVarUb(0));
                    cout << "Dmax R(Iq) for X[" << xIndex <<"] is: " << DMax << endl;
                    SpreadConsX.setDomVarLb(xIndex, originalLowerBound); // reset original lower bound
                    maxShift = ceil(SpreadConsX.getDomVarLb(xIndex) + DMax);
                    if (maxShift < SpreadConsX.getDomVarUb(xIndex))
                    {
                        // SpreadConsX.setDomVarUb(xIndex, maxShift);
                        XUB[xIndex] = maxShift;
                        changed = true;
                        cout << "Success: Pruned X[" << xIndex << "]'s upper bound to " <<maxShift << endl;
                    }
                    else
                    {
                        XUB[xIndex] = SpreadConsX.getDomVarUb(xIndex);
                    }
                    // Recreate the intervals to original
                    SpreadConsX.createInterval();
                } // end of if(pruneUseR)

                // If X wasn't in R
                if(stopPruning == false)
                {
                    for (int i = 0; i < SpreadConsX.I[originalIndexQSelected]->m; i++)
                    {
                        // If this X is part of M
                        if (xIndex == SpreadConsX.I[originalIndexQSelected]->M[i])
                        {
                            pruneUseM = true;
                            stopPruning = true; // don't setting default values for XUB for this X
                        }
                    }
                }
                // Step 2: Prune this X using M(Iq)
                if (pruneUseM)
                {
                    cout << "Pruning upperbound using M(Iq)" << endl;

                    v = 0;
                    // Save the originalLowerBound of Current X as it will change
                    originalLowerBound = SpreadConsX.getDomVarLb(xIndex);
                    // Since you can shift all the min to at least up to V, shift it
                    // First , get the indexQSelected without shifting
                    // If in here, m won't be 0, so can definitely calculate for v
                    // Calculate for v
                    v = (qSelected - SpreadConsX.I[originalIndexQSelected]->ES)/((double) SpreadConsX.I[originalIndexQSelected]->m);
                    DMax = v - SpreadConsX.getDomVarLb(xIndex); // can shift by at least v - current lower bound
                    // Shift current xIndex to v
                    SpreadConsX.setDomVarLb(xIndex, floor(v));
                    SpreadConsX.createInterval(); // create the new intervals
                    // Find the new intervals for indexQSelected
                    for (int i = 0; i< SpreadConsX.numI; i++)
                    {
                        if ((qSelected >= SpreadConsX.I[i]->Vlb) && (qSelected <= SpreadConsX.I[i]->Vub))
                        {
                            indexQSelected = i;
                            break;
                        }
                    }
                    DMax +=  SpreadConsX.findDMax(qSelected,indexQSelected, xIndex, SpreadConsS.getDomVarUb(0));
                    cout << "Dmax M(Iq) for X[" << xIndex <<"] is: " << DMax << endl;
                    SpreadConsX.setDomVarLb(xIndex, originalLowerBound); // reset original lower bound
                    maxShift = ceil(SpreadConsX.getDomVarLb(xIndex) + DMax);

                    if (maxShift < SpreadConsX.getDomVarUb(xIndex))
                    {
                        //SpreadConsX.setDomVarUb(xIndex, maxShift);
                        XUB[xIndex] = maxShift;
                        changed = true;
                        cout << "Success: Pruned X[" << xIndex << "]'s upper bound to " <<maxShift << endl;
                    }
                    else
                    {
                        XUB[xIndex] = SpreadConsX.getDomVarUb(xIndex);
                    }
                    // Reset the  intervals to original
                    SpreadConsX.createInterval();
                }// end of if(pruneUseM)

                // If did not prune using R or M for this X, set to original value
                if(stopPruning == false)
                {
                    XUB[xIndex] = SpreadConsX.getDomVarUb(xIndex);
                }

                //--------------------------------------------------------------------------------------
                // Prune Lower Bound of X using L(Iq) and M(Iq)
                //--------------------------------------------------------------------------------------
                // Find qSelect and indexQSelect using DMIN
                // Loop through every interval and keep track of which interval has the largest dMaxTemp

                // Note: Each X can either belong to (R or M) & (L or M) but not both
                // Stop pruning is used to not bother looking for X in M if already prune for R for DMax case
                pruneUseL = false;
                pruneUseM = false;
                stopPruning = false;

                for(int t = 0; t < SpreadConsX.numI; t++)
                {
                    qArgumentLB = SpreadConsX.I[t]->Vlb;
                    qArgumentUB = SpreadConsX.I[t]->Vub;
                    // Case 1: Refer to left end point (lower bound)
                    // Note: If you sketch the graph, you can have 6 different scenarios,
                    // each scenario below covers 2 of the 6 by these inequality tests
                    // Scenario 1:  (qArgumentLB >= QMin) && (qArgumentLB <= QMax)
                                    // set qArgumentLB to qArgumentLB
                    // Scenario 2: ( qArgumentLB<= QMin) && (qArgumentUB >= QMin)
                                    // set qArgumentLB to QMin
                    // Scenario 3:  (qArgumentUB < QMin) || (qArgumentLB > QMax)
                                    // continue; cause nothing can be done (including for Case 2)
                    // Scenario 1
                    if((qArgumentLB >= QMin) && (qArgumentLB <= QMax))
                    {
                        // Do nothing cause already assigned correctly
                    }
                    // Scenario 2
                    else if((qArgumentLB<= QMin) && (qArgumentUB >= QMin))
                    {
                        qArgumentLB = QMin;
                    }
                    // Scenario 3
                    else if ((qArgumentUB < QMin) || (qArgumentLB > QMax))
                    {
                        continue; // go to next iteration of lSpreadConsXoop
                    }
                    else
                    {
                        cout << " SHOULDN'T BE IN HERE" << endl;
                        continue;
                    }
                    dMinTemp = SpreadConsX.DMIN(qArgumentLB,XMax,optMax,0);
                    if (dMinTemp > dMinFinal)
                    {
                        dMinFinal = dMinTemp;
                        qSelected = qArgumentLB;
                    }
                    // Case 2: Refer to right end point (upper bound)
                    // Scenario 1:  (qArgumentUB >= QMin) && (qArgumentUB <= QMax)
                                    // set qArgumentUB to qArgumentUB
                    // Scenario 2: (qArgumentUB >= QMax) && (qArgumentLB <= QMax)
                                    // set qArgumentLB to QMin
                    // Scenario 3:  (qArgumentUB < QMin) || (qArgumentLB > QMax)

                    // Scenario 1
                    if ((qArgumentUB >= QMin) && (qArgumentUB <= QMax))
                    {
                        // do nothing, qArgumentUB is already the right value
                    }
                    // Scenario 2
                    else if((qArgumentUB >= QMax) && (qArgumentLB <= QMax))
                    {
                        qArgumentUB = QMax;
                    }
                    // Scenario 3
                    else if ((qArgumentUB < QMin) || (qArgumentLB > QMax))
                    {
                        cout << " Should never be in here since above Case 1 should have continue to next iteration" << endl;
                        continue;
                    }
                    else
                    {
                        cout << " Definitely should never be in here" << endl;
                        continue;
                    }
                    dMinTemp = SpreadConsX.DMIN(qArgumentUB,XMax,optMax,1);
                    if (dMinTemp > dMinFinal)
                    {
                        dMinFinal = dMinTemp;
                        qSelected = qArgumentUB;
                    }
                } // End of current interval
                // Here, have finish assigning dMinFinal and qSelected for L(Iq) and M(Iq)
                // To prune lower bound of X
                // Need calculate dMin from qSelected to prune L(Iq) and M(Iq)
                // First, get QIndex
                cout << "Using DMin to prune lower bound of current X" << endl;
                cout << "xIndex is " << xIndex << endl;
                cout << "qSelected is " << qSelected << endl;
                // Since you can shift all the min to at least up to V, shift it
                // First , get the indexQSelected without shifting
                for (int j = 0; j < SpreadConsX.numI; j++)
                {
                    if ((qSelected >= SpreadConsX.I[j]->Vlb) && (qSelected <= SpreadConsX.I[j]->Vub))
                    {
                        indexQSelected = j; // get new indexQSelected
                        break;
                    }
                }
                cout << "indexQSelected is " << indexQSelected << endl;
                originalIndexQSelected = indexQSelected;
                for (int i = 0; i < SpreadConsX.I[originalIndexQSelected]->l; i++)
                {
                    // If this X is part of L
                    if (xIndex == SpreadConsX.I[originalIndexQSelected]->L[i])
                    {
                        pruneUseL = true;
                        stopPruning = true; // don't bother pruning M for this X
                    }
                }
                // Step 3: Prune this X using L(Iq) if X is in L(Iq)
                if (pruneUseL)
                {
                    cout << "Pruning lowerbound using L(Iq)" << endl;
                    // Save the originalLowerBound of Current X as it will change
                    DMin = 0;
                    originalUpperBound = 0;
                    maxShift = 0;
                    // Prune current X on L(I) upper bound
                    // Save the originalUpperBound of Current X as it will be changed
                    originalUpperBound = SpreadConsX.getDomVarUb(xIndex);
                    DMin =  SpreadConsX.findDMin(qSelected,originalIndexQSelected, xIndex, SpreadConsS.getDomVarUb(0));
                    cout << "Dmin L(Iq) for X[" << xIndex <<"] is: " << DMin << endl;
                    SpreadConsX.setDomVarUb(xIndex, originalUpperBound); // reset original upper bound
                    maxShift = floor(SpreadConsX.getDomVarUb(xIndex) - DMin);
                    if (maxShift > SpreadConsX.getDomVarLb(xIndex))
                    {
                        // SpreadConsX.setDomVarLb(xIndex, maxShift);
                        XLB[xIndex] = maxShift;
                        changed = true;
                        cout << "Success: Pruned X[" << xIndex << "]'s lower bound to " << maxShift << endl;
                    }
                    else
                    {
                        XLB[xIndex] = SpreadConsX.getDomVarLb(xIndex);
                    }
                    // Recreate the intervals to original
                    SpreadConsX.createInterval();
                } // end of if(pruneUseR)

                // If X wasn't in R
                if(stopPruning == false)
                {
                    for (int i = 0; i < SpreadConsX.I[originalIndexQSelected]->m; i++)
                    {
                        // If this X is part of M
                        if (xIndex == SpreadConsX.I[originalIndexQSelected]->M[i])
                        {
                            pruneUseM = true;
                            stopPruning = true; // don't setting default values for XUB for this X
                        }
                    }
                }
                // Step 4: Prune this X using M(Iq) if X is in M(Iq)
                if (pruneUseM)
                {
                    cout << "Pruning lowerbound using M(Iq)" << endl;
                    v = 0;
                    // Save the originalUpperBound of Current X as it will change
                    originalUpperBound = SpreadConsX.getDomVarUb(xIndex);
                    // Since you can shift all the max to at least down to V, shift it
                    // First , get the indexQSelected without shifting
                    // Note: If in here, definitely m isn't 0, so can calculate for v normally
                    // Calculate for v
                    v = (qSelected - SpreadConsX.I[originalIndexQSelected]->ES)/((double) SpreadConsX.I[originalIndexQSelected]->m);
                    DMin = SpreadConsX.getDomVarUb(xIndex) - v; // can shift by at least current upper bound - v to v
                    // Shift current xIndex to v
                    SpreadConsX.setDomVarUb(xIndex, ceil(v));
                    SpreadConsX.createInterval(); // create the new intervals
                    // Find the new intervals for indexQSelected
                    for (int i = 0; i < SpreadConsX.numI; i++)
                    {
                        if ((qSelected >= SpreadConsX.I[i]->Vlb) && (qSelected <= SpreadConsX.I[i]->Vub))
                        {
                            indexQSelected = i;
                            break;
                        }
                    }
                    DMin +=  SpreadConsX.findDMin(qSelected,indexQSelected, xIndex, SpreadConsS.getDomVarUb(0));
                    cout << "Dmin M(Iq) for X[" << xIndex <<"] is: " << DMin << endl;
                    SpreadConsX.setDomVarUb(xIndex, originalUpperBound); // reset original upper bound
                    maxShift = floor(SpreadConsX.getDomVarUb(xIndex) - DMin);
                    if (maxShift  > SpreadConsX.getDomVarLb(xIndex))
                    {
                        //SpreadConsX.setDomVarUb(xIndex, maxShift);
                        XLB[xIndex] = maxShift;
                        changed = true;
                        cout << "Success: Pruned X[" << xIndex << "]'s lower bound to " <<maxShift << endl;
                    }
                    else
                    {
                        XLB[xIndex] = SpreadConsX.getDomVarLb(xIndex);
                    }
                    // If X is in M(Iq), will need to prune both ways
                    // Reset the  intervals to original
                    SpreadConsX.createInterval();
                }// end of if(pruneUseM)

                // If did not prune using L or M for this X, set to original value
                if(stopPruning == false)
                {
                    XLB[xIndex] = SpreadConsX.getDomVarLb(xIndex);
                }
            } // End of current X

            // Note: Important : TODO:
            // Doing this below that you only update values after checking all X
            // is incorrect as you may over-prune but it makes the code a lot faster than
            // simply re-calculating everthing as soon as you find one single X that gets to be pruned.

            // After done checking all X
            // Update all bound values
            for(i = 0; i < SpreadConsX.getDomN(); i++)
            {
                SpreadConsX.setDomVarLb(i, XLB[i]);
                SpreadConsX.setDomVarUb(i, XUB[i]);
            }
            // Recreate the new intervals
            SpreadConsX.createInterval();
            // Re-print latest bounds
            cout<<"Latest Bounds:" << endl;
            SpreadConsX.printDomVar();
            SpreadConsU.printDomVar();
            SpreadConsS.printDomVar();
        } // End of Pruning X by General Case of Mean
        //------------------------------------------------------------
    }
    cout << " SOON SOON SOON SOON" << endl;
    // End of while loop, here, means can't propagate anymore cause nothing is changing
    SpreadConsX.printDomVar();
    SpreadConsS.printDomVar();
    SpreadConsU.printDomVar();
    // End of Propagation algorithm
    // In SCIP , just copy paste everything, and change the main() function to the CONSPROP() propagation function
    return 0;
}

//--------------------------------------------------------------------------------------------------------
// Previous code on based on errors that exist on paper that is no longer needed,
// left here for reference.
//--------------------------------------------------------------------------------------------------------
/* REMOVED THIS BELOW FROM PAPER AS ALREADY INCLUDED AS PART OF ALGORITHM ABOVE
// Note: This is already check in Prune U using X above
// Check feasible solution using X on U
if (Slb >  qSelected || Sub < qSelected)
{
    cout << "Slb: " << Slb << " Sub: " << Sub << " qSelected: " << qSelected << endl;
    cout << "No feasible solution cause mean is above or above extreme values" <<endl;
    return -1;
}
*/
//--------------------------------------------------------------------------------------------------------
/*
 OLD IMPLEMENTATION WHERE ASSUME DMAX(q) WAS CONCAVE AND DERIVABLE FROM PAPER'S CLAIMS.
        AS PAPER IS WRONG SINCE DMAX(q) IS CONVEX AND NOT CONTINUOUS AS VERIFIED FROM
        MATLAB AS WELL AS AUTHOR HIMSELF VIA EMAIL, I'VE WRITTEN A WHOLE NEW CODE
        THAT WORKS WITH THE NEW FINDINGS
        ALSO, THIS IMPLEMENTATION IS FILLED WITH BUGS THAT I DID NOT CORRECT
        I LEFT IT HERE JUST TO BE USED TO UNDERSTAND THE ALGORITHM I DEVELOPED FOR A FUTURE PROBLEM
        WHERE A CONCAVE GRAPH IS USED
        else // OLD VERSION THAT IS WRONG
        {
            cout << "Pruning X using GENERAL U" << endl;
            double XMin = 0; // current XMin used to calculate DMax for L(Iq), R(Iq)
            double XMax = 0; // current XMax used to calculate DMin for L(Iq), R(Iq)
            double QMax = 0; // Qmax from current U's upper bound
            double QMin = 0; // Qmin from current U's lower bound
            double stdDevUb = 0;

            QMax = SpreadConsU.getDomVarUb(0) * SpreadConsX.getDomN();
            QMin = SpreadConsU.getDomVarLb(0) * SpreadConsX.getDomN();
            stdDevUb = SpreadConsS.getDomVarUb(0);
            // Loop through each variable
            for(int s = 0; s < SpreadConsX.getDomN(); s++)
            {
                XMin = SpreadConsX.getDomVarLb(s);
                XMax = SpreadConsX.getDomVarUb(s);
                qSelected = -1; // Initialize to a value it can never be
                double qLast = 0; // The optimal value of q from the last interval
                double qLastVer1 = 0, qLastVer2 = 0, qLastVer3 = 0; // for debugging
                int lastCheckedIndex = 0; // The last checked index of the Intervals
                                         // will definitely be updated to 0 or something else
                                         // as one of the intervals must be checked since U's Q is a subspace of X's

                // Not needed as paper is wrong and graph is NOT concave
                // Variables to calculate qSelected or qLast
                double a = 0;
                double b = 0;
                double c = 0;
                double d = 0;

                double m = 0;  double n = 0;
                double inRoot = 0;

                // Note: e & g on paper gets removed to 0 on 2nd derivative

                // Loop through every interval
                for(int t = 0; t < SpreadConsX.numI; t++)
                {
                    // Make sure interval encompasses part of Q1 and Q2 from U
                    if ((SpreadConsX.I[t]->Ilb <= QMax)  && (SpreadConsX.I[t]->Iub >= QMin))
                    {
                        lastCheckedIndex = t;
                        // Calcualte for the optimal qLast
                        m = SpreadConsX.I[t]->m;
                        n = SpreadConsX.getDomN();
//  cout << " CHECKING WORK" << endl;
//  cout << "M: " << m << " N: " << n << endl;

                        a = ((-1.0)/m + ((1.0 + 1.0/m)*(1.0/n)));
                        b = ((-2.0)/m) * (XMin + SpreadConsX.I[t]->ES);
                        c = pow((double)XMin, 2.0) + ((2.0*(XMin)*SpreadConsX.I[t]->ES)/(double)m)
                        - (1.0+1.0/m)*(SpreadConsX.I[t]->C - (pow((double) stdDevUb, 2.0)*SpreadConsX.getDomN()))
                         - (pow( (double)SpreadConsX.I[t]->ES, 2.0)/m);
                        d = 1.0/m;
                        inRoot = ((-1.0)*b*b + 4.0*a*c)*((-1.0)*d*d + a);
                        qLastVer1 = ((d* sqrt(inRoot) + a*b - b*d*d)/(2.0*((-1.0)*a*a + a*d*d)));
                        if (qLastVer1 < 0)
                        {
                            qLastVer1 = (-1.0)*qLastVer1; // Version 1: My theory with Matlab Solving the algebraic eqn
                        }                           // Note: If  use version 1, need account for infinity answers
            cout << "Version 1 qLast: " << qLastVer1 << endl;

// Using Purely Matlab Solution
// % Version 2: Solution 1
// % (ES*n^2 - Xjmin*n^2 + n*(-(m - n)*(C + C*m - C*n - OptMax^2*m + OptMax^2*n - Xjmin^2*m + Xjmin^2*n + ES^2 - OptMax^2 - 2*ES*Xjmin))^(1/2) - ES*m*n + Xjmin*m*n)/(m^2 - 2*m*n + m + n^2 - n)

// % Version 3: Solution 2
// % -(Xjmin*n^2 - ES*n^2 + n*(-(m - n)*(C + C*m - C*n - OptMax^2*m + OptMax^2*n - Xjmin^2*m + Xjmin^2*n + ES^2 - OptMax^2 - 2*ES*Xjmin))^(1/2) + ES*m*n - Xjmin*m*n)/(m^2 - 2*m*n + m + n^2 - n)
                        double ES, C, OptMax, Xjmin;
                        Xjmin = (double) 1.0*XMin;
                        OptMax = (double) 1.0*stdDevUb*stdDevUb*n;
                        C =  (double) SpreadConsX.I[t]->C;
                        ES = (double) SpreadConsX.I[t]->ES;

                        // Solution 1
                        qLastVer2 = ((ES*n*n - Xjmin*n*n + n*sqrt(-(m - n)*(C + C*m - C*n - OptMax*OptMax*m + OptMax*OptMax*n - Xjmin*Xjmin*m + Xjmin*Xjmin*n + ES*ES - OptMax*OptMax - 2*ES*Xjmin)) - ES*m*n + Xjmin*m*n)/((double)(1.0*(m*m - 2*m*n + m + n*n - n))));
            cout << "Version 2 qLast: " << qLastVer2 << endl;

                        // Solution 2
                        qLastVer3 = (-1.0*(Xjmin*n*n - ES*n*n + n*sqrt(-(m - n)*(C + C*m - C*n - OptMax*OptMax*m + OptMax*OptMax*n - Xjmin*Xjmin*m + Xjmin*Xjmin*n + ES*ES - OptMax*OptMax - 2.0*ES*Xjmin)) + ES*m*n - Xjmin*m*n)/((double)((1.0)*(m*m - 2*m*n + m + n*n - n))));
            cout << "Version 3 qLast: " << qLastVer3 << endl;

                        if(qLastVer2 > 0)
                        {
                            qLast = qLastVer2;
                        }
                        else if (qLastVer3 > 0)
                        {
                            qLast = qLastVer3;
                        }
                        else
                        {
                            qLast = qLastVer1;
                        }
            cout << "Chosen qLast: " << qLast << endl;

                        if ((qLast >= SpreadConsX.I[t]->Vlb) && (qLast <= SpreadConsX.I[t]->Vub))
                        {
                               // qLast is in this interval, assign it to proper X Values

                            // Case 1: Mean  Bounds (Q0 and Q1) covers qlast
                            if ((qLast >= QMax)&& (qLast <= QMin))
                            {
                                qSelected = qLast;
                            }
                            else if (qLast <= QMin)
                            {
                                qSelected = QMin;
                            }
                            else // if (qLast >= QMax);
                            {
                                qSelected = QMax;
                            }
                        }
                        else
                        {
                            // if qLast does not belong to this interval,
                            continue;
                        }
                    }
                    // else, skip this interval and move on to next interval
                    continue;
                    // NOTE: A backup fix would be to loop through all q from U and pick the highest DMax.
                }
                // If q wasn't assigned
                if(qSelected < 0)
                {
                    // Look at last interval's qSelected
                    if (qLast <= SpreadConsX.I[lastCheckedIndex]->Vlb)
                    {
                        qSelected = QMin;
                    }
                    else if ( qLast >= SpreadConsX.I[lastCheckedIndex]->Vlb)
                    {
                        qSelected = QMax;
                    }
                    else
                    {
                        // Should never be here
                        cout << "ERROR: ENTER a place it should never be in" << endl;
                    }
                }


                // Here, qSelected is assigned for DMax and XMin
                // Time to propagate X based on R(I) or M(I)
                // get indexQSelected
                cout << "From General Mean for calculating q0 of DMax" << endl <<"qSelected: " << qSelected << endl;
                xIndex = s;
                indexQSelected = 0;
                for (int i = 0; i < SpreadConsX.numI; i++)
                {
                    if ((qSelected >= SpreadConsX.I[i]->Vlb) && (qSelected <= SpreadConsX.I[i]->Vub))
                    {
                        indexQSelected = i;
                        break;
                    }
                }
                originalIndexQSelected = indexQSelected;
                // Check R(Iq)
                for (int i = 0; i < SpreadConsX.I[originalIndexQSelected]->r; i++ )
                {
                    if (xIndex == SpreadConsX.I[originalIndexQSelected]->R[i])
                    {
    cout << "indexQSelected: " << originalIndexQSelected << endl;
    cout << " XINDEX IS IN R(I) WHERE XINDEX IS: " << xIndex << endl;
                        // Try pruning domains of X using R(I)
                        originalLowerBound = SpreadConsX.getDomVarLb(xIndex);
                        DMax =  SpreadConsX.findDMax(qSelected,originalIndexQSelected, xIndex, SpreadConsS.getDomVarUb(0));
    cout << "Dmax R(Iq) for X[" << xIndex <<"] is: " << DMax << endl;
                        SpreadConsX.setDomVarLb(xIndex, originalLowerBound); // reset original lower bound
                        maxShift = ceil(SpreadConsX.getDomVarLb(xIndex) + DMax);
                        if (maxShift < SpreadConsX.getDomVarUb(xIndex))
                        {
                            SpreadConsX.setDomVarUb(xIndex, maxShift);
                            changed = true;
                            cout << "Success: Pruned X[" << xIndex << "]'s upper bound to " <<maxShift << endl;
                        }
                        // Recreate the new intervals regardless
                        SpreadConsX.createInterval();
                        // Determine which interval qSelected is in
                        indexQSelected = 0;
                        for (int p = 0; p < SpreadConsX.numI; p++)
                        {
                            if ((qSelected >= SpreadConsX.I[p]->Vlb) && (qSelected <= SpreadConsX.I[p]->Vub))
                            {
                                indexQSelected = p;
                                break;
                            }
                        }
                        originalIndexQSelected = indexQSelected;
                        break; // don' tneed loop others
                    }
                }
                //Can skip this forloop if above executes,
                        // note: May be able to skip L[i] too  if above code executes
                        // based on Prof. Beck's logic of if q was same value
                for (int i = 0; i < SpreadConsX.I[originalIndexQSelected]->m; i++ )
                {
                    if (xIndex == SpreadConsX.I[originalIndexQSelected]->M[i])
                    {
                        cout << " INDEXQSELECTED IS " << originalIndexQSelected << endl;
                        cout << " XINDEX IS IN M(I) WHERE XINDEX IS: " << xIndex << endl;
                        v = 0;
                        // Save the originalLowerBound of Current X as it will change
                        originalLowerBound = SpreadConsX.getDomVarLb(xIndex);
                        // Since you can shift all the min to at least up to V, shift it
                        // First , get the indexQSelected without shifting
                        for (int j = 0; j < SpreadConsX.numI; j++)
                        {
                            if ((qSelected >= SpreadConsX.I[j]->Vlb) && (qSelected <= SpreadConsX.I[j]->Vub))
                            {
                                indexQSelected = j;
                                break;
                            }
                        }
                        // Calculate for v
                        v = (qSelected - SpreadConsX.I[indexQSelected]->ES)/((double) SpreadConsX.I[indexQSelected]->m);
                        DMax = v - SpreadConsX.getDomVarLb(xIndex); // can shift by at least v - current lower bound
                        // Shift current xIndex to v
                        SpreadConsX.setDomVarLb(xIndex, floor(v));
                        SpreadConsX.createInterval(); // create the new intervals
                        // Find the new intervals for indexQSelected
                        for (int j = 0; j < SpreadConsX.numI; j++)
                        {
                            if ((qSelected >= SpreadConsX.I[j]->Vlb) && (qSelected <= SpreadConsX.I[j]->Vub))
                            {
                                indexQSelected = j;
                                break;
                            }
                        }
                        DMax +=  SpreadConsX.findDMax(qSelected,indexQSelected, xIndex, SpreadConsS.getDomVarUb(0));
        cout << "Dmax M(Iq) for X[" << xIndex <<"] is: " << DMax << endl;
                        SpreadConsX.setDomVarLb(xIndex, originalLowerBound); // reset original lower bound
                        maxShift = ceil(SpreadConsX.getDomVarLb(xIndex) + DMax);
                        if (maxShift < SpreadConsX.getDomVarUb(xIndex))
                        {
                            SpreadConsX.setDomVarUb(xIndex, maxShift);
                            changed = true;
cout << "Success: Pruned X[" << xIndex << "]'s upper bound to " <<maxShift << endl;

                        }
                        // Recreate the new intervals regardless
                        SpreadConsX.createInterval();
                        // Determine which interval qSelected is in
                        indexQSelected = 0;
                        for (int p = 0; p < SpreadConsX.numI; p++)
                        {
                            if ((qSelected >= SpreadConsX.I[p]->Vlb) && (qSelected <= SpreadConsX.I[p]->Vub))
                            {
                                indexQSelected = p;
                                break;
                            }
                        }
                        originalIndexQSelected = indexQSelected;
                        break; // don't need loop others
                    }
                }

                // GET PROPER qSelect and Qindex for DMIN
                // Assuming above is done, below should be right
                for (int i = 0; i < SpreadConsX.I[originalIndexQSelected]->l; i++ )
                {
                    if (xIndex == SpreadConsX.I[originalIndexQSelected]->L[i])
                    {
                        cout << " INDEXQSELECTED IS " << originalIndexQSelected << endl;
                        cout << " XINDEX IS IN L(I) WHERE XINDEX IS: " << xIndex << endl;
                        DMin = 0;
                        originalUpperBound = 0;
                        xIndex = SpreadConsX.I[originalIndexQSelected]->L[i];
                        originalUpperBound = SpreadConsX.getDomVarUb(xIndex);
                        DMin =  SpreadConsX.findDMin(qSelected,originalIndexQSelected, xIndex, SpreadConsS.getDomVarUb(0));
                cout << "Dmin L(Iq) for X[" << xIndex <<"] is: " << DMax << endl;
                        SpreadConsX.setDomVarUb(xIndex, originalUpperBound); // reset original lower bound
                        maxShift = floor(SpreadConsX.getDomVarUb(xIndex) - DMin);
                        if (maxShift > SpreadConsX.getDomVarLb(xIndex))
                        {
                            SpreadConsX.setDomVarLb(xIndex, maxShift);
                            changed = true;
cout << "Success: Pruned X[" << xIndex << "]'s lower bound to " <<maxShift << endl;
                        }
                        // Recreate the new intervals regardless
                        SpreadConsX.createInterval();
                        // Determine which interval qSelected is in
                        indexQSelected = 0;
                        for (int p = 0; p < SpreadConsX.numI; p++)
                        {
                            if ((qSelected >= SpreadConsX.I[p]->Vlb) && (qSelected <= SpreadConsX.I[p]->Vub))
                            {
                                indexQSelected = p;
                                break;
                            }
                        }
                        originalIndexQSelected = indexQSelected;
                        break;
                    }
                }
            } // End of looping through each X variable using int s
        }
*/
//--------------------------------------------------------------------------------------------------------
// Trial 1 Code Explanation
//        cout << " a: " << a << " b: " << b << " c: " << c << " d: " << d << endl;
                        // calculate for qLast
                        // Note: Check if need to calculate for both cases + and -

                // Note: Matlab and Wolfrom Alpha gave different answers


                 // Wolfrom Alpha's Answer (may make no sense cause result in (-) values and NaN
                //        inRoot = ((4.0*a*a*c*d*d) - (a*b*b*d*d) - (4.0*a*c*d*d*d*d )+ (b*b*d*d*d*d));
                 //       qLast = ((((1.0)*sqrt(inRoot)) - (a*b) + (b*d*d))/((2.0)*(a*(a-d*d))));
                 // Matlab's Answer (makes more sense since result can always be (+) values and INFINITY
     // Note: % => Matlab Syntax for code
     //  % answer1 Code In Matlab = (d*((- b^2 + 4*a*c)*(- d^2 + a))^(1/2) + a*b - b*d^2)/(2*(- a^2 + a*d^2))
     //  % answer2 Code In Matlab = -(d*((- b^2 + 4*a*c)*(- d^2 + a))^(1/2) - a*b + b*d^2)/(2*(- a^2 + a*d^2))
     // Note: answer2 = -(answer1)
                        // inRoot = ((-1.0)*b*b + 4.0*a*c)*((-1.0)*d*d + a);
                        //   qLast = ((d* sqrt(inRoot) + a*b - b*d*d)/(2.0*((-1.0)*a*a + a*d*d)));


          //              if (a*a-a*d*d != 0) // only account for cases where qLast is not Infinity
// {
// cout << "inRoot: " << inRoot << " qLast: " << qLast << endl;
//--------------------------------------------------------------------------------------------------------
/*
 Results from Matlab

// 1. Solution for dmax(q)

// This is the positive solution (should be using this)
 -(ES - q + Xjmin*m - m*((OptMax^2*n - ES^2*n - C*n + m*q^2 - n*q^2 + q^2 + 2*ES*Xjmin*n - C*m*n + 2*ES*n*q - 2*Xjmin*n*q + OptMax^2*m*n + Xjmin^2*m*n)/(m*n))^(1/2))/(m + 1)

// This is the negative solution
 -(ES - q + Xjmin*m + m*((OptMax^2*n - ES^2*n - C*n + m*q^2 - n*q^2 + q^2 + 2*ES*Xjmin*n - C*m*n + 2*ES*n*q -
                          2*Xjmin*n*q + OptMax^2*m*n + Xjmin^2*m*n)/(m*n))^(1/2))/(m + 1)

*/
//--------------------------------------------------------------------------------------------------------
