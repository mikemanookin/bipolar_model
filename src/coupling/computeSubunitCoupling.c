/*=================================================================
 *
 * COMPUTESUBUNITCOUPLING.C	Sample .MEX file corresponding to YPRIME.M
 *	        Solves simple 3 body orbit problem
 *
 * The calling syntax is:
 *
 *		out = computeSubunitCoupling(currents, exp(-sDistances/lambda)*gain, numSubunits, numT)
 *
 *
 * This is a MEX-file for MATLAB.
 * Written by Mike Manookin 3/2/2021
 *
 *=================================================================*/

#include <math.h>
#include "mex.h"

/* Input Arguments */
#define CURRENTS prhs[0]
#define DISTANCES prhs[1]

/* Output Arguments */
#define YP_OUT plhs[0]

static void computeCoupling(double* x, double yp[], double* distances, int numSubunits, int numT) {
    int ns, nt, j;
    double currentDelta, drivingForce;
    
    for (nt = 0; nt < numT; nt++) {
        /* Loop through the subunits. */
        for (ns = 0; ns < numSubunits; ns++) {
            currentDelta = 0;
            for (j = 0; j < numSubunits; j++) {
                if (ns != j) {
                    drivingForce = x[nt*numSubunits+j] - x[nt*numSubunits+ns];
                    currentDelta += drivingForce * distances[ns*numSubunits+j];
                }
            }
            yp[nt*numSubunits+ns] = x[nt*numSubunits+ns] + currentDelta;
        }
    }
    
    return;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    double* yp;
    double *distances, *currents;
    
    size_t m, n;

    /* Check the dimensions of the input currents. */
    m = mxGetM(CURRENTS);
    n = mxGetN(CURRENTS);

    /* Create a matrix for the return argument */
    YP_OUT = mxCreateDoubleMatrix((mwSize)m, (mwSize)n, mxREAL);

    /* Assign pointers to the various parameters */
    yp = mxGetPr(YP_OUT);
    distances = mxGetPr(DISTANCES);
    currents = mxGetPr(CURRENTS);
    
    /* Do the actual computations in a subroutine */
    computeCoupling(currents, yp, distances, (int)m, (int)n);
    
    return;
}


