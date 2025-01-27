/*=================================================================
 *
 * getZNICIM.C
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"

/* Input Arguments */

#define Y_IN0 prhs[0]
#define Y_IN1 prhs[1]




/* Output Arguments */

#define znics_OUT plhs[0]

int ifexists(double *z,  int u, double v[], int W)
{
    int i, j;
    for (j=0; j<W;j++)
        for (i=0; i<u;i++)
            if (*(z+i)==v[j]) return (1);
    return (0);
}



static void getZNICIM(double znics[], double nodes[], int intW, int m) {
		
	  
	double logNchoosen ;
	double nodesi[intW];
	double *unqNodes =(double*) mxCalloc(m, sizeof(double));
    int i,k, i1;
   for(i1=0;i1<intW;i1++ )
      *(unqNodes+i1)=nodes[i1];

    k=intW;
	znics[0]=1;
    for(i1=1;i1<m;i1++ )
         znics[i1]=0;	
	
    for (i=intW;i<m;i+=intW)
    {
		for(i1=0;i1<intW;i1++ )
            nodesi[i1]=nodes[i+i1];		
        if(!ifexists(unqNodes, k,nodesi, intW))
        {
            for(i1=0;i1<intW;i1++ )
                *(unqNodes+k+i1)=nodes[i+i1];
   
            k+=intW;
			znics[i]=1;

        }

    }
	

    return;
}
 
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])

{

    double* znics;
    double  *nodes,  *y2, *y3;
    size_t m, intW;
    m = mxGetM(Y_IN0);


    /* Create a matrix for the return argument */
    znics_OUT = mxCreateDoubleMatrix((mwSize)m, (mwSize)1, mxREAL);

    /* Assign pointers to the various parameters */
    znics = mxGetPr(znics_OUT);
	


    nodes = mxGetPr(prhs[0]);
    y2 = mxGetPr(prhs[1]);
	intW= mxGetM(Y_IN1);
	
    /* Do the actual computations in a subroutine */
    getZNICIM(znics, nodes,intW, m);
    return;
}

/* LocalWords:  znics maxlhs
 */

