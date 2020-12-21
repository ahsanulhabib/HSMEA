#include <malloc.h>
#include <assert.h>
#include <stdlib.h>         // qsort()
#include <math.h>           // floor()

#include "mex.h"

double *f;                      // pointer points to the objective matrix
int N, M;                       // dimensions of the input objective matrix

// this function will read the indices i and j from a and b,
// and rank a and b based on the comparison among f[a][0..M-1] and 
// f[b][0..M-1]
int compare(const void *a, const void *b)
{
    int i =  *(int*)a - 1;
    int j =  *(int*)b - 1;
    int k;  

    for (k = 0; k < M; k++)
    {
      if (f[k*N+i] < f[k*N+j])
          return -1;
      else if (f[k*N+i] > f[k*N+j])
          return 1;
    }

    return 0;     // a = b
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, j, j2;
    int kmin, kmax, k, x, count;
    int* indices;                   // original indices of    
    int** id_fronts;                // store the ids of each fronts
    int* no_ids_in_fronts;          // store number of ids in each front
    
    unsigned long NumComp = 0;      // number of comparisons 
    int NumFronts = 0;              // number of fronts
    mxArray *temp_array;            // temporary array to store the front ids of a front
    double* temp_ptr;               // pointer points to temp_array
    
    /* Check if proper number of arguments */
	if(nrhs != 1) 
    {
		mexErrMsgTxt("Only one input is required.");
	} 
    
    if(nlhs != 3) {
		mexErrMsgTxt("Three outputs are required.");
	}
    
    /* Input is the F Objective Matrix */
	N = (int)mxGetM(prhs[0]);
	M = (int)mxGetN(prhs[0]);
	if(!mxIsDouble(prhs[0]) || (N <= 0) || (M <= 0)) {
		mexErrMsgTxt("There is error in input matrix!");
	}
    
    // set the pointer to the objective matrix
    f = (double *)mxGetPr(prhs[0]);
//     for (i = 0; i < N; i++)
//     {
//         for (j = 0; j < M; j++)
//         {
//             printf("%d %f ", j*N + i, f[j*N + i]);
//         }
//         printf("\n");
//     }
    
    // initialize the indices matrix
    indices = (int *)malloc(N * sizeof(int));
    for (i = 0; i < N; i++)
        indices[i] = i+1;
    
    // sort the indices based on the first objective
    qsort (indices, N, sizeof(int), compare);
        
    // initialize the id_fronts: there are possible N fronts, so N pointers
    // are needed at the beginning
    id_fronts = (int **)malloc(N * sizeof(int*));
    no_ids_in_fronts = (int *)calloc (N, sizeof(int));
    
    // alocate memory for the first front
    id_fronts[0] = (int *)malloc(N * sizeof(int));
    id_fronts[0][0] = indices[0];       // assign the first index to the first front
    no_ids_in_fronts[0]++;
    NumFronts++;
    
	
	if (M <= 3)
	{
		// for each solution in the sorted population, do binary search
		for (i = 1; i < N; i++)
		{
			kmin = -1;                                      // the lower bound for checking
			kmax = NumFronts-1;                             // the upper bound for checking
			k = (int)floor(((double)(kmax + kmin)) / 2 + 0.5);   // the front now checked
			while (1)
			{
				// compare the objectives of solution indices[i] with those 
				// of the solutions in k(th) front
				// starting from the last solution and ending with the first one of k(th) front
				for (j = no_ids_in_fronts[k]-1; j >=0; j--)
				{
					// do non-dominated comparison
					// x = 0 means the two solutions are non-dominated
					// x = 1 means the former one dominates the latter one (but this won't happen)
					// x = 2 means the latter one dominates the former one
					x = 2;
                    count = 0;
					for (j2 = 0; j2 < M; j2++)
					{
						if ((f[j2*N + indices[i]-1] < f[j2*N + id_fronts[k][j]-1]))
						{    
							x = 0;
							break;
						}
                        else if (f[j2*N + indices[i]-1] == f[j2*N + id_fronts[k][j]-1])
                            count++;
					} // end for j2
                    if (count == M)
                    {
                        x = 0;
                        break;
                    }
					NumComp++;
					if (x == 2 || M == 2)
						break;
				} // end for j 

				if (x != 2)   // if front k has no solution dominating solution indices[i]
				{
					if (k == kmin + 1)
					{
						id_fronts[k][no_ids_in_fronts[k]] = indices[i];   // add solution indices[i] to front k
						no_ids_in_fronts[k]++;
						break;
					}
					else
					{
						kmax = k;
						k = (int)floor(((double)(kmax + kmin)) / 2 + 0.5);
					}
				}
				else // front k has at least one soultion dominating solution indices[i]
				{
					kmin = k;
					if (kmax == kmin + 1 && kmax < NumFronts-1)
					{
						// add solution indices[i] to front kmax
						id_fronts[kmax][no_ids_in_fronts[kmax]] = indices[i];
						no_ids_in_fronts[kmax]++;
						break;
					}
					else if (kmin == NumFronts-1)
					{
						// add solution indices[i] to the newest front
						id_fronts[NumFronts] = (int *)malloc(N * sizeof(int));
						id_fronts[NumFronts][0] = indices[i];
						no_ids_in_fronts[NumFronts]++;
						NumFronts++;
						break;
					}
					else
						k = (int)floor(((double)(kmax + kmin)) / 2 + 0.5);
				} 
			} // end while
		} // end for i
	}
	else
	{
		// for each solution in the sorted population, do sequential search
		for (i = 1; i < N; i++)
		{
			k = 0;	// the front is being checked
			while (1)
			{
				// compare the objectives of solution indices[i] with those 
				// of the solutions in k(th) front
				// starting from the last solution and ending with the first one of k(th) front
				for (j = no_ids_in_fronts[k]-1; j >=0; j--)
				{
					// do non-dominated comparison
					// x = 0 means the two solutions are non-dominated
					// x = 1 means the former one dominates the latter one (but this won't happen)
					// x = 2 means the latter one dominates the former one
					x = 2;
					count = 0;
					for (j2 = 0; j2 < M; j2++)
					{
						if ((f[j2*N + indices[i]-1] < f[j2*N + id_fronts[k][j]-1]))
						{    
							x = 0;
							break;
						}
                        else if (f[j2*N + indices[i]-1] == f[j2*N + id_fronts[k][j]-1])
                            count++;
					} // end for j2
                    if (count == M)
                    {
                        x = 0;
                        break;
                    }
					NumComp++;
					if (x == 2 || M == 2)
						break;
				} // end for j
				
				if (x != 2)   // if front k has no solution dominating solution indices[i]
				{
					id_fronts[k][no_ids_in_fronts[k]] = indices[i];   // add solution indices[i] to front k
					no_ids_in_fronts[k]++;
					break;
				}
				else
				{
					if (k < NumFronts-1)
						k++;
					else
					{
						// add solution indices[i] to the newest front
						id_fronts[NumFronts] = (int *)malloc(N * sizeof(int));
						id_fronts[NumFronts][0] = indices[i];
						no_ids_in_fronts[NumFronts]++;
						NumFronts++;
						break;
					}
				}
			} // end while
		} // end for i
	}
   
    // prepare the outputs
    plhs[0] = mxCreateCellMatrix(NumFronts,1);
    plhs[1] = mxCreateDoubleScalar((double) NumComp);
    plhs[2] = mxCreateDoubleScalar((double) NumFronts);
    
    for (i = 0; i < NumFronts; i++)
    {
        temp_array = mxCreateDoubleMatrix(no_ids_in_fronts[i], 1, mxREAL);
        temp_ptr = (double *)mxGetPr(temp_array);
        
        for (j = 0; j < no_ids_in_fronts[i]; j++)
            temp_ptr[j] = (double)id_fronts[i][j];
        
        mxSetCell(plhs[0], i, mxDuplicateArray(temp_array));
        mxDestroyArray(temp_array); // release memory
        free(id_fronts[i]);         // release memory
    }
    
    // release the memory
    free(indices);
    free(no_ids_in_fronts);
    free(id_fronts);
}