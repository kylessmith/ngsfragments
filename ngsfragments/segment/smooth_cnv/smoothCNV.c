#include <stdio.h>
#include <stdlib.h>
#include "smoothCNV.h"


static int cmpfunc(const void * a, const void * b)
{
	if (*(double*)a > *(double*)b)
	{
		return 1;
	}
	else if (*(double*)a < *(double*)b)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}

void smoothRL(int n, double *gdat, int nchr, long *cfrq, double *sgdat, int k, double oSD, double sSD)
{
	int i, j, ilo, ihi, k1, k2, j1, ic, cilo, cihi;
	double mnnbd, mxnbd, distij, xmed;
	
    k2 = 2*k + 1;
	// temporary array for finding median
	//double xnbhd[k2];
	double *xnbhd = (double *)malloc(sizeof(double)*k2);
	
	//initial values for start and end of chromosomes
	cilo = 0;
	cihi = 0;
	
	//loop over chromomsomes
	for (ic=0; ic<nchr; ic++)
	{
		//end of the current chromosome
		cihi = cihi + cfrq[ic];
		for (i=cilo; i<cihi; i++)
		{
			//range of the neighborhood
			ilo = MAX(cilo, i-k);
			ihi = MIN(cihi, i+k);
			
			//check if ith observation is an outlier
			//initialize the distances to be large
			mxnbd = 100 * oSD;
			mnnbd = 100 * oSD;
			for (j=ilo; j<ihi; j++)
			{
				if (j != i)
				{
					//calculate distance from between ith and jth obsn
					distij = gdat[i] - gdat[j];
					
					//if distance is less than oSD no smoothing necessary
					if (fabs(distij) <= oSD)
					{
						sgdat[i] = gdat[i];
						break;
					}//otherwise calculate distances from above and below
					else
					{
						//mxnbd is distance from above
						if (distij < mxnbd) { mxnbd = distij; }
						//mnnbd is distance from below
						if (-distij < mnnbd) { mnnbd = -distij; }
					}
				}
			}

			//if distance is less than oSD no smoothing necessary
			if (fabs(distij) <= oSD)
			{
				continue;
			}

			//if all the points in the nbhd are above than mxnbd will be negative
			//and mnnbd will be greater than oSD. Vice versa if all points below
			// 
			//If both are negative then the ith point is singleton in the middle 
			//but distance oSD away from all points in the nbhd. No smoothing done.
			
			if ((mxnbd <= 0) && (mnnbd <= 0))
			{
				sgdat[i] = gdat[i];
			}
			else
			{
				//calculate the median of the nbhd
				//number of points in the nbhd
				k1 = ihi - ilo + 1;
				//get the data into temporary array
				for (j=ilo; j<ihi; j++)
				{
					xnbhd[j-ilo+1] = gdat[j];
				}
				
				//sort the data
				//qsort(xnbhd, k1, sizeof(xnbhd[0]), cmpfunc);
				qsort(xnbhd, k1, sizeof(xnbhd[0]), cmpfunc);
				
				//median is the midpoint if n is odd and average of the two if even
				j1 = k1 / 2;
				if (k1 == (2*j1))
				{
					xmed = (xnbhd[j1] + xnbhd[j1+1]) / 2;
				}
				else
				{
					xmed = xnbhd[j1+1];
				}
				
				//if point is above the nbhd bring it down
				if (mxnbd > 0)
				{
					sgdat[i] = xmed + sSD;
				}
				//if point is below the nbhd bring it up
				if (mnnbd > 0)
				{
					sgdat[i] = xmed - sSD;
				}
			}
		}
		//beginning of next chromosome
		cilo = cilo + cfrq[ic];
	}
	
	free(xnbhd);
	
	return;
}