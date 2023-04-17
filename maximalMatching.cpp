#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include "maximalMatching.h"


// helper function used in Serial Karp-Sipser initialization
void findMateS(long u, graph* G, long* flag,long* mate, long* degree)
{
	if(flag[u] != 0) return;
    flag[u] = 1;
	long *endVertex = G->endV;
	long *edgeStart = G->vtx_pointer;
	
	long neighbor_first = edgeStart[u];
	long neighbor_last = edgeStart[u+1];
	for(long j=neighbor_first; j<neighbor_last; j++)
	{
		long v = endVertex[j];
		if(flag[v] == 0) // if I can lock then this v node is unmatched
		{
            flag[v] = 1;
			mate[u] = v;
			mate[v] = u;
			// update degree
			long neighborFirstU = edgeStart[v];
			long neighborLastU = edgeStart[v+1];
			for(long k=neighborFirstU; k< neighborLastU; k++)
			{
				long nextU = endVertex[k];
                degree[nextU]--;
				if( degree[nextU] == 1)
				{
					
					findMateS(nextU,G,flag,mate,degree);
                    
				}
				
			}
			break;
		}
	}
}

// Serial Karp-Sipser maximal matching
long KarpSipserInitS(graph* G, long* unmatchedU,  long* mate)
{
	long nrows = G->nrows;
	long * degree = (long*) malloc(sizeof(long) * nrows);
	long* degree1Vtx = (long*) malloc(sizeof(long) * nrows);
	
	long *endVertex = G->endV;
	long *edgeStart = G->vtx_pointer;
	long nrowsV = G->n;
	long numUnmatchedU = 0;
	long* flag = (long*) malloc(sizeof(long) * nrowsV);
	
	double timeStart = omp_get_wtime();
	
	for(long i=0; i< nrowsV; i++)
	{
		flag[i] = 0;
		mate[i] = -1;
	}
	
	long degree1Tail = 0;
	long degree1Count = 0;
	
	
	
	for(long u=0; u<nrows; u++)
	{
		degree[u] = edgeStart[u+1] - edgeStart[u];
		if(degree[u] == 1)
		{
			degree1Vtx[degree1Count++] = u;
		}
	}
	
	
	
	for(long u=0; u<degree1Count; u++)
	{
            findMateS(degree1Vtx[u],G,flag,mate,degree);
	}
	
	
	for(long u=0; u<nrows; u++)
	{
		if(flag[u] == 0 && degree[u]>0)
			findMateS(u,G,flag,mate,degree);
	}
	
	double timeInit = omp_get_wtime()-timeStart;
	for(long u=0; u<nrows; u++)
	{
		
		if(mate[u] == -1 && (edgeStart[u+1] > edgeStart[u]))
		{
			unmatchedU[numUnmatchedU++] = u;
		}
	}
    
    long matched_rows = 0;
    for(long u=0; u<nrows; u++)
    {
        if(mate[u]!=-1) matched_rows++;
    }

	
    printf("===========================================\n");
    printf("Serial Karp-Sipser Initialization\n");
    printf("===========================================\n");
    printf("Matched Rows        = %ld (%.2lf%%)\n", matched_rows, (double)100.0*(matched_rows)/nrows);
    printf("Computation time    = %lf\n", timeInit);
    printf("===========================================\n");
	//printf("%lf %lf \n", 100.0 * (nrows - numUnmatchedU)/nrows, timeInit);
	free(degree1Vtx);
	free(degree);
	free(flag);
	return numUnmatchedU;
}



// helper function used in Karp-Sipser initialization 
void findMate(long u, graph* G, long* flag,long* mate, long* degree)
{
	if(__sync_fetch_and_add(&flag[u],1) != 0) return;
	long *endVertex = G->endV;
	long *edgeStart = G->vtx_pointer;
	
	long neighbor_first = edgeStart[u];
	long neighbor_last = edgeStart[u+1];   
	for(long j=neighbor_first; j<neighbor_last; j++) 
	{
		long v = endVertex[j];
		if(__sync_fetch_and_add(&flag[v],1) == 0) // if I can lock then this v node is unmatched
		{
			mate[u] = v;
			mate[v] = u;
			// update degree
			long neighborFirstU = edgeStart[v];
			long neighborLastU = edgeStart[v+1];
			for(long k=neighborFirstU; k< neighborLastU; k++)
			{
				long nextU = endVertex[k];
				if( __sync_fetch_and_add(&degree[nextU],-1) == 2)
				{
					
					findMate(nextU,G,flag,mate,degree);
				}
				
			}
			break;
		} 
	}
}


// Multithreaded  Karp-Sipser maximal matching
long KarpSipserInit(graph* G, long* unmatchedU,  long* mate)
{
	long nrows = G->nrows;
	long * degree = (long*) malloc(sizeof(long) * nrows);
	long* degree1Vtx = (long*) malloc(sizeof(long) * nrows);
	
	long *endVertex = G->endV;
	long *edgeStart = G->vtx_pointer;
	long nrowsV = G->n;
	long numUnmatchedU = 0;
	long* flag = (long*) malloc(sizeof(long) * nrowsV);
	
	double timeStart = omp_get_wtime();
	
#pragma omp parallel for default(shared) schedule(static)
	for(long i=0; i< nrowsV; i++)
	{
		flag[i] = 0;
		mate[i] = -1;
	}
	
	long degree1Tail = 0;
	long degree1Count = 0;  
	
	
	
	// populate degree and degree1Vtx
#pragma omp parallel for default(shared) schedule(static)//schedule(dynamic)
	for(long u=0; u<nrows; u++)
	{      
		degree[u] = edgeStart[u+1] - edgeStart[u];
		if(degree[u] == 1)
		{
			degree1Vtx[__sync_fetch_and_add(&degree1Count,1)] = u;
			//flag[u] = 1; // means already taken 
		}
	}
	
	
	
#pragma omp parallel for default(shared) //schedule(dynamic,100)
	for(long u=0; u<degree1Count; u++)
	{
		//findMate1(degree1Vtx[u],G,flag,mate,degree,degree1Vtx,&degree1Head);
		findMate(degree1Vtx[u],G,flag,mate,degree);		  
	}
	
	
	// process other vertices 
#pragma omp parallel for default(shared) schedule(dynamic,100)//schedule(dynamic)
	for(long u=0; u<nrows; u++)
	{
		if(flag[u] == 0 && degree[u]>0)
			findMate(u,G,flag,mate,degree);	  
	}
	
	double timeInit = omp_get_wtime()-timeStart;
#pragma omp parallel for default(shared) 
	for(long u=0; u<nrows; u++)
	{
		
		if(mate[u] == -1 && (edgeStart[u+1] > edgeStart[u]))
		{
			unmatchedU[__sync_fetch_and_add(&numUnmatchedU, 1)] = u;
		}
	}
    
    long matched_rows = 0;
#pragma omp parallel
    {
        long tmatched = 0; //thread private variables
#pragma omp for
        for(long u=0; u<nrows; u++)
        {
            if(mate[u]!=-1) tmatched++;
        }
        __sync_fetch_and_add(&matched_rows,tmatched);
    }

	
    printf("===========================================\n");
    printf("Karp-Sipser Initialization\n");
    printf("===========================================\n");
    printf("Matched Rows        = %ld (%.2lf%%)\n", matched_rows, (double)100.0*(matched_rows)/nrows);
    printf("Computation time    = %lf\n", timeInit);
    //printf("%lf %lf", 100.0 * (nrows - numUnmatchedU)/nrows, timeInit);
    printf("===========================================\n");
	
	free(degree1Vtx);
	free(degree);
	free(flag);
	return numUnmatchedU;
}







