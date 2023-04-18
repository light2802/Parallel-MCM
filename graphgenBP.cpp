#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <alloca.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iterator>
#include "graphgenBP.h"
#include "ThreadedMMReader.h"

using namespace std;




void free_graph( graph* bGraph)
{
	delete [] bGraph->vtx_pointer;
	delete [] bGraph->endV;
    if(bGraph->weight) delete [] bGraph->weight;
}


graph* swap_side(graph* bGraph)
{
    long* tvtx_pointer = new long[bGraph->n+1];
    double* tweight = new double[bGraph->m];
    long* tendV = new long[bGraph->m];
    bGraph->vtx_pointer[0]=0;
    
    long last_rp = bGraph->vtx_pointer[bGraph->nrows];
    for(long i=0, j=last_rp; i<last_rp; i++, j++)
    {
        tendV[j] = bGraph->endV[i] - bGraph->nrows;
        tweight[j] = bGraph->weight[i];
    }
    for(long i=0, j=last_rp; i<last_rp; i++, j++)
    {
        tendV[i] = bGraph->endV[j] + bGraph->nrows;
        tweight[j] = bGraph->weight[i];
    }
    
    for(long i=0, j=bGraph->nrows; i<=bGraph->nrows; i++, j++)
    {
        tvtx_pointer[j] = bGraph->vtx_pointer[i] + last_rp;
    }
    for(long i=0, j=bGraph->nrows; i<bGraph->nrows; i++, j++)
    {
        tvtx_pointer[i] = bGraph->vtx_pointer[j] - last_rp;
    }

    graph* tGraph = (graph *) malloc(sizeof(graph));
    tGraph->m = bGraph->m;
    tGraph->n = bGraph->n;
    tGraph->nrows = bGraph->n - bGraph->nrows;
    tGraph->vtx_pointer = tvtx_pointer;
    tGraph->weight = tweight;
    tGraph->endV = tendV;
    return tGraph;
}

bool isEqual(graph* bGraph1, graph* bGraph2)
{
    if(bGraph1->n != bGraph2->n)
    {
        cout << "bGraph1->n != bGraph2->n \n";
        return false;
    }
    else if(bGraph1->nrows != bGraph2->nrows)
    {
        cout << "bGraph1->nrows != bGraph2->nrows \n";
        return false;
    }
    else if(bGraph1->m != bGraph2->m)
    {
        cout << "bGraph1->m != bGraph2->m \n";
        return false;
    }
    
    for(long i=0; i<=bGraph1->n ;i++)
    {
        if(bGraph1->vtx_pointer[i] != bGraph2->vtx_pointer[i])
        {
            cout << i << " : bGraph1->vtx_pointer[i] != bGraph2->vtx_pointer[i] \n" ;
            return false;
        }
    }
    
    for(long i=0; i<bGraph1->m ;i++)
    {
        if(bGraph1->endV[i] != bGraph2->endV[i])
        {
            cout << i << " : " << bGraph1->endV[i] << " & " << bGraph2->endV[i] << " : bGraph1->endV[i] != bGraph2->endV[i] \n";
            return false;
        }
    }

    
    cout <<  "******* The graphs are equal **********\n";
    return true;
    
}

/*
 Multithreaded prefix sum
 Inputs:
    in: an input array
    size: the length of the input array "in"
    nthreads: number of threads used to compute the prefix sum
 
 Output:
    return an array of size "size+1"
    the memory of the output array is allocated internallay 
 
 Example:
 
    in = [2, 1, 3, 5]
    out = [0, 2, 3, 6, 11]
 */
template <typename T>
T* prefixsum(T* in, int size, int nthreads)
{
    vector<T> tsum(nthreads+1);
    tsum[0] = 0;
    T* out = new T[size+1];
    out[0] = 0;
    T* psum = &out[1];
    
#pragma omp parallel
    {
        int ithread = omp_get_thread_num();
        T sum = 0;
#pragma omp for schedule(static)
        for (int i=0; i<size; i++)
        {
            sum += in[i];
            psum[i] = sum;
        }
        
        tsum[ithread+1] = sum;
#pragma omp barrier
        T offset = 0;
        for(int i=0; i<(ithread+1); i++)
        {
            offset += tsum[i];
        }
#pragma omp for schedule(static)
        for (int i=0; i<size; i++)
        {
            psum[i] += offset;
        }
        
    }
    return out;
}




void fast_mtx_read_build(char *fname, graph* bGraph)
{
    vector<long> allrows;
    vector<long> allcols;
    vector<double> allvals;
    long nrows, ncols;
    bool isWeighted;
    
    double time_start = omp_get_wtime();
    ThreadedMMReader(fname, allrows, allcols, allvals, nrows, ncols, isWeighted);
    cout << "ThreadedMMReader read in " << omp_get_wtime() - time_start << "  seconds"<< endl;
    
    
    time_start = omp_get_wtime();
    
    bGraph->m = allrows.size() * 2; // save edges in both direction
    bGraph->n = nrows + ncols;
    bGraph->nrows = nrows;
    // This will be allocated in the prefixsum function
    //bGraph->vtx_pointer = new long[bGraph->n+1];
    bGraph->weight = new double[bGraph->m];
    bGraph->endV = new long[bGraph->m];
    //bGraph->vtx_pointer[0]=0;
    long* degrees = new long[bGraph->n];
    
    
#pragma omp parallel for
    for(long i =0; i<bGraph->n; i++)
    {
        degrees[i] = 0;
    }

#pragma omp parallel for
    for(long i =0; i< allrows.size(); i++)
    {
        __sync_fetch_and_add(degrees + allrows[i],1);
        __sync_fetch_and_add(degrees + nrows + allcols[i],1);
    }
    
    long nthreads=1;
#pragma omp parallel
    {
        nthreads = omp_get_num_threads();
    }
    
    bGraph->vtx_pointer = prefixsum(degrees, bGraph->n, nthreads);
    
    long* cur_vtx_pointer = degrees; // just reuse the same memory
#pragma omp parallel for
    for(long i =0; i<bGraph->n; i++)
    {
        cur_vtx_pointer[i] = bGraph->vtx_pointer[i];
    }
    
    
#pragma omp parallel for
    for(long i =0; i<allrows.size(); i++)
    {
        long rowIdx = __sync_fetch_and_add(cur_vtx_pointer + allrows[i],1);
        long colIdx = __sync_fetch_and_add(cur_vtx_pointer + nrows + allcols[i],1);
        //if(rowIdx < bGraph->vtx_pointer[allrows[i]+1] && colIdx < bGraph->vtx_pointer[allcols[i]+nrows+1])
        {
            bGraph->endV[rowIdx] = nrows + allcols[i];
            bGraph->endV[colIdx] = allrows[i];
            if(isWeighted)
            {
                bGraph->weight[rowIdx] = allvals[i];
                bGraph->weight[colIdx] = allvals[i];
            }
            
        }
        /*
        else
        {
            cout << "Error in creating graph from matrix market file\n";
            exit(1);
        }*/
    }
  
    cout << "Graph generated in " << omp_get_wtime() - time_start << "  seconds"<< endl;
    delete [] degrees;
}
