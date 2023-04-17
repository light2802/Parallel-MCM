#include <stdio.h>



// check the mate array and return maching size
// nrows  - number of u vertices (row vertice )
long checkMatching(long nrows,long* mate)
{
    long count = 0;
    long umSize = 0;
#pragma omp parallel for default(shared)
    for(long i=0; i<nrows; i++)
    {
        if(mate[i]!=-1)
        {
            if(mate[mate[i]]!= i)
            {
                __sync_fetch_and_add(&count,1);
            }
        }
        else
            __sync_fetch_and_add(&umSize,1);
        
    }
    printf("Matching Size = %d\n", nrows - umSize);
    if(count !=0)
        printf("!!!!!!!!!!!!!error count : %d\n",count);
    return nrows-umSize;
}


