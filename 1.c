#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[])
{
    const int N = 20;
    int nthreads, threadid, i;
    double a[N], b[N], c[N];

    // Initialize
    for (i=0; i < N; i++) {
        a[i] = 1.0*i;
        b[i] = 2.0*i;
    }

    int chunk = 3;

#pragma omp parallel shared(a,b,c,nthreads,chunk) private(i,threadid)
    {
        threadid = omp_get_thread_num();
        if (threadid == 0) {
            nthreads = omp_get_num_threads();
            printf("Number of threads = %d\n", nthreads);
        }
        printf(" My threadid %d\n",threadid);

#pragma omp for schedule(static,chunk)
        for (i=0; i<N; i++) {
            c[i] = a[i] + b[i];
            printf(" Thread id: %d working on index %d\n",threadid,i);
        }

    } // join

    printf(" TEST c[19] = %g\n",c[19]);

    return 0;
}

