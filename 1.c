#include <stdio.h>
#include <omp.h>

int main ()
{
    const int n = 30;
    int   i,chunk;
    double a[n], b[n], result = 0.0;
    double results[10];

    /* Some initializations */
    chunk = 5;
    for (i=0; i < n; i++) {
        a[i] = i * 3.14;
        b[i] = i * 6.67;
        results[i] = 0.0;
    }

#pragma omp parallel for default(shared) private(i) schedule(static,chunk) reduction(+:results[:10])
    for (i=0; i < n; i++){
        int j;
        for (j=0; j < n; j++){
            ++results[j];
        }
    }

    printf("Final result= %f\n",result);
}