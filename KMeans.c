/*
 * KMeans algorithm with openmp and c
 */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <float.h>

int numPoints;
int numCenters;
int dimension;
int numIteration;
int bindThreads;
char *pointFile;
char *centerFile;
char *resultsFile;


int main(int argc, int **argv[]){

    int numproc, myid;
    int ret = parse_args(argc,argv);

    int localNumPoints;
    double *centers = malloc(sizeof(double) * numCenters * dimension);

#pragma omp parallel private(myid)
    {
        numproc = omp_get_num_threads();
        myid = omp_get_thread_num();
        localNumPoints = numPoints/numproc;
        double *localPoints =  malloc(sizeof(double) * localNumPoints * dimension);

        if(myid == 0){
            printf("Reading centers");
            FILE *c = fopen(centerFile, "rb");
            fread(centers, sizeof(double), numCenters * dimension, c);
            fclose(c);
        }

        int startIdx = myid*localNumPoints;
        FILE *f = fopen(pointFile, "rb");
        printf("Reading points %d %d\n", myid, startIdx);

        fseek(f, startIdx * dimension * sizeof(double), SEEK_SET);
        fread(localPoints, sizeof(double), localNumPoints * dimension, f);


        int i;
        for(i = 0; i < localNumPoints; ++i){
            int points_offset = i * dimension;
            
        }
    }


    return 0;
}

double euclidean_distance(double *points1, double* points2, int offset1,
                          int offset2, int dim) {
    double d = 0.0;
    double tmp;
    int i;
    for (i = 0; i < dim; ++i) {
        tmp = points1[i + offset1] - points2[i + offset2];
        d += tmp * tmp;
    }
    return sqrt(d);
}

int find_nearest_center(double *points, double *centers, int num_centers,
                         int dim, int points_offset) {
    double min_dist = DBL_MAX;
    int min_dist_idx = -1;
    int i;
    for (i = 0; i < num_centers; ++i) {
        double dist = euclidean_distance(points, centers, points_offset,
                                         i * dim, dim);
        if (dist < min_dist) {
            min_dist = dist;
            min_dist_idx = i;
        }
    }
    return min_dist_idx;
}

int parse_args(int argc, char **argv) {
    int index;
    int c;

    opterr = 0;
    while ((c = getopt(argc, argv, "n:d:k:m:o:c:p:b:v")) != -1)
        switch (c) {
            case 'n':
                numPoints = atoi(optarg);
                break;
            case 'd':
                dimension = atoi(optarg);
                break;
            case 'k':
                numCenters = atoi(optarg);
                break;
            case 'm':
                numIteration = atoi(optarg);
                break;
            case 'o':
                resultsFile = optarg;
                break;
            case 'c':
                centerFile   = optarg;
                break;
            case 'p':
                pointFile = optarg;
                break;
            case 'b':
                bindThreads = atoi(optarg);
                break;
            default:
                abort();
        }
        printf("Program Arguments\n");
        printf(
                " n = %d\n d = %d\n k = %d\n m = %d\n o = %s\n c = %s\n p = %s\n b = %d\n",
                numPoints, dimension, numCenters, numIteration, resultsFile, centerFile, pointFile, bindThreads);

    for (index = optind; index < argc; index++)
        printf("Non-option argument %s\n", argv[index]);
    return 0;
}