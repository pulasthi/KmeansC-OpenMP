/*
 * KMeans algorithm with openmp and c
 */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <jmorecfg.h>


int numPoints;
int numCenters;
int dimension;
int numIteration;
int bindThreads;
char *pointFile;
char *centerFile;
char *resultsFile;

double euclidean_distance(double *points1, double* points2, int offset1,
                          int offset2, int dim);
void addToSum(double *points, int pointOffset, double *sums, int nearestCenter, int dimension);
void resetToZero(double *array, int length);
void resetToZeroInt(int *array, int length);

int main(int argc, int **argv[]){
    boolean PrintResults = TRUE;

    int numproc, myid;
    int ret = parse_args(argc,argv);
    printf("######################################## Starting Kmeans ########################################\n");
    clock_t start = clock(), diff;
    int localNumPoints;
    double *centers = malloc(sizeof(double) * numCenters * dimension);

    int *centerCountsTotal = malloc((sizeof(int))*numCenters);
    double *centerSumsTotal = malloc((sizeof(double))*numCenters*dimension);

    resetToZero(centerSumsTotal,numCenters*dimension);
    resetToZeroInt(centerCountsTotal,numCenters);

#pragma omp parallel private(myid)
    {
        numproc = omp_get_num_threads();
        myid = omp_get_thread_num();
        localNumPoints = numPoints/numproc;
        double *localPoints =  malloc(sizeof(double) * localNumPoints * dimension);
        int *centerCounts = malloc((sizeof(int))*numCenters);
        double *centerSums = malloc((sizeof(double))*numCenters*dimension);
        resetToZero(centerSums,numCenters*dimension);
        resetToZeroInt(centerCounts,numCenters);

        if(myid == 0){
            printf("Reading centers \n");
            FILE *c = fopen(centerFile, "rb");
            fread(centers, sizeof(double), numCenters * dimension, c);
            fclose(c);
        }

        int startIdx = myid*localNumPoints;
        FILE *f = fopen(pointFile, "rb");
        printf("Reading points %d %d\n", myid, startIdx);

        fseek(f, startIdx * dimension * sizeof(double), SEEK_SET);
        fread(localPoints, sizeof(double), localNumPoints * dimension, f);
        fclose(f);

        while(numIteration > 0){
            int i;
            for(i = 0; i < localNumPoints; ++i){
                int points_offset = i * dimension;
                int nearest_center = find_nearest_center(localPoints, centers, numCenters,
                                                         dimension, points_offset);
                ++centerCounts[nearest_center];
                addToSum(localPoints, points_offset, centerSums, nearest_center, dimension);

            }


            #pragma omp critical
            {
                int j;
                for (j = 0; j < numCenters; ++j) {
                    centerCountsTotal[j] += centerCounts[j];
                    int k;
                    for (k = 0; k < dimension; ++k) {
                        centerSumsTotal[j*dimension + k] += centerSums[j*dimension + k];
                    }
                }
            }


            #pragma omp barrier
            if(myid == 0){
                numIteration -= 1;
                int j;
                for (j = 0; j < numCenters; ++j) {
                    int k;
                    for (k = 0; k < dimension; ++k) {
                        centers[j*dimension + k] = centerSumsTotal[j*dimension + k]/centerCountsTotal[j];
                    }
                }
                resetToZero(centerSumsTotal,numCenters*dimension);
                resetToZeroInt(centerCountsTotal,numCenters);
            }



            #pragma omp barrier
            resetToZero(centerSums,numCenters*dimension);
            resetToZeroInt(centerCounts,numCenters);
        }

    }

    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    int i;
    int sum = 0;

    for(i = 0; i < numCenters; ++i){
        printf("centers %d x  y z value %f %f %f \n",i,centers[i*dimension],centers[i*dimension+1],centers[i*dimension+2]);
    }
    printf("######################################## Ending Kmeans ########################################\n");

    if(PrintResults == TRUE){
        int localNumPoints = numPoints;
        double *localPoints =  malloc(sizeof(double) * numPoints * dimension);
        int startIdx = 0;
        FILE *fall = fopen(pointFile, "rb");
        FILE *fout;
        fout = fopen("/home/pulasthi/ClionProjects/hpc/assigments/hpcproject/data/ouput.txt", "w+");


        printf("Reading points %d %d\n", myid, startIdx);

        fseek(fall, startIdx * dimension * sizeof(double), SEEK_SET);
        fread(localPoints, sizeof(double), localNumPoints * dimension, fall);
        fclose(fall);

        for(i = 0; i < localNumPoints; ++i){
            int points_offset = i * dimension;
            int nearest_center = find_nearest_center(localPoints, centers, numCenters,
                                                     dimension, points_offset);

            fprintf(fout, "%d %f %f %f %d %d\n",i, localPoints[i*dimension],localPoints[i*dimension + 1], localPoints[i*dimension + 2], nearest_center ,nearest_center);
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

void addToSum(double *points, int pointOffset, double *sums, int nearestCenter, int dimension){
    int i;
    for (i = 0; i < dimension; ++i) {
        sums[nearestCenter*dimension + i] += points[pointOffset + i];
    }
}
void resetToZero(double *array, int length) {
    int i;
    for (i = 0; i < length; ++i) {
        array[i] = 0.0;
    }
}

void resetToZeroInt(int *array, int length) {
    int i;
    for (i = 0; i < length; ++i) {
        array[i] = 0;
    }
}

int find_nearest_center(double *points, double *centers, int num_centers,
                         int dim, int points_offset) {
    double min_dist = DBL_MAX;
    int min_dist_idx = -1;
    int i;
    for (i = 0; i < num_centers; ++i) {
        double dist = euclidean_distance(points, centers, points_offset,
                                         i * dimension, dimension);
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