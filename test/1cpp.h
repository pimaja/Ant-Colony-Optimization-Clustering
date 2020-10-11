
#include <iostream>
#include <time.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
//#include <conio.h>
#include <limits.h>
#include <iomanip>
#include <algorithm>
#include <list>
#include <iterator>
#include <set>
#include <math.h>
#include <string.h>

using namespace std;

class klasa
{
    int R, maxIterations;
    string s;
    double **pointerMatrix;

public:
    klasa( int r, int mi, string podaci) { R=r; maxIterations=mi; s=podaci;};

    int clusterWithMaxPheromoneTrail(double** F, int element, int K);

    int clusterChosenStochastically(double** F, int element, int K);

    void generateSolution(int** S, double** F, int agent, int N, int K, int** W);

    void centersOfClusters(int** W, double** attributes, double** centers, int N, int K, int n);

    double calculateFitness(int** W, double** attributes, double** centers, int N, int K, int n);

    void findBestL(double* fitness, double percent , int R, int* bestInd);

    void localSearch(int** S, int N, int* best, int L, int K, double* fitness, double** attributes, int n);

    void evaporatePheromone(double** F, int N, int K);

    void addPheromone(int* best, int L, double** F, int N, int K, int** S, double* fitness);

    int main2(double &bfit, double &kfit);

    friend class MainWindow;

};

