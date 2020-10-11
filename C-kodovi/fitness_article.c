#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include <errno.h>

void printMatrixInt(int** M, int m, int n){
    int i,j;
    for (i = 0; i <  m; i++){ //print matrix
        for (j = 0; j < n; j++)
            printf("%d ", M[i][j]);
        printf("\n");
    }
    printf("\n");
    return;
}

void printMatrixDouble(double** M, int m, int n){
    int i,j;
    for (i = 0; i <  m; i++){ //print matrix
        for (j = 0; j < n; j++)
            printf("%f ", M[i][j]);
        printf("\n");
    }
    printf("\n");
    return;
}




void centersOfClusters(int** W, double** attributes, double** centers, int N, int K, int n){
    int i, j, v;
    //printMatrixDouble(attributes, N, n);
    for(j=0; j<K; j++)
        for(v=0; v<n; v++){
            double top=0.0, bottom=0.0;
            for (i=0; i<N; i++)
                if(W[i][j]){
                    top += (W[i][j]*attributes[i][v]);
                    bottom += W[i][j];
                }
            centers[j][v] = (top/bottom);
        }
}

double calculateFitness(int** W, double** attributes, double** centers, int N, int K, int n){
    int i, j, v;
    double fitness = 0.0;
    //printMatrixDouble(attributes, N, n);
    for(j=0; j<K; j++)
        for(i=0; i<N; i++)
            for(v=0; v<n; v++)
                if(W[i][j])
                   fitness += (attributes[i][v] - centers[j][v])*(attributes[i][v] - centers[j][v]); //dal je ona norma iz formule samo apsolutna vrijednost na kvadrat, tj. ti brojevi na kvadrat??
    return fitness;                                                                                 //mislim da da...
}








const char* getfield(char* line, int num)
{
	const char* tok;
	for (tok = strtok(line, ",");
			tok && *tok;
			tok = strtok(NULL, ",\n"))
	{
		if (!--num)
			return tok;
	}
	return NULL;
}

int main(){
    //srand((unsigned)time(NULL));
    int N=8, K=3, n=4; //Mislim da u clanku nije definiran R
    int i, j, k, brojac;

    int *W[N]; //allocate weight matrix s N elemenata i K clustera
    for (i=0; i<N; i++)
         W[i] = (int *)malloc(K * sizeof(int));

    double *centers[K]; //allocate centers of clusters matrix s K clustera and n attributes
    for (i=0; i<K; i++)
         centers[i] = (double *)malloc(n * sizeof(double));

    double fitness;

    double *attributes[N]; //allocate attribute matrix with N elements and n attributes
    for (i=0; i<N; i++)
         attributes[i] = (double *)malloc(n * sizeof(double));

    FILE* stream = fopen("article.csv", "r");
	char line[1024];
	i=0;
	fgets(line, 1024, stream);
    while (fgets(line, 1024, stream))
	{
		char* tmp = strdup(line);
		for(j=0; j<n; j++)
        {
            attributes[i][j]=atof(getfield(tmp, j+1));
            tmp = strdup(line);
		}
        i++;
        free(tmp);
	}
    //printMatrixDouble(attributes, N, n);
    //printf("%f %f\n" , attributes[0][0] , attributes[214][4]);
    //1. iteracija
    //printf("1. iteracija:\n");
    int num;
    char c;
    FILE *in = fopen("output_article.txt" , "r");
    for (k = 0; k < N; k++) //initialize weight matrix for every agent
        for (j = 0; j < K; j++)
            W[k][j]= 0; //kasnije mijenjamo u 1 one koje treba
    for(i=0 ; i<N ; i++)
    {
        fscanf(in, "%d " , &num);
        W[i][num-1]=1;
    }
    centersOfClusters(W, attributes, centers, N, K, n);
    fitness = calculateFitness(W, attributes, centers, N, K, n);
    printf("Fitness nakon metode k-means: %f " , fitness);


    for(i=0 ; i<N ; i++)
    {
        free(attributes[i]);
        free(W[i]);
    }
    for(i=0 ; i<K ; i++)
        free(centers[i]);
    return 0;
}
