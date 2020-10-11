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

int clusterWithMaxPheromoneTrail(double** F, int element, int K){
    int j, maxPheromone = 1;
    double maxTrail=F[element][0];
    for(j=1; j<K; j++)
    {
        if(F[element][j] > maxTrail)
        {
            maxPheromone = j+1;
            maxTrail=F[element][j];
        }
        else if(F[element][j] == maxTrail)
        //printf("Index: %d\n" , j);
            if(rand()%2)
            {
                maxPheromone = j+1;
                maxTrail=F[element][j];
            }
    }
    return maxPheromone;
}

int clusterChosenStochastically(double** F, int element, int K){
    double* normalizedPheromone = (double*) malloc(K * sizeof(double)); // array of normalized pheromone probabilities for given element
    double sumPheromones = 0.0, lastSum;
    int i;
    for(i=0; i<K; i++)
        sumPheromones += F[element][i]; // first sum all pheromone values
    for(i=0; i<K; i++)
        normalizedPheromone[i] = F[element][i]/sumPheromones; //calculate normalized pheromone array for given element
    double randomNumber = (double)rand()/(double)RAND_MAX; //generate random number
    if(randomNumber < normalizedPheromone[0])
        return 1; //if 1st cluster
    else
    {
        sumPheromones = normalizedPheromone[0] + normalizedPheromone[1]; // second sum pheromone values for if condition
        lastSum = normalizedPheromone[0];
    }
    for(i=1; i<K-1; i++)
        if((lastSum <= randomNumber) && (randomNumber < sumPheromones))
            return i+1;
        else
        {
            lastSum = sumPheromones;
            sumPheromones += normalizedPheromone[i+1];
        }
    return K;
}

void generateSolution(int** S, double** F, int agent, int N, int K, int** W){
    int i, j, cluster;
    double q0 = 0.01; //u clanku je 0.98!! ali sa smanjenjem q0 dobijemo puno bolje rezultate, znaci ako ne pratimo toliko feromone... u ovom slucaju bas zglavi s q0=0.98
    double* randomNumbers = (double *) malloc(N * sizeof(double)); // generate N random numbers between 0 and 1
    for(i = 0; i<N; i++)
        randomNumbers[i] = (double)rand()/(double)RAND_MAX;
    /*for(i = 0; i<N; i++)
        printf("%f ", randomNumbers[i]);
    printf("\n");*/

    for(i = 0; i<N; i++)
        if(randomNumbers[i] < q0){
            cluster = clusterWithMaxPheromoneTrail(F, i, K);
            //printf("ClusterWithMaxPheromoneTrail: %d\n" , cluster);
            S[agent][i] = cluster;
            W[i][cluster-1]=1;
        }
        else{
            cluster = clusterChosenStochastically(F, i, K);
            //printf("ClusterChosenStochastically: %d\n" , cluster);
            S[agent][i] = cluster;
            W[i][cluster-1]=1;
        }
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

void findBestL(double* fitness, double percent , int R, int* bestInd)
{
    int L = (int)(percent*R);
    //printf("Postotak: %f\n" , percent);
    //printf("Izracunati L: %d\n" , L);
    double fit[R] , sum = 0.0;
    int i, j, tempInd;
    double temp;
    for(i=0 ; i<R ; i++)
    {
        fit[i]=fitness[i];
        sum+=fit[i];
    }
    for(i=0 ; i<L ; i++)
    {
        temp = sum; //Sigurno ce postojat neki manji od ukupne sume, da stavim fitness[0] ne valja(jer mozda je fit[0]=-1 pa necemo nac drugi minimalni
        for(j=0 ; j<R ; j++)
        {
            if(fit[j]<temp && fit[j]!=-1)
            {
                temp = fit[j];
                tempInd = j;
            }
        }
        fit[tempInd]=-1; //Nacin da vec nadjeni minimalni ne gledamo vise kako bi nasli ostale minimalne
        bestInd[i]=tempInd;
    }
}

void localSearch(int** S, int N, int* best, int L, int K, double* fitness, double** attributes, int n)
{
    double pls = 0.01;
    double randomNumber, tempArray[N];
    int i, j, random;
    double *centers[K];
    for (i=0; i<K; i++)
         centers[i] = (double *)malloc(n * sizeof(double));
    int *WTemp[N];
    for (i=0; i<N; i++)
         WTemp[i] = (int *)malloc(K * sizeof(int));
    for(i=0 ; i<N ; i++)
        for(j=0 ; j<K ; j++)
            WTemp[i][j] = 0;
    for(i=0 ; i<L ; i++)
    {
        //printf("Broj agenta: %d\n", best[i]+1);
        for(j = 0; j<N; j++)
        {
            randomNumber = (double)rand()/(double)RAND_MAX;
            //printf("Random number: %f\n", randomNumber);
            if(randomNumber <= pls)
            {
                random = (rand() % K) + 1; //Random number in range [1,K]
                while(random == S[best[i]][j])
                    random = (rand() % K) + 1; //Mora random bit razlicit od trenutnog clustera
                //printf("Local search cluster number: %f\n", random);
                tempArray[j] = random;
                WTemp[j][random-1] = 1;
            }
            else
            {
                tempArray[j] = S[best[i]][j];
                WTemp[j][(S[best[i]][j]-1)] = 1;
            }
        }
        centersOfClusters(WTemp, attributes, centers, N, K, n);
        double fitness2 = calculateFitness(WTemp, attributes, centers, N, K, n);
        if(fitness2 < fitness[best[i]])
        {
            for(j=0 ; j<N ; j++)
            {
                S[best[i]][j] = tempArray[j];
                fitness[best[i]] = fitness2;
            }
        }
    }
}

void evaporatePheromone(double** F, int N, int K)
{
    int i, j;
    double ro = 0.01;
    for(i=0 ; i<N ; i++)
        for(j=0 ; j<K ; j++)
            F[i][j] = (1-ro)*(F[i][j]);
}

void addPheromone(int* best, int L, double** F, int N, int K, int** S, double* fitness)
{
    int i, j;
    for(i=0 ; i<L ; i++)
    {
        for(j=0 ; j<N ; j++)
        {
            F[j][S[best[i]][j]-1] += 1/(fitness[best[i]]);
        }
    }
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
    int R=50, N=150, K=3, n=4, maxIterations=1300; //Mislim da u clanku nije definiran R
    int i, j, k, brojac;
    double percentForLocalSearch = 0.4, percentForPheromonUpdate = 0.03;
    int L_ls = (int)(R*percentForLocalSearch);
    int L_pu = (int)(R*percentForPheromonUpdate);
    int* best_ls = (int *)malloc(L_ls*sizeof(int));
    int* best_pu = (int *)malloc(L_pu*sizeof(int));

    int *S[R]; //allocate solution matrix s R agenata i N elemenata
    for (i=0; i<R; i++)
         S[i] = (int *)malloc(N * sizeof(int));

    double *F[N]; //allocate pheromone matrix s N elemenata i K clustera
    for (i=0; i<N; i++)
         F[i] = (double *)malloc(K * sizeof(double));
    for (i = 0; i <  N; i++) //initialize pheromone matrix
      for (j = 0; j < K; j++)
        F[i][j]= 0.01; //samo nju iniciliziramo izvan for petlje

    int *W[N]; //allocate weight matrix s N elemenata i K clustera
    for (i=0; i<N; i++)
         W[i] = (int *)malloc(K * sizeof(int));

    double *centers[K]; //allocate centers of clusters matrix s K clustera and n attributes
    for (i=0; i<K; i++)
         centers[i] = (double *)malloc(n * sizeof(double));

    double* fitness = (double*)malloc(R*sizeof(double));

    double *attributes[N]; //allocate attribute matrix with N elements and n attributes
    for (i=0; i<N; i++)
         attributes[i] = (double *)malloc(n * sizeof(double));

    FILE* stream = fopen("Iris.csv", "r");
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

    //1. iteracija
    //printf("1. iteracija:\n");
    int cluster;
    for (i = 0; i <  R; i++) //initialize solution matrix
          for (j = 0; j < N; j++)
            S[i][j]=0;
        for(i=0; i<R; i++)
        {
            for (k = 0; k < N; k++) //initialize weight matrix for every agent
                for (j = 0; j < K; j++)
                    W[k][j]= 0; //kasnije mijenjamo u 1 one koje treba
            for(j = 0; j<N; j++) // u 1. iteraciji clustere dodjeljujemo nasumicno
            {
                cluster=(rand() % K) + 1; //Random number in range [1,K]
                S[i][j] = cluster;
                W[j][cluster-1]=1;
            }
            centersOfClusters(W, attributes, centers, N, K, n);
            fitness[i] = calculateFitness(W, attributes, centers, N, K, n);
            //printf("Fitness, %d. agent: %f\n" , i+1 , fitness[i]);
        }
        findBestL(fitness, percentForLocalSearch , R, best_ls);
        //for(i=0 ; i<L ; i++) printf("Indeks: %d \n" , best[i]+1);
        //printMatrixInt(S,R,N);
        localSearch(S, N, best_ls, L_ls, K, fitness, attributes, n);
        //for(i=0 ; i<L ; i++) printf("Indeks: %d \n" , best[i]+1);
        //printMatrixInt(S,R,N);
        /*for(i=0 ; i<L ; i++){
            //printf("Indeks: %d \n" , best[i]+1);
            printf("Fitness nakon local search: %f\n", fitness[best[i]]);
        }
        printf("\n");*/
        printf("Najbolji fitness 1. iteracije: %f\n", fitness[best_ls[0]]);
        //printf("Feromoni prije updatea: \n");
        //printMatrixDouble(F,N,K);
        findBestL(fitness, percentForPheromonUpdate , R, best_pu);
        evaporatePheromone(F, N, K);
        addPheromone(best_pu, L_pu, F, N, K, S, fitness);
        //printf("Feromoni poslije updatea: \n");
        //printMatrixDouble(F,N,K);


    //ostale iteracije
    for(brojac=1; brojac<maxIterations; brojac++)
    {
        //printf("%d. iteracija:\n", brojac+1);
        for (i = 0; i <  R; i++) //initialize solution matrix
          for (j = 0; j < N; j++)
            S[i][j]=0;
        for(i=0; i<R; i++)
        {
            for (k = 0; k < N; k++) //initialize weight matrix for every agent
                for (j = 0; j < K; j++)
                    W[k][j]= 0; //kasnije mijenjamo u 1 one koje treba
            generateSolution(S, F, i, N, K, W);
            centersOfClusters(W, attributes, centers, N, K, n);
            fitness[i] = calculateFitness(W, attributes, centers, N, K, n);
            //printf("Fitness, %d. agent: %f\n" , i+1 , fitness[i]);
        }
        findBestL(fitness, percentForLocalSearch , R, best_ls);
        //for(i=0 ; i<L ; i++) printf("Indeks: %d \n" , best[i]+1);
        //printMatrixInt(S,R,N);
        localSearch(S, N, best_ls, L_ls, K, fitness, attributes, n);
        /*for(i=0 ; i<L ; i++) printf("Indeks: %d \n" , best[i]+1);
        printMatrixInt(S,R,N);
        for(i=0 ; i<L ; i++){
            //printf("Indeks: %d \n" , best[i]+1);
            printf("Fitness nakon local search: %f\n", fitness[best[i]]);
        }
        printf("\n");
        printf("Feromoni prije updatea: \n");
        printMatrixDouble(F,N,K);*/
        printf("Najbolji fitness %d. iteracije: %f\n", brojac+1, fitness[best_ls[0]]);
        findBestL(fitness, percentForPheromonUpdate , R, best_pu);
        evaporatePheromone(F, N, K);
        addPheromone(best_pu, L_pu, F, N, K, S, fitness);
        //printf("Feromoni poslije updatea: \n");
        //printMatrixDouble(F,N,K);
    }
    //printf(" Najbolji fitness zadnje iteracije: %f\n", fitness[best[0]]);
    //printMatrixInt(S,R,N);
    //printMatrixDouble(F,N,K);

    for(i=0 ; i<N ; i++)
        printf("%d " , S[best_ls[0]][i]);
    printf("\n");
    printf("%d " , best_ls[0]);

    for(i=0 ; i<N ; i++)
    {
        free(attributes[i]);
        free(F[i]);
        free(W[i]);
    }
    for(i=0 ; i<R ; i++)
        free(S[i]);
    for(i=0 ; i<K ; i++)
        free(centers[i]);
    free(best_ls);
    free(best_pu);
    return 0;
}

/*double attributes[150][4]={{5.1,3.5,1.4},{4.9,3,1.4,.2},{4.7,3.2,1.3,.2},{4.6,3.1,1.5,.2}, {5,3.6,1.4,.2}, {5.4,3.9,1.7,.4}, {4.6,3.4,1.4,.3}, {5,3.4,1.5,.2}, {4.4,2.9,1.4,.2}, {4.9,3.1,1.5,.1},
        {5.4,3.7,1.5,.2}, {4.8,3.4,1.6,.2}, {4.8,3,1.4,.1}, {4.3,3,1.1,.1}, {5.8,4,1.2,.2}, {5.7,4.4,1.5,.4}, {5.4,3.9,1.3,.4}, {5.1,3.5,1.4,.3}, {5.7,3.8,1.7,.3}, {5.1,3.8,1.5,.3}, {5.4,3.4,1.7,.2},
        {5.1,3.7,1.5,.4}, {4.6,3.6,1,.2}, {5.1,3.3,1.7,.5}, {4.8,3.4,1.9,.2}, {5,3,1.6,.2}, {5,3.4,1.6,.4}, {5.2,3.5,1.5,.2}, {5.2,3.4,1.4,.2}, {4.7,3.2,1.6,.2}, {4.8,3.1,1.6,.2}, {5.4,3.4,1.5,.4},
        {5.2,4.1,1.5,.1}, {5.5,4.2,1.4,.2}, {4.9,3.1,1.5,.2}, {5,3.2,1.2,.2}, {5.5,3.5,1.3,.2}, {4.9,3.6,1.4,.1}, {4.4,3,1.3,.2}, {5.1,3.4,1.5,.2}, {5,3.5,1.3,.3}, {5,3.5,1.3,.3}, {4.4,3.2,1.3,.2},
        {5,3.5,1.6,.6}, {5.1,3.8,1.9,.4}, {4.8,3,1.4,.3}, {5.1,3.8,1.6,.2}, {4.6,3.2,1.4,.2}, {5.3,3.7,1.5,.2}, {5,3.3,1.4,.2}, {7,3.2,4.7,1.4}, {6.4,3.2,4.5,1.5}, {6.9,3.1,4.9,1.5}, {5.5,2.3,4,1.3},
        {6.5,2.8,4.6,1.5}, {5.7,2.8,4.5,1.3}, {6.3,3.3,4.7,1.6}, {4.9,2.4,3.3,1}, {6.6,2.9,4.6,1.3}, {5.2,2.7,3.9,1.4},{5,2,3.5,1},{5.9,3,4.2,1.5},{6,2.2,4,1},{6.1,2.9,4.7,1.4},{5.6,2.9,3.6,1.3},
        {6.7,3.1,4.4,1.4},{5.6,3,4.5,1.5},{5.8,2.7,4.1,1},{6.2,2.2,4.5,1.5},{5.6,2.5,3.9,1.1},{5.9,3.2,4.8,1.8},{6.1,2.8,4,1.3},{6.3,2.5,4.9,1.5},{6.1,2.8,4.7,1.2},{6.4,2.9,4.3,1.3},{6.6,3,4.4,1.4},
        {6.8,2.8,4.8,1.4},{6.7,3,5,1.7},{6,2.9,4.5,1.5},{5.7,2.6,3.5,1},{5.5,2.4,3.8,1.1},{5.5,2.4,3.7,1},{5.8,2.7,3.9,1.2},{6,2.7,5.1,1.6},{5.4,3,4.5,1.5},{6,3.4,4.5,1.6},{6.7,3.1,4.7,1.5},{6.3,2.3,4.4,1.3},
        {5.6,3,4.1,1.3},{5.5,2.5,4,1.3},{5.5,2.6,4.4,1.2},{6.1,3,4.6,1.4},{5.8,2.6,4,1.2},{5,2.3,3.3,1},{5.6,2.7,4.2,1.3},{5.7,3,4.2,1.2},{5.7,2.9,4.2,1.3},{6.2,2.9,4.3,1.3},{5.1,2.5,3,1.1},{5.7,2.8,4.1,1.3},
        {6.3,3.3,6,2.5},{5.8,2.7,5.1,1.9},{7.1,3,5.9,2.1},{6.3,2.9,5.6,1.8},{6.5,3,5.8,2.2},{7.6,3,6.6,2.1},{4.9,2.5,4.5,1.7},{7.3,2.9,6.3,1.8},{6.7,2.5,5.8,1.8},{7.2,3.6,6.1,2.5},{6.5,3.2,5.1,2},{6.4,2.7,5.3,1.9},
        {6.8,3,5.5,2.1},{5.7,2.5,5,2},{5.8,2.8,5.1,2.4},{6.4,3.2,5.3,2.3},{6.5,3,5.5,1.8},{7.7,3.8,6.7,2.2},{7.7,2.6,6.9,2.3},{6,2.2,5,1.5},{6.9,3.2,5.7,2.3},{5.6,2.8,4.9,2},{7.7,2.8,6.7,2},{6.3,2.7,4.9,1.8},
        {6.7,3.3,5.7,2.1},{7.2,3.2,6,1.8},{6.2,2.8,4.8,1.8},{6.1,3,4.9,1.8},{6.4,2.8,5.6,2.1},{7.2,3,5.8,1.6},{7.4,2.8,6.1,1.9},{7.9,3.8,6.4,2},{6.4,2.8,5.6,2.2},{6.3,2.8,5.1,1.5},{6.1,2.6,5.6,1.4},{7.7,3,6.1,2.3},
        {6.3,3.4,5.6,2.4},{6.4,3.1,5.5,1.8},{6,3,4.8,1.8},{6.9,3.1,5.4,2.1},{6.7,3.1,5.6,2.4},{6.9,3.1,5.1,2.3},{5.8,2.7,5.1,1.9},{6.8,3.2,5.9,2.3},{6.7,3.3,5.7,2.5},{6.7,3,5.2,2.3},{6.3,2.5,5,1.9},{6.5,3,5.2,2},
        {6.2,3.4,5.4,2.3},{5.9,3,5.1,1.8}};*/

