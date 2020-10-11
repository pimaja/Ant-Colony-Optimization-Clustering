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
        if(F[element][j] > maxTrail) //Ako ih je vise jednakih uzima se za najveceg prvi kojeg nadjemo
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
    double q0 = 0.001; //u clanku je 0.98!! ali sa smanjenjem q0 dobijemo puno bolje rezultate, znaci ako ne pratimo toliko feromone...
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
    int R=50, N=215, K=3, n=5, maxIterations=4600; //Mislim da u clanku nije definiran R, za R=10 se dobiju najbolji rezultati, isprobavala sam do R=100
    int i, j, k, brojac;
    double percentForLocalSearch = 0.7, percentForPheromonUpdate = 0.03;
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

    FILE* stream = fopen("Thyroid.csv", "r");
	char line[1024];
	i=0;
    while (fgets(line, 1024, stream))
	{
		char* tmp = strdup(line);
		for(j=0; j<n; j++)
        {
            attributes[i][j]=atof(getfield(tmp, j+2));
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
                W[j][cluster]=1;
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
        printf("Fitness: %f\n", fitness[best_ls[0]]);
        //printf("Feromoni prije updatea: \n");
        //printMatrixDouble(F,N,K);
        findBestL(fitness, percentForPheromonUpdate , R, best_pu);
        evaporatePheromone(F, N, K);
        addPheromone(best_pu, L_pu, F, N, K, S, fitness);
        //printf("Feromoni poslije updatea: \n");
        //printMatrixDouble(F,N,K);


    //ostale iteracije
    for(brojac=0; brojac<maxIterations; brojac++)
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
        //printf("Fitness: %f\n", fitness[best[0]]);
        findBestL(fitness, percentForPheromonUpdate , R, best_pu);
        evaporatePheromone(F, N, K);
        addPheromone(best_pu, L_pu, F, N, K, S, fitness);
        //printf("Feromoni poslije updatea: \n");
        //printMatrixDouble(F,N,K);
        printf("Najbolji fitness %d. iteracije: %f\n", brojac+1, fitness[best_pu[0]]);
    }
    //printf("Fitness: %f\n", fitness[best[0]]);
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

/*double attributes[215][5]={{1,107,10.1,2.2,0.9,2.7},{113,9.9,3.1,2.0,5.9},{127,12.9,2.4,1.4,0.6},{109,5.3,1.6,1.4,1.5},{105,7.3,1.5,1.5,-0.1},{105,6.1,2.1,1.4,7.0},{110,10.4,1.6,1.6,2.7},{114,9.9,2.4,1.5,5.7},{106,9.4,2.2,1.5,0.0},{107,13.0,1.1,0.9,3.1},{106,4.2,1.2,1.6,1.4},{110,11.3,2.3,0.9,3.3},{116,9.2,2.7,1.0,4.2},{112,8.1,1.9,3.7,2.0},{122,9.7,1.6,0.9,2.2},{109,8.4,2.1,1.1,3.6},{111,8.4,1.5,0.8,1.2},{114,6.7,1.5,1.0,3.5},{119,10.6,2.1,1.3,1.1},{115,7.1,1.3,1.3,2.0},{101,7.8,1.2,1.0,1.7},{103,10.1,1.3,0.7,0.1},{109,10.4,1.9,0.4,-0.1},{102,7.6,1.8,2.0,2.5},{121,10.1,1.7,1.3,0.1},{100,6.1,2.4,1.8,3.8},{106,9.6,2.4,1.0,1.3},{116,10.1,2.2,1.6,0.8},{105,11.1,2.0,1.0,1.0}
,{110,10.4,1.8,1.0,2.3},{120,8.4,1.1,1.4,1.4},{116,11.1,2.0,1.2,2.3},{110,7.8,1.9,2.1,6.4},{90,8.1,1.6,1.4,1.1},{117,12.2,1.9,1.2,3.9},{117,11.0,1.4,1.5,2.1},{113,9.0,2.0,1.8,1.6},{106,9.4,1.5,0.8,0.5},{130,9.5,1.7,0.4,3.2},{100,10.5,2.4,0.9,1.9},{121,10.1,2.4,0.8,3.0},{110,9.2,1.6,1.5,0.3},{129,11.9,2.7,1.2,3.5},{121,13.5,1.5,1.6,0.5},{123,8.1,2.3,1.0,5.1},{107,8.4,1.8,1.5,0.8},{109,10.0,1.3,1.8,4.3},{120,6.8,1.9,1.3,1.9},{100,9.5,2.5,1.3,-0.2},{118,8.1,1.9,1.5,13.7},{100,11.3,2.5,0.7,-0.3},{103,12.2,1.2,1.3,2.7},{115,8.1,1.7,0.6,2.2},{119,8.0,2.0,0.6,3.2},{106,9.4,1.7,0.9,3.1},{114,10.9,2.1,0.3,1.4},{93,8.9,1.5,0.8,2.7},{120,10.4,2.1,1.1,1.8},{106,11.3,1.8,0.9,1.0},{110,8.7,1.9,1.6,4.4},{103,8.1,1.4,0.5,3.8},{101,7.1,2.2,0.8,2.2},{115,10.4,1.8,1.6,2.0},{116,10.0,1.7,1.5,4.3},{117,9.2,1.9,1.5,6.8},{106,6.7,1.5,1.2,3.9},{118,10.5,2.1,0.7,3.5},{97,7.8,1.3,1.2,0.9},{113,11.1,1.7,0.8,2.3},{104,6.3,2.0,1.2,4.0},{96,9.4,1.5,1.0,3.1},{120,12.4,2.4,0.8,1.9},{133,9.7,2.9,0.8,1.9}
,{126,9.4,2.3,1.0,4.0},{113,8.5,1.8,0.8,0.5},{109,9.7,1.4,1.1,2.1},{119,12.9,1.5,1.3,3.6},{101,7.1,1.6,1.5,1.6},{108,10.4,2.1,1.3,2.4},{117,6.7,2.2,1.8,6.7},{115,15.3,2.3,2.0,2.0},{91,8.0,1.7,2.1,4.6},{103,8.5,1.8,1.9,1.1},{98,9.1,1.4,1.9,-0.3},{111,7.8,2.0,1.8,4.1},{107,13.0,1.5,2.8,1.7},{119,11.4,2.3,2.2,1.6},{122,11.8,2.7,1.7,2.3},{105,8.1,2.0,1.9,-0.5},{109,7.6,1.3,2.2,1.9},{105,9.5,1.8,1.6,3.6},{112,5.9,1.7,2.0,1.3},{112,9.5,2.0,1.2,0.7},{98,8.6,1.6,1.6,6.0},{109,12.4,2.3,1.7,0.8},{114,9.1,2.6,1.5,1.5},{114,11.1,2.4,2.0,-0.3},{110,8.4,1.4,1.0,1.9},{120,7.1,1.2,1.5,4.3},{108,10.9,1.2,1.9,1.0},{108,8.7,1.2,2.2,2.5},{116,11.9,1.8,1.9,1.5}
,{113,11.5,1.5,1.9,2.9},{105,7.0,1.5,2.7,4.3},{114,8.4,1.6,1.6,-0.2},{114,8.1,1.6,1.6,0.5},{105,11.1,1.1,0.8,1.2},{107,13.8,1.5,1.0,1.9},{116,11.5,1.8,1.4,5.4},{102,9.5,1.4,1.1,1.6},{116,16.1,0.9,1.3,1.5},{118,10.6,1.8,1.4,3.0},{109,8.9,1.7,1.0,0.9},{110,7.0,1.0,1.6,4.3},{104,9.6,1.1,1.3,0.8},{105,8.7,1.5,1.1,1.5},{102,8.5,1.2,1.3,1.4},{112,6.8,1.7,1.4,3.3},{111,8.5,1.6,1.1,3.9},{111,8.5,1.6,1.2,7.7},{103,7.3,1.0,0.7,0.5},{98,10.4,1.6,2.3,-0.7},{117,7.8,2.0,1.0,3.9},{111,9.1,1.7,1.2,4.1},{101,6.3,1.5,0.9,2.9},{106,8.9,0.7,1.0,2.3},{102,8.4,1.5,0.8,2.4},{115,10.6,0.8,2.1,4.6},{130,10.0,1.6,0.9,4.6},{101,6.7,1.3,1.0,5.7},{110,6.3,1.0,0.8,1.0},{103,9.5,2.9,1.4,-0.1},{113,7.8,2.0,1.1,3.0},{112,10.6,1.6,0.9,-0.1},{118,6.5,1.2,1.2,1.7},{109,9.2,1.8,1.1,4.4},{116,7.8,1.4,1.1,3.7},{127,7.7,1.8,1.9,6.4},{108,6.5,1.0,0.9,1.5},{108,7.1,1.3,1.6,2.2},{105,5.7,1.0,0.9,0.9},{98,5.7,0.4,1.3,2.8},{112,6.5,1.2,1.2,2.0},{118,12.2,1.5,1.0,2.3}
,{94,7.5,1.2,1.3,4.4},{126,10.4,1.7,1.2,3.5},{114,7.5,1.1,1.6,4.4},{111,11.9,2.3,0.9,3.8},{104,6.1,1.8,0.5,0.8},{102,6.6,1.2,1.4,1.3},{139,16.4,3.8,1.1,-0.2},{111,16.0,2.1,0.9,-0.1},{113,17.2,1.8,1.0,0.0},{65,25.3,5.8,1.3,0.2},{88,24.1,5.5,0.8,0.1},{65,18.2,10.0,1.3,0.1},{134,16.4,4.8,0.6,0.1},{110,20.3,3.7,0.6,0.2},{67,23.3,7.4,1.8,-0.6},{95,11.1,2.7,1.6,-0.3},{89,14.3,4.1,0.5,0.2},{89,23.8,5.4,0.5,0.1},{88,12.9,2.7,0.1,0.2},{105,17.4,1.6,0.3,0.4},{89,20.1,7.3,1.1,-0.2},{99,13.0,3.6,0.7,-0.1},{80,23.0,10.0,0.9,-0.1},{89,21.8,7.1,0.7,-0.1},{99,13.0,3.1,0.5,-0.1},{68,14.7,7.8,0.6,-0.2},{97,14.2,3.6,1.5,0.3},{84,21.5,2.7,1.1,-0.6}
,{84,18.5,4.4,1.1,-0.3},{98,16.7,4.3,1.7,0.2},{94,20.5,1.8,1.4,-0.5},{99,17.5,1.9,1.4,0.3},{76,25.3,4.5,1.2,-0.1},{110,15.2,1.9,0.7,-0.2},{144,22.3,3.3,1.3,0.6},{105,12.0,3.3,1.1,0.0},{88,16.5,4.9,0.8,0.1},{97,15.1,1.8,1.2,-0.2},{106,13.4,3.0,1.1,0.0},{79,19.0,5.5,0.9,0.3},{92,11.1,2.0,0.7,-0.2},{125,2.3,0.9,16.5,9.5},{120,6.8,2.1,10.4,38.6},{108,3.5,0.6,1.7,1.4},{120,3.0,2.5,1.2,4.5},{119,3.8,1.1,23.0,5.7},{141,5.6,1.8,9.2,14.4},{129,1.5,0.6,12.5,2.9},{118,3.6,1.5,11.6,48.8},{120,1.9,0.7,18.5,24.0},{119,0.8,0.7,56.4,21.6},{123,5.6,1.1,13.7,56.3},{115,6.3,1.2,4.7,14.4},{126,0.5,0.2,12.2,8.8},{121,4.7,1.8,11.2,53.0}
,{131,2.7,0.8,9.9,4.7},{134,2.0,0.5,12.2,2.2},{141,2.5,1.3,8.5,7.5},{113,5.1,0.7,5.8,19.6},{136,1.4,0.3,32.6,8.4},{120,3.4,1.8,7.5,21.5},{125,3.7,1.1,8.5,25.9},{123,1.9,0.3,22.8,22.2},{112,2.6,0.7,41.0,19.0},{134,1.9,0.6,18.4,8.2},{119,5.1,1.1,7.0,40.8},{118,6.5,1.3,1.7,11.5},{139,4.2,0.7,4.3,6.3},{103,5.1,1.4,1.2,5.0},{97,4.7,1.1,2.1,12.6},{102,5.3,1.4,1.3,6.7}};*/
