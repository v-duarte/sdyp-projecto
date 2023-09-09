#include <stdio.h>
#include <stdlib.h>  
#include <sys/time.h>
#include <float.h>
#include <limits.h>
#include<pthread.h>
#include<math.h>
double *A, *B, *C, *D; //Matrices iniciadas con valores aleatorios
double *P, *R; //Matrices que utilizan ecuaciones para inicializarse
double *E, *F, *G, *H; //Matrices temporales E=AB, F=EC, G=DC, H=GB, X=MaxD(F), Y=MinA(H)
double maxD = DBL_MIN; //MaxD
double minA = DBL_MAX; //MinA
double promP = 0;
int N;
int threads;
int bs;
double suma = 0;
pthread_mutex_t maxMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t minMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t promMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t sumaMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_barrier_t barrera;
pthread_barrier_t *barreras;

double dwalltime(){
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

/* genera un valor tipo double entre el rango min y max*/
double randouble(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

void* calculoTotal(void *arg){
	double *ablk, *bblk, *cblk, *dblk, *eblk, *fblk, *gblk, *hblk, *pblk, *rblk;
	double minlocal=DBL_MAX, maxlocal= DBL_MIN;
	int I, J, K;    
	int i, j, k;
	
	int id = *(int*) arg;
	int inicio = id*(N/threads);
	int fin = inicio +((N/threads));
	//1er for
	for(I = inicio; I < fin; I += bs)
	{
		for(J = 0; J < N; J += bs)
		{
			eblk = &E[I*N + J];
			gblk = &G[I*N + J];
			for(K = 0; K < N; K += bs)
			{
				ablk = &A[I*N + K];
				bblk = &B[J*N + K];
				cblk = &C[J*N + K];
				dblk = &D[I*N + K];		
				for (i = 0; i < bs; i++)
				{
					for (j = 0; j < bs; j++)
					{
						for  (k = 0; k < bs; k++)
						{
							//AB
							eblk[i*N + j] += ablk[i*N + k] * bblk[j*N + k];
							//DC
							gblk[i*N + j] += dblk[i*N + k] * cblk[j*N + k];
						}
						//MinA
						if (ablk[i*N + j] < minlocal)
						{
							minlocal = ablk[i*N + j];
						}
						//MaxD
						if (dblk[i*N + j] > maxlocal)
						{
							maxlocal = dblk[i*N + j];
						}
						
					}
				}
			}
		}
	}
	//Accedo a la variable compartida minimo
	pthread_mutex_lock(&minMutex);
	if ( minA > minlocal )
		minA = minlocal;
	pthread_mutex_unlock(&minMutex);
	//Accedo a la variable compartida maximo
	pthread_mutex_lock(&maxMutex);
	if ( maxD < maxlocal )
		maxD = maxlocal;
	pthread_mutex_unlock(&maxMutex);
	//2do for
	for(I = inicio; I < fin; I += bs)
	{
		for(J = 0; J < N; J += bs)
		{
			fblk = &F[I*N + J];
			hblk = &H[I*N + J];
			for(K = 0; K < N; K += bs)
			{
			eblk = &E[I*N + K];
			gblk = &G[I*N + K];
			cblk = &C[J*N + K];
			bblk = &B[J*N + K];
			for (i = 0; i < bs; i++)
				{
					for (j = 0; j < bs; j++)
					{
						for  (k = 0; k < bs; k++)
						{
							//ABC
							fblk[i*N + j] += eblk[i*N + k] * cblk[j*N + k];
							//DCB
							hblk[i*N + j] += gblk[i*N + k] * bblk[j*N + k];
						}
					}
				}
			}
		}
	}
	//3er for (Calculo P con su promedio)
	for(I = inicio; I < fin; I ++)
	{
		double sumalocal = 0;
		for(J = 0; J < N; J ++)
		{
			//P= MaxD(ABC) + MinA(DCB)
			P[I*N + J] = (F[I*N + J] * maxD) + (H[I*N + J] * minA);
			sumalocal += P[I*N + J];
		}
		pthread_mutex_lock(&sumaMutex); //accedo a la variable compartida suma
		suma += sumalocal;
		pthread_mutex_unlock(&sumaMutex);
	}
	pthread_mutex_lock(&promMutex); //calculo el promedio
	promP = (suma/(N*N));
	pthread_mutex_unlock(&promMutex);
	pthread_barrier_wait(&barrera); //Espero a que todos esten listos ya que debo usar las matrices INTERMEDIAS
	//4to for
	for(I = inicio; I < fin; I ++)
	{
		for(J = 0; J < N; J ++)
		{
			//R= PromP(P)
			R[I*N + J] = (P[I*N + J] * promP) ;
		}
    }
  }
  
int main(int argc, char *argv[]){
  int i, j;
  char exito=1;

  double timetick;

  // Chequeo de parámetros
	if ( (argc != 4) || ((N = atoi(argv[1])) < 1024) || ((bs = atoi(argv[2])) <= 0) || ((N % bs) != 0) || (((atoi(argv[3])) != 4) && ((atoi(argv[3])) != 8)) ){
		printf("Error en los parámetros. Usar: ./%s N (minimo 1024), BS (N debe ser multiplo de BS), T hilos (usar 4 u 8)\n", argv[0]);
		exit(1);
	}
	threads = atoi(argv[3]);
  // Alocar memoria
	A = (double *) malloc(N*N*sizeof(double));
	B = (double *) malloc(N*N*sizeof(double));
	C = (double *) malloc(N*N*sizeof(double));
	D = (double *) malloc(N*N*sizeof(double));
	E = (double *) malloc(N*N*sizeof(double));
	F = (double *) malloc(N*N*sizeof(double));
	G = (double *) malloc(N*N*sizeof(double));
	H = (double *) malloc(N*N*sizeof(double));
	P = (double *) malloc(N*N*sizeof(double));
	R = (double *) malloc(N*N*sizeof(double));

  // Inicializacion
	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
			A[i*N + j] = 1.0; //Inicializo por fila
			B[j*N + i] = 2.0; //Inicializo por columna
			C[j*N + i] = 3.0; //Inicializo por columna
			D[i*N + j] = 4.0; //Inicializo por fila
			
			E[i*N + j] = 0.0;
			F[i*N + j] = 0.0;
			G[i*N + j] = 0.0;
			H[i*N + j] = 0.0;
			P[i*N + j] = 0.0;
			R[i*N + j] = 0.0;
		}
	}


	printf("Realizando operaciones con matrices de %d x %d en bloques de %d x %d con %d hilos\n", N, N, bs, bs, threads);
	int threads_id[threads];
	pthread_t hilos[threads];
  pthread_barrier_init(&barrera, NULL, threads);
  //Inicializar barreras
  for (i= log2(threads) ; i< 0; i--){
  	pthread_barrier_init(&barreras[log2(threads) -i], NULL, pow(2, (i-1)));
  }
  timetick = dwalltime();
  for(i=0; i<threads; i++){
	 	threads_id[i]=i;
	 	pthread_create(&hilos[i], NULL, &calculoTotal, (void*)&threads_id[i]);
	}
	
  for (i=0; i<threads; i++){
		pthread_join(hilos[i],NULL);
	}
double totalTime = dwalltime() - timetick;

	printf("Tiempo en segundos %f\n",totalTime);
  // Validando
  printf("Validando P y R...\n");
  //La validacion asume que la matriz esta cargada con el mismo numero
  double validarP = N*N*A[0]*B[0]*C[0]*D[0]*2;
  double validarR = validarP*validarP;
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      if (P[i*N + j] != validarP)
      {
        printf("Error en P[%d, %d], valor: %f\n", i, j, P[i*N + j]);
		if(exito)
				exito=0;
      }
      if (R[i*N + j] != validarR)
      {
        printf("Error en R[%d, %d], valor: %f\n", i, j, R[i*N + j]);
		if(exito)
			exito=0;
      }
    }
  }
  if (exito)
	printf("Validacion exitosa.\n");
  else
	  printf("Validacion fallida.\n");
 

	free(A);
	free(B);
	free(C);
	free(D);
	free(E);
	free(F);
	free(G);
	free(H);
	free(P);
	free(R);
	pthread_mutex_destroy(&minMutex);
	pthread_mutex_destroy(&maxMutex);
	pthread_mutex_destroy(&promMutex);
	pthread_mutex_destroy(&sumaMutex);
	pthread_barrier_destroy(&barrera);
 
	return 0;
}

void* reduce(void *arg){
    int tid=*(int*)arg;
    int index_inf = tid * N/T;
    int index_sup = (tid + 1) * N/T;
    int indice = 0;
    //Realiza el conteo
    partial[tid] = 0;
    for(int i=index_inf; i<index_sup; i++){
        if (A[i] < partial[tid]){
            partial[tid] = A[i];
        }
    }
    int mid_T = T / 2;
    int barrier_id;
    printf("Hilo id:%d\n",tid);
    while (mid_T > 0){
        barrier_id = tid % mid_T;
        pthread_barrier_wait(&barreras[barrier_id]);
        if (tid < mid_T){
            if (partial[tid] > partial[tid+mid_T]){
                partial[tid] = partial[tid+mid_T];
            }
            mid_T = mid_T / 2;
        }
        else {
            break;
        }
        pthread_barrier_wait(&barreras[indice]);
        indice++;
    }
    if (tid == 0) {
        ocurrencias = partial[tid];
    } 
    pthread_exit(NULL);
    pthread_barrier_wait(&barrera);
    
}
