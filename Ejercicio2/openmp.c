#include <stdio.h>
#include <stdlib.h>  
#include <sys/time.h>
#include <float.h>
#include <limits.h>
#include<pthread.h>
#include <omp.h>
double *A, *B, *C, *D; //Matrices iniciadas con valores aleatorios
double *P, *R;         //Matrices que utilizan ecuaciones para inicializarse
double *E, *F, *G, *H; //Matrices temporales E=AB, F=EC, G=DC, H=GB
double maxD = DBL_MIN; //MaxD
double minA = DBL_MAX; //MinA
double promP = 0;
int N;
int threads;
int bs;
double suma = 0;


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

void calculoTotal(){
	double *ablk, *bblk, *cblk, *dblk, *eblk, *fblk, *gblk, *hblk, *pblk, *rblk;
	//double minlocal=DBL_MAX, maxlocal= DBL_MIN;
	int I, J, K;    
	int i, j, k;
	//1er for
	
	#pragma omp for reduction(max:maxD) reduction(min:minA)            //Solo es for, porq en este punto llegan los 4 hilos, el for indica que divide el indice del for q viene en la cant de hilos disp.
	for(I = 0; I < N; I += bs)
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
						if (ablk[i*N + j] < minA)
						{
							minA = ablk[i*N + j];
						}
						//MaxD
						if (dblk[i*N + j] > maxD)
						{
							maxD = dblk[i*N + j];
						}
						
					}
				}
			}
		}
	}
	
	//2do for
	#pragma omp for
	for(I = 0; I < N; I += bs)
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
	#pragma omp for reduction (+:suma)
	for(I = 0; I < N; I ++)
	{
		for(J = 0; J < N; J ++)
		{
			//P= MaxD(ABC) + MinA(DCB)
			P[I*N + J] = (F[I*N + J] * maxD) + (H[I*N + J] * minA);
			//#pragma omp critical
			suma += P[I*N + J];
		}
	}
	//#pragma omp barrier //Espero a que todos esten listos ya que debo usar las matrices INTERMEDIAS
	//#pragma omp critical
	promP = (suma/(N*N));
	//4to for
	#pragma omp for
	for(I = 0; I < N; I ++)
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
	if ( (argc != 4) || ((N = atoi(argv[1])) < 1024) || ((bs = atoi(argv[2])) <= 0) || ((N % bs) != 0) || (((atoi(argv[3])) != 4) && ((atoi(argv[3])) != 8))){
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
	
	timetick = dwalltime();
	// Creo los hilos con los que voy a trabajar.
	#pragma omp parallel num_threads(threads)
		calculoTotal();
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
 
	return 0;
}
