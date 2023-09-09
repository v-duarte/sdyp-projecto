#include <stdio.h>
#include <stdlib.h>  
#include <sys/time.h>
#include <float.h>
#include <limits.h>
#include<pthread.h>
#include<math.h>
double *A, *B, *C, *D;  //Matrices iniciadas con valores aleatorios
double *P, *R; 			//Matrices que utilizan ecuaciones para inicializarse
double *E, *F, *G, *H;  //Matrices temporales E=AB, F=EC, G=DC, H=GB, X=MaxD(F), Y=MinA(H)
double *maxD;    		//Vector para max
double *minA;  			//vector para min 
double *suma;			//vector para suma
double promP = 0;
int N;
int threads;
int bs;
//Definicion de barreras
pthread_barrier_t barrera;        // Para todos los Hilos.
pthread_barrier_t *barreras;	  // Para los que van quedando del reduce.
pthread_barrier_t *barreras_dos;  // Para los pares de posiciones a reducir.

double dwalltime(){
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

void* reduce_valores(int id, double minlocal, double maxlocal){		// Id del hilo
    int indice = 0;					//Para manejo de vueltas en reduccion
    int mid_T = threads / 2;				//Para dividir el vector a la mitad.
    int barrier_id;					//Para identificar posicion para barrera de pares
	
	minA[id] = minlocal;
	maxD[id] = maxlocal;
	
    printf("Hilo id:%d\n",id);		//imprime id del hilo.
    while (mid_T > 0){
        barrier_id = id % mid_T;							//El hilo obtiene la posicion para la barrera de par.
        pthread_barrier_wait(&barreras_dos[barrier_id]);	//El hilo espera que su par halla depositado su valor.
        if (id < mid_T){									//Solo ingresa el hilo que se encuentra de la mitad inferior del vector
            if (minA[(id + mid_T)] < minA[id]){				//Evalua cual de los 2 valores es menor, y se queda con el menor.		
                minA[id] = minA[(id + mid_T)];
            }
			if (maxD[(id + mid_T)] > maxD[id]){				//Evalua cual de los 2 valores es mayor, y se queda con el mayor.		
                maxD[id] = maxD[(id + mid_T)];
            }
            mid_T = mid_T / 2;								//Reduzco el vector a la mitad.
        }
        else {
            break;											//Los hilos que se encuentra de la mitad superior del vector son liberados
        }
        pthread_barrier_wait(&barreras[indice]);			//Barrera para que todos los hilos que se encuentren en la mitad inferior esperen que todos terminen
        indice++;											//Para que la proxima vuelta esperen la mitad.
    }
}

void* reduce_suma(int id, double sumalocal ){		// Id del hilo
    int indice = 0;					//Para manejo de vueltas en reduccion
    int mid_T = threads / 2;				//Para dividir el vector a la mitad.
    int barrier_id;					//Para identificar posicion para barrera de pares
	
	suma[id]= sumalocal; 			//Guarda su dato en su posicion.
	
    printf("Hilo id:%d\n",id);		//imprime id del hilo.
    while (mid_T > 0){				
        barrier_id = id % mid_T;							//El hilo obtiene la posicion para la barrera de par.
        pthread_barrier_wait(&barreras_dos[barrier_id]);	//El hilo espera que su par halla depositado su valor.
        if (id < mid_T){									//Solo ingresa el hilo que se encuentra de la mitad inferior del vector		
            suma[id] += suma[id + mid_T];						//Suma su valor, mas el del par y lo guarda en su posicion
            mid_T = mid_T / 2;								//Reduzco el vector a la mitad.
        }
        else {
            break;											//Los hilos que se encuentra de la mitad superior del vector son liberados
        }
        pthread_barrier_wait(&barreras[indice]);			//Barrera para que todos los hilos que se encuentren en la mitad inferior esperen que todos terminen
        indice++;											//Para que la proxima vuelta esperen la mitad.
    }
	if (id == 0 ){
		promP = (suma[0]/(N*N));							//Ya calculamos el promedio
	}
}

//Realiza todos los calculos, 4 Secciones.
void* calculoTotal(void *arg){
	double *ablk, *bblk, *cblk, *dblk, *eblk, *fblk, *gblk, *hblk, *pblk, *rblk;
	double minlocal=DBL_MAX, maxlocal= DBL_MIN;
	int I, J, K;    
	int i, j, k;
	
	int id = *(int*) arg;
	int inicio = id*(N/threads);
	int fin = inicio +((N/threads));
	
	//1er seccion obtiene AB DC minA maxD
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
	
	// Se llama a reducir el maximo y minimo dejandolos en maxD y minA.
	reduce_valores(id, minlocal, maxlocal);
	//2da seccion obtiene ABC DCB
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
	
	// MaxD y MinA debe estar listo     No los espere antes, dado que para la 2da seccion no son requeridos.
	pthread_barrier_wait (&barrera);
	
	//3era seccion obtiene (Calculo P con su promedio)
	double sumalocal = 0;
	for(I = inicio; I < fin; I ++)
	{
		for(J = 0; J < N; J ++)
		{
			//P= MaxD[0](ABC) + MinA[0](DCB)
			P[I*N + J] = (F[I*N + J] * maxD[0]) + (H[I*N + J] * minA[0]);
			sumalocal += P[I*N + J];
		}	
	}
	reduce_suma(id, sumalocal);     // reduce la suma y calcula promP
	pthread_barrier_wait(&barrera); //Espero a que todos esten listos ya que debo usar las matrices INTERMEDIAS
	
    //4ta seccion obtiene R = [promedio (p)]*p
	for(I = inicio; I < fin; I ++)
	{
		for(J = 0; J < N; J ++)
		{
			//R= PromP(P)
			R[I*N + J] = (P[I*N + J] * promP) ;
		}
    }
	pthread_exit(NULL);
 }
 
  
int main(int argc, char *argv[]){
	int i, j;
	char exito=1;

	double timetick;
	
	// Chequeo de parametros
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
	//Aloco memoria para vectores del reduce
	maxD = (double*)malloc(sizeof(double)*threads);
	minA = (double*)malloc(sizeof(double)*threads);
	suma = (double*)malloc(sizeof(double)*threads);

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
	//Aloco memoria para las barreras
	barreras_dos = (pthread_barrier_t*)malloc(sizeof(pthread_barrier_t)*(threads/2));
	barreras = (pthread_barrier_t*)malloc(sizeof(pthread_barrier_t)*((int)log2(threads)));
	printf("Realizando operaciones con matrices de %d x %d en bloques de %d x %d con %d hilos\n", N, N, bs, bs, threads);
	int threads_id[threads];
	pthread_t hilos[threads];
	
	//Inicializacion de barreras
	pthread_barrier_init(&barrera, NULL, threads);
	for(int i=0;i<threads/2;i++){
        pthread_barrier_init(&barreras_dos[i], NULL, 2);   
    }
	for (i= (int)log2(threads) ; i< 0; i--){
		pthread_barrier_init(&barreras[(int)(log2(threads)) -i], NULL, (int)pow(2, (i-1)));
	}
	
	//empiezo a contar el tiempo
	timetick = dwalltime();
	
	//Creacion de hilos
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
	free(barreras_dos);
	free(barreras);
	pthread_barrier_destroy(&barrera);
	pthread_barrier_destroy(barreras_dos);
	pthread_barrier_destroy(barreras);
 
	return 0;
}