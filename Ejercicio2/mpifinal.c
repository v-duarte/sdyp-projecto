#include <stdio.h>
#include <stdlib.h>  
#include <sys/time.h>
#include <float.h>
#include <limits.h>
#include <mpi.h>

double dwalltime(){
	double sec; 
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

void init(double *a, double *b, double *c, double *d, double *ab, double *dc, double *abc, double *dcb, double *p, double *r, int N, int sizePart){
	 int i, j;
	// Inicializacion
	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
			a[i*N + j] = 1.0; //Inicializo por fila
			b[j*N + i] = 2.0; //Inicializo por columna
			c[j*N + i] = 3.0; //Inicializo por columna
			d[i*N + j] = 4.0; //Inicializo por fila
			
			ab[i*N + j] = 0.0;
			abc[i*N + j] = 0.0;
			dc[i*N + j] = 0.0;
			dcb[i*N + j] = 0.0;
			p[i*N + j] = 0.0;
			r[i*N + j] = 0.0;
		}
	}
 }
//Obtenemos AB CD minA maxD
void mult1(double *a, double *b, double *c, double *d, double *ab, double *dc, double *minlocal, double *maxlocal, int N, int bs, int sizePart){
    double *ablk, *bblk, *cblk, *dblk, *abblk, *dcblk;
	int I, J, K, i, j ,k , min = 999, max= 0;  

    for(I = 0; I < sizePart/N; I += bs)
    {
		for(J = 0; J < N; J += bs)
		{
			//Asignacion de bloques
			abblk = &ab[I*N + J];
			dcblk = &dc[I*N + J];
			
            for(K = 0; K < N; K += bs)
			{
			    ablk = &a[I*N + K];
			    dblk = &d[I*N + K];
			    bblk = &b[J*N + K];
			    cblk = &c[J*N + K];
			
			   for (i = 0; i < bs; i++)
				   {
					for (j = 0; j < bs; j++)
					{
						for  (k = 0; k < bs; k++)
						{
						//AB
						abblk[i*N + j] += ablk[i*N + k] * bblk[j*N + k];
						//DC
						dcblk[i*N + j] += dblk[i*N + k] * cblk[j*N + k];
						}
						//MinA
                        if (ablk [i*N + j] < min)
                        {
                            min = ablk [i*N +j];
                        }
						//MaxD
						if (dblk[i*N + j] > max)
						{
							max = dblk[i*N + j];
                   
						}
					}
				}
			}
		}
	}
    *minlocal = min;
    *maxlocal = max;
}
//Obtenemos ABC y DCB
void mult2(double* ab, double* c, double* dc, double* b, double* abc, double* dcb, int N, int bs, int sizePart){
    double *bblk, *cblk, *abblk, *dcblk, *abcblk, *dcbblk; 
    int I, J, K, i, j, k;
for(I = 0; I < sizePart/N; I += bs)
	{
		for(J = 0; J < N; J += bs)
		{
			abcblk = &abc[I*N + J];
			dcbblk = &dcb[I*N + J];
			for(K = 0; K < N; K += bs)
			{
				abblk = &ab[I*N + K];
				dcblk = &dc[I*N + K];
				bblk = &b[J*N + K];
				cblk = &c[J*N + K];
			
				for (i = 0; i < bs; i++)
				{
					for (j = 0; j < bs; j++)
					{
						for  (k = 0; k < bs; k++)
						{
							abcblk[i*N + j] += abblk[i*N + k] * cblk[j*N + k];
							dcbblk[i*N + j] += dcblk[i*N + k] * bblk[j*N + k];
						}
					}
				}
			}
		}
	}
}
// Calculamos p
void calcularP(double* abc, double* dcb, double minA, double maxD,double* P ,int N,int sizePart){
  int i,j;
  for (i=0; i<sizePart/N; i++){        
        for (j=0; j<N; j++){
            abc[i*N + j] *= maxD;
            dcb[i*N + j] *= minA;
            P[i*N + j]= abc[i*N + j] + dcb[i*N + j];
        }
    }
}
// Calculamos suma de P
void sumatoriaP(double *p, double *suma, int N,int sizePart){
int i, j;
double s=0;
  for(i = 0; i < sizePart/N; i++)
  {
      for(j = 0; j < N; j++)
      {
        s += p[i*N + j];
      }
  }
  *suma = s;
}
void calcularR(double *p, double prom, double *r, int N, int sizePart){
int i, j;
  for(i = 0; i < sizePart/N; i++)
  {
      for(j = 0; j < N; j++)
      {
        r[i*N + j]= p[i*N + j] *prom;
      }
  }
}

void Proceso0(int N, int T, int bs, int sizePart){
        double *A, *B, *C, *D;          //Matrices iniciadas con valores
	    double *AB, *DC, *ABC, *DCB;    //matrices temporales
	    double *P, *R;                  //Matrices que utilizan ecuaciones para inicializarse
	    
	    //Buffers
	    double *A_buf, *D_buf;
	    double maxD,maxlocal = DBL_MIN; //MaxD
	    double minA,minlocal = DBL_MAX; //MinA
	    double promP = 0;
	    double suma = 0, sumaP;
	    int i, j, k;
	    int I, J, K;
	    char exito=1;
	    int miID; int cantidadDeProcesos;
	    int sizeMatrix = N*N; // Cantidad total de datos matriz 	
        
        //Alocacion de memoria
	    A = (double *) malloc(sizeMatrix*sizeof(double));
	    D = (double *) malloc(sizeMatrix*sizeof(double));
	    B = (double *) malloc(sizeMatrix*sizeof(double));
	    C = (double *) malloc(sizeMatrix*sizeof(double));
	    AB = (double *) malloc(sizeMatrix*sizeof(double));
	    DC = (double *) malloc(sizeMatrix*sizeof(double));
	    ABC = (double *) malloc(sizeMatrix*sizeof(double));
	    DCB = (double *) malloc(sizeMatrix*sizeof(double));
        P = (double *) malloc(sizeMatrix*sizeof(double));
        R = (double *) malloc(sizeMatrix*sizeof(double));
        A_buf= (double *)malloc(sizeof(double)*sizePart); 
        D_buf= (double *)malloc(sizeof(double)*sizePart); 
        
         // Inicializacion
        init(A, B, C, D, AB, DC, ABC, DCB, P, R, N, sizePart);

        double t0 = dwalltime();
        // Distribuye los elementos de la matriz A y D
	    MPI_Scatter(A,sizePart,MPI_DOUBLE,A_buf,sizePart,MPI_DOUBLE,0,MPI_COMM_WORLD);
  	    MPI_Scatter(D,sizePart,MPI_DOUBLE,D_buf,sizePart,MPI_DOUBLE,0,MPI_COMM_WORLD);

  	    //Envio la matriz B y C completas a todos los procesos
        MPI_Bcast(B, sizeMatrix, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(C, sizeMatrix, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Llamo a multiplicar AB CD y obtener MinA MaxD de cada proceso.
        mult1(A_buf, B, C, D_buf, AB, DC, &minlocal, &maxlocal, N, bs, sizePart);

        //Se obtiene el MinA y MaxD total.
        MPI_Allreduce(&minlocal, &minA, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&maxlocal, &maxD, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        // Llamo a multiplicar ABC DCB y obtener MinA MaxD de cada proceso.
        mult2(AB, C, DC, B, ABC, DCB, N, bs, sizePart);
        //Se obtiene ABC y DCB Total.

        //Calculamos P y en gather combinamos todos
        calcularP(ABC, DCB, minA, maxD, P, N, sizePart);
        MPI_Gather(P, sizePart, MPI_DOUBLE, P, sizePart, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        //Calculamos la suma de P y en allreduce combinamos todos
        sumatoriaP(P, &suma, N, sizePart);
        MPI_Allreduce(&suma, &sumaP, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // Calculamos promedio de P
        promP=sumaP/(N*N);

        //Obtenemos R=PromP * P y en gather combinamos todos
        calcularR(P, promP, R, N, sizePart);
        MPI_Gather(R, sizePart, MPI_DOUBLE, R, sizePart, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        double t1 = dwalltime();
    
        // Tiempo de ejecucion
        printf("Tiempo en segundos: %f\n", t1-t0);

        // verificacion
        printf("Maximo de D es %f\n" ,maxD);
        printf("Min de A es %f\n" ,minA);
        printf("PromP %f\n", promP);
		printf("Validando matrices P y R.\n");
        // Validacion
        double validarP = N*N*A[0]*B[0]*C[0]*D[0]*2;
        double validarR = validarP*validarP;
        for (i = 0; i < N; i++)
        {
          for (j = 0; j < N; j++)
          {
            if (P[i*N + j] != validarP)
            {
              printf("Error en %d, %d, valor: %f\n", i, j, P[i*N + j]);
            if(exito)
		        exito=0;
            }
            if (R[i*N + j] != validarR)
            {
              printf("Error en %d, %d, valor: %f\n", i, j, R[i*N + j]);
            if(exito)
	          exito=0;
            }
          }
        }
        if (exito)
            printf("Validacion exitosa.\n");
        else
             printf("Validacion fallida.\n");
        
        // liberar memoria      
        free(A);
        free(B);
        free(C);
        free(D);
        free(P);
        free(R);
        free(AB);
        free(DC);
        free(A_buf);
        free(D_buf);
        free(ABC);
        free(DCB);
}

void Proceso1(int N, int T, int bs, int sizePart){
        double *A_buf, *B, *C, *D_buf; //Matrices iniciadas con valores
	    double *AB, *DC, *ABC, *DCB;		//matrices temporales
	    double *P, *R; //Matrices que utilizan ecuaciones para inicializarse
	    double *vacio;
	    //Buffers
	    double maxD,maxlocal = DBL_MIN; //MaxD
	    double minA,minlocal = DBL_MAX; //MinA
	    double promP = 0;
	    double suma = 0, sumaP;
	    int i, j, k;
	    int I, J, K;
	    int miID;
	    int sizeMatrix = N*N; // Cantidad total de datos matriz 	
        
        //Alocacion de memoria
	    A_buf= (double*)malloc(sizeof(double)*sizePart); 
        D_buf= (double*)malloc(sizeof(double)*sizePart);
	    B = (double *) malloc(sizeMatrix*sizeof(double));
	    C = (double *) malloc(sizeMatrix*sizeof(double));
	    AB = (double*)malloc(sizeof(double)*sizePart);
	    DC = (double*)malloc(sizeof(double)*sizePart);
	    ABC = (double*)malloc(sizeof(double)*sizePart);
	    DCB = (double*)malloc(sizeof(double)*sizePart);
        P = (double*)malloc(sizeof(double)*sizePart);
        R = (double*)malloc(sizeof(double)*sizePart);
        
        // Distribuye los elementos de la matriz A y D
	    MPI_Scatter(vacio,0,MPI_DOUBLE,A_buf,sizePart,MPI_DOUBLE,0,MPI_COMM_WORLD);
  	    MPI_Scatter(vacio,0,MPI_DOUBLE,D_buf,sizePart,MPI_DOUBLE,0,MPI_COMM_WORLD);

  	    //Envio la matriz B y C completas a todos los procesos
        MPI_Bcast(B, sizeMatrix, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(C, sizeMatrix, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Llamo a multiplicar AB CD y obtener MinA MaxD de cada proceso.
        mult1(A_buf, B, C, D_buf, AB, DC, &minlocal, &maxlocal, N, bs, sizePart);

        //Se obtiene el MinA y MaxD total.
        MPI_Allreduce(&minlocal, &minA, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&maxlocal, &maxD, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        // Llamo a multiplicar ABC DCB y obtener MinA MaxD de cada proceso.
        mult2(AB, C, DC, B, ABC, DCB, N, bs, sizePart);
        //Se obtiene ABC y DCB Total.

        //Calculamos P y en gather combinamos todos
        calcularP(ABC, DCB, minA, maxD, P, N, sizePart);
        MPI_Gather(P, sizePart, MPI_DOUBLE, vacio, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        //Calculamos la suma de P y en allreduce combinamos todos
        sumatoriaP(P, &suma, N, sizePart);
        MPI_Allreduce(&suma, &sumaP, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // Calculamos promedio de P
        promP=sumaP/(N*N);

        //Obtenemos R=PromP * P y en gather combinamos todos
        calcularR(P, promP, R, N, sizePart);
        MPI_Gather(R, sizePart, MPI_DOUBLE, vacio, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        // liberar memoria  
        free(A_buf);
        free(D_buf);    
        free(C);
        free(B);
        free(AB);
        free(DC);
        free(ABC);
        free(DCB);
        free(P);
        free(R);      
}
int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);

    int bs;
    int N;
	int miID; 
    int cantidadDeProcesos;
    int sizeMatrix;	
	int sizePart; // Cantidad de elementos por proceso
  
	MPI_Comm_rank(MPI_COMM_WORLD,&miID);
	MPI_Comm_size(MPI_COMM_WORLD,&cantidadDeProcesos);

     // Chequeo de parametros
	if ( (argc != 3) || ((N = atoi(argv[1])) < 1) || ((bs = atoi(argv[2])) <= 0) || ((N % bs) != 0) ){
		printf("Error en los parametros. Usar: ./%s N (minimo 1024), BS (N debe ser multiplo de BS)\n", argv[0]);
		exit(1);
	}
    sizeMatrix= N*N;
 	sizePart=sizeMatrix/cantidadDeProcesos; // Cantidad de elementos por bloque x Cantidad de bloques por proceso
    

    //computo
    if (miID==0){
        printf("Realizando operaciones con matrices de %d x %d en bloques de %d x %d con %d procesos\n", N, N, bs, bs, cantidadDeProcesos);
        Proceso0(N, cantidadDeProcesos, bs, sizePart);

    }
    else if (miID >= 1 && miID < cantidadDeProcesos){
        
        Proceso1(N, cantidadDeProcesos, bs, sizePart);

    }		 

	MPI_Finalize();
 
	return 0;

}   
    
    



















