#include <stdio.h>
#include <stdlib.h>  
#include <sys/time.h>
#include <float.h>
#include <limits.h>

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

/*La siguiente funcion realiza las siguientes operaciones:
AB
DC
ABC
DCB
MaxD
MinA
*/
void multibloques(double *a, double *b, double *c, double *d, double *e, double *f, double *g, double *h, double *minA, double *maxD, int n, int bs){
	double *ablk, *bblk, *cblk, *dblk, *eblk, *fblk, *gblk, *hblk;
	double minlocal=DBL_MAX, maxlocal=DBL_MIN;
	int I, J, K;    
	int i, j, k;
	for(I = 0; I < n; I += bs)
	{
		for(J = 0; J < n; J += bs)
		{
			eblk = &e[I*n + J];
			gblk = &g[I*n + J];
			for(K = 0; K < n; K += bs)
			{
				ablk = &a[I*n + K];
				bblk = &b[J*n + K];
				cblk = &c[J*n + K];
				dblk = &d[I*n + K];
				
				for (i = 0; i < bs; i++)
				{
					for (j = 0; j < bs; j++)
					{
						for  (k = 0; k < bs; k++)
						{
							//AB
							eblk[i*n + j] += ablk[i*n + k] * bblk[j*n + k];
							//DC
							gblk[i*n + j] += dblk[i*n + k] * cblk[j*n + k];
						}
						//MinA
						if (ablk[i*n + j] < minlocal)
						{
							minlocal = ablk[i*n + j];
						}
						//MaxD
						if (dblk[i*n + j] > maxlocal)
						{
							maxlocal = dblk[i*n + j];
						}
					}
				}
			}
		}
	}
	for(I = 0; I < n; I += bs)
	{
		for(J = 0; J < n; J += bs)
		{
			fblk = &f[I*n + J];
			hblk = &h[I*n + J];
			for(K = 0; K < n; K += bs)
			{
			eblk = &e[I*n + K];
			gblk = &g[I*n + K];
			cblk = &c[J*n + K];
			bblk = &b[J*n + K];
			for (i = 0; i < bs; i++)
				{
					for (j = 0; j < bs; j++)
					{
						for  (k = 0; k < bs; k++)
						{
							//ABC
							fblk[i*n + j] += eblk[i*n + k] * cblk[j*n + k];
							//DCB
							hblk[i*n + j] += gblk[i*n + k] * bblk[j*n + k];
						}
					}
				}
			}
		}
	}
	*maxD = maxlocal;
	*minA = minlocal;
}

/*La siguiente funcion realiza las siguientes operaciones:
MaxD(ABC)
MinA(DCB)
P = MaxD(ABC) + MinA(DCB)
PromP
*/
void sumabloques(double *f, double *h, double *p, double minA, double maxD, double *promP, int n){
	double suma = 0;
	int I, J, K; 
	for(I = 0; I < n; I++)
	{
		for(J = 0; J < n; J++)
		{

			//P= MaxD(ABC) + MinA(DCB)
			p[I*n + J] = (f[I*n + J] * maxD) + (h[I*n + J] * minA);
			suma += p[I*n + J];
		}
	}
	*promP = (suma/(n*n));
}
/*La siguiente funcion realiza las siguientes operaciones:
R= PromP(P)
*/
void calcularR(double *p, double *r, double promP, int n){
	double suma = 0;
	int I, J, K; 	
	for(I = 0; I < n; I++)
	{
		for(J = 0; J < n; J++)
		{

			//R= PromP(P)
			r[I*n + J] = (p[I*n + J] * promP) ;
		}
	}
}


// Multiplicación de matrices por bloques
void matmulblks(double *a, double *b, double *c, int n, int bs){
double *ablk, *bblk, *cblk;
int I, J, K;    
int i, j, k; 
 
  for(I = 0; I < n; I += bs)
  {
    for(J = 0; J < n; J += bs)
    {
		cblk = &c[I*n + J];
		for(K = 0; K < n; K += bs)
		{
		ablk = &a[I*n + K];
		bblk = &b[J*n + K];
		
		for (i = 0; i < bs; i++)
			{
				for (j = 0; j < bs; j++)
				{
					for  (k = 0; k < bs; k++)
					{
					cblk[i*n + j] += ablk[i*n + k] * bblk[j*n + k];
					}
				}
			}
		}
    }
  }
}

// Hallar el minimo de matrices por bloques
double matminblks(double *a, int n, int bs){
double *ablk, minlocal=DBL_MAX;
int I, J;    
int i, j; 
 
  for(I = 0; I < n; I += bs)
  {
    for(J = 0; J < n; J += bs)
    {
		ablk = &a[I*n + J];		
		for (i = 0; i < bs; i++)
		{
			for (j = 0; j < bs; j++)
			{
				if (ablk[i*n + j] < minlocal)
				{
					minlocal = ablk[i*n + j];
				}
			}
		}
    }
  }
  return minlocal;
}

// Hallar el maximo de matrices por bloques
double matmaxblks(double *a, int n, int bs){
double *ablk, maxlocal=DBL_MIN;
int I, J;    
int i, j; 
 
  for(I = 0; I < n; I += bs)
  {
    for(J = 0; J < n; J += bs)
    {
		ablk = &a[I*n + J];		
		for (i = 0; i < bs; i++)
		{
			for (j = 0; j < bs; j++)
			{
				if (ablk[i*n + j] > maxlocal)
				{
					maxlocal = ablk[i*n + j];
				}
			}
		}
    }
  }
  return maxlocal;
}

// Hallar el promedio de una matriz por bloques
double matpromblks(double *a, int n, int bs){
double *ablk, suma=0;
int I, J;    
int i, j; 
 
  for(I = 0; I < n; I += bs)
  {
    for(J = 0; J < n; J += bs)
    {
		ablk = &a[I*n + J];		
		for (i = 0; i < bs; i++)
		{
			for (j = 0; j < bs; j++)
			{
				suma += ablk[i*n + j];
			}
		}
    }
  }
  return (suma/(n*n));
}

// Multiplicación de una matriz por un escalar por bloques
void matescblks(double *a, double *b, double esc, int n, int bs){
double *ablk, *bblk;
int I, J;    
int i, j; 
 
  for(I = 0; I < n; I += bs)
  {
    for(J = 0; J < n; J += bs)
    {
		bblk = &b[I*n + J];
		ablk = &a[I*n + J];
		for (i = 0; i < bs; i++)
		{
			for (j = 0; j < bs; j++)
			{
			
				bblk[i*n + j] = ablk[i*n + j] * esc;
			}
		}
    }
  }
}

// Suma de matrices por bloques
void matsumablks(double *a, double *b, double *c, int n, int bs){
double *ablk, *bblk, *cblk;
int I, J, K;    
int i, j, k; 
 
  for(I = 0; I < n; I += bs)
  {
    for(J = 0; J < n; J += bs)
    {
		cblk = &c[I*n + J];
		for(K = 0; K < n; K += bs)
		{
		ablk = &a[I*n + K];
		bblk = &b[J*n + K];
		
		for (i = 0; i < bs; i++)
			{
				for (j = 0; j < bs; j++)
				{
					for  (k = 0; k < bs; k++)
					{
					cblk[i*n + j] = ablk[i*n + k] + bblk[j*n + k];
					}
				}
			}
		}
    }
  }
}

int main(int argc, char *argv[]){
  double *A, *B, *C, *D; //Matrices iniciadas con valores aleatorios
  double *P, *R; //Matrices que utilizan ecuaciones para inicializarse
  double *E, *F, *G, *H; //Matrices temporales E=AB, F=EC, G=DC, H=GB, I=MaxD(F), J=MinA(H)
  double maxD = DBL_MAX; //MaxD
  double minA = DBL_MIN; //MinA
  double promP = 0;
  int n, bs, i, j;
  char exito=1;

  double timetick;

  // Chequeo de parámetros
	if ( (argc != 3) || ((n = atoi(argv[1])) < 1) || ((bs = atoi(argv[2])) <= 0) || ((n % bs) != 0)){
		printf("Error en los parámetros. Usar: ./%s N (minimo 1024), BS (N debe ser multiplo de BS)\n", argv[0]);
		exit(1);
	}

  // Alocar memoria
	A = (double *) malloc(n*n*sizeof(double));
	B = (double *) malloc(n*n*sizeof(double));
	C = (double *) malloc(n*n*sizeof(double));
	D = (double *) malloc(n*n*sizeof(double));
	E = (double *) malloc(n*n*sizeof(double));
	F = (double *) malloc(n*n*sizeof(double));
	G = (double *) malloc(n*n*sizeof(double));
	H = (double *) malloc(n*n*sizeof(double));
	P = (double *) malloc(n*n*sizeof(double));
	R = (double *) malloc(n*n*sizeof(double));

  // Inicializacion
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			A[i*n + j] = 1.0; //Inicializo por fila
			B[j*n + i] = 2.0; //Inicializo por columna
			C[j*n + i] = 3.0; //Inicializo por columna
			D[i*n + j] = 4.0; //Inicializo por fila
			
			E[i*n + j] = 0.0;
			F[i*n + j] = 0.0;
			G[i*n + j] = 0.0;
			H[i*n + j] = 0.0;
			P[i*n + j] = 0.0;
			R[i*n + j] = 0.0;
		}
	}


	printf("Realizando operaciones con matrices de %d x %d en bloques de %d x %d\n", n, n, bs, bs);
  
  timetick = dwalltime();
 
  multibloques(A, B, C, D, E, F, G, H, &minA, &maxD, n, bs);
  sumabloques(F, H, P, minA, maxD, &promP, n);
  calcularR(P, R, promP, n);
  
  
  double totalTime = dwalltime() - timetick;
	printf("Tiempo en segundos %f\n",totalTime);
  // Validando
  printf("Validando P y R...\n");
  //La validacion asume que la matriz esta cargada con el mismo numero
  double validarP = n*n*A[0]*B[0]*C[0]*D[0]*2;
  double validarR = validarP*validarP;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      if (P[i*n + j] != validarP)
      {
        printf("Error en P[%d, %d], valor: %f\n", i, j, P[i*n + j]);
		if(exito)
				exito=0;
      }
      if (R[i*n + j] != validarR)
      {
        printf("Error en R[%d, %d], valor: %f\n", i, j, R[i*n + j]);
		if(exito)
			exito=0;
      }
    }
  }
  if (exito)
	printf("Validacion exitosa.\n");
  else
	  printf("Validacion fallida.\n");
	//printf("Tiempo en segundos %f\n",totalTime);
 

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
