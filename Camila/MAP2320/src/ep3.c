/*
 * Arquivo: ep.c
 * -------------
 * Este programa resolve a equação
 * de Poisson 2D no quadrado unitário.
 * O objetivo do programa é modelar
 * numericamente o campo eletrostático
 * entre duas placas paralelas.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*
 * Constantes
 * ----------
 * ...
 */

#define L 1.0
#define TOL 1.0e-5
#define MAX_ITER 100000
#define INI_EXP 3 // expoente inicial de 2 para malha mesh
#define FIN_EXP 9 // expoente final de 2 para malha mesh
#define ALPHA 7 // último algarismo do número usp
#define DENSITY_1 8.84*10e-12
#define DENSITY_2 75*10e-12

/* Protótipos das funções */

static void set_axis_mesh(double h, int n, double *x);
static double f_eval(double x, double y);
static double ue_eval(double x, double y);
static double permiss_eval(double x, double y);
static void set_dirichlet_bound(int n, double *x, double *y, double **u);
static void set_u_exact(int n, double *x, double *y, double **u);
static void set_init_cond(int n, double **u);
static int jacobi(double h, int n, double *x, double *y, double **u);
static int gauss_seidel(double h, int n, double *x, double *y, double **u);
static int sor(double h, int n, double *x, double *y, double **u);
static double norm_h_eval(double h, int n, double **u, double **u_prev);
static double norm_inf_eval(int n, int m, double **u, double **ue);
static double *vector_alloc(int n);
static double **matrix_alloc(int m, int n);

static void print_matrix(int n, int m, double **matrix); //debug: for gdb

/* Programa princial */
int main()
{
	int i, n, mesh, iter;
	double h, *x, *y, **u, **ue, norm, norm_prev, elapsed_t;
	clock_t start_t;

	/*
	 * O índice i é usado para selecionar os
	 * elementos do array norm, que armazena
	 * as normas sucessivas da diferença entre
	 * as matrizes u (aproximação) e ue (solução
	 * exata).
	 */
	//i = 0;

	for (mesh = INI_EXP; mesh <= FIN_EXP; mesh++) {
		n = pow(2, mesh);
		h = L / n;
	
		x = vector_alloc(n + 1);
		y = vector_alloc(n + 1);
		u = matrix_alloc(n + 1, n + 1);
		ue = matrix_alloc(n + 1, n + 1);

		set_axis_mesh(h, n, x);
		set_axis_mesh(h, n, y);

		set_dirichlet_bound(n, x, y, u);
		set_init_cond(n, u);

		start_t = clock();
		//iter = jacobi(h, n, x, y, u);
		//iter = gauss_seidel(h, n, x, y, u);
		iter = sor(h, n, x, y, u);
		elapsed_t = (double) (clock() - start_t) / CLOCKS_PER_SEC; 
		set_u_exact(n, x, y, ue);
		norm = norm_inf_eval(n, n, u, ue);
		
		/*
		 * Imprime tabela de convergência com tempo
		 * de resolução
		 */
		if (mesh == INI_EXP) {
			printf("%d\t%f\t-\t\t%d\t%f\n", n, norm, iter, elapsed_t);
		} else {
			printf("%d\t%f\t%f\t%d\t%f\n", n, norm, norm_prev / norm, iter, elapsed_t);
		}
		norm_prev = norm;

		free(x);
		free(y);
		free(u);
		free(ue);
	}

}

/*
 * Função: set_axis_mesh
 * Uso: set_axis_mesh(h, x);
 * -------------------------
 * Define o pontos da malha na coordenada x.
 */

static void set_axis_mesh(double h, int n, double *x)
{
	int i;

	for (i = 0; i <= n; i++) {
		x[i] = i*h;
	}
}

static void set_dirichlet_bound(int n, double *x, double *y, double **u)
{
	int i, j;

	for (i = 0; i <= n; i++) {
		// Tarefa 3.1
		u[i][0] = ue_eval(x[i], y[0]);
		u[i][n] = ue_eval(x[i], y[n]);

		// Tarefa 3.2
		//u[i][0] = 0;
		//u[i][n] = 110;
	}
	for (j = 0; j <= n; j++) {
		// Tarefa 3.1
		u[0][j] = ue_eval(x[0], y[j]);
		u[n][j] = ue_eval(x[n], y[j]);

		// Tarefa 3.2
		//u[0][j] = 110 * sin(M_PI / 2 * y[j]);
		//u[n][j] = u[0][j];
	}
}

/*
 * Função: ue_eval
 * Uso: ue_ij = ue_eval(x, y);
 * ---------------------------
 * Retorna o valor da solulção exata
 * no ponto (x, y).
 */

static double ue_eval(double x, double y) {
	// Tarefa 3.1.a
	return (ALPHA * exp(x) * sin(y));

	// Tarefa 3.2.b
	//return (ALPHA * cos(M_PI * x) * sin(M_PI * y));
}

/*
 * Função: f_eval
 * Uso: result = f_eval(x, y);
 * ---------------------------
 * Retorna o termo forçante da equação
 * de Poisson no ponto (x, y).
 */

static double f_eval(double x, double y)
{
	// Tarefa 3.1.a
	return (0);

	// Tarefa 3.1.b
	//return ( 2 * ALPHA * M_PI * M_PI * cos(M_PI * x) * sin(M_PI * y) );

	// Tarefa 3.2.a
	//return (permiss_eval(x, y) / DENSITY_1);

	// Tarefa 3.2.b
	//return (permiss_eval(x, y) / DENSITY_2);
}

/*
 * Função: permiss_eval
 * Uso: result = permiss_eval(x, y);
 * ---------------------------------
 * Retorna a permissividade de um campo
 * elétrico bi-dimensional no ponto
 * (x, y).
 */

static double permiss_eval(double x, double y)
{
	// Tarefa 3.2.a
	//return (100 * 10e-12);

	// Tarefa 3.2.b
	//return ( 10 * sin(M_PI * (x + y)) * 10e-08 );
}

/*
 * Função: set_u_exact
 * Uso: set_u_exact(n, x, y, u);
 * -----------------------------
 * Define a solução exata nos pontos
 * internos da malha.
 */

static void set_u_exact(int n, double *x, double *y, double **ue)
{
	int i, j;

	for (i = 1; i < n; i++) {
		for (j = 1; j < n; j++) {
			ue[i][j] = ue_eval(x[i], y[j]);
		}
	}
}

/*
 * Função: jacobi
 * Uso: ite3r = jacobi(h, n, x, y, u);
 * -----------------------------------
 * Implementa o algoritmo de Jacobi
 * para o resolver o sistema linear
 * oriundo do método numérico para
 * aproximar a solução da equação
 * de Poisson.
 */

static int jacobi(double h, int n, double *x, double *y, double **u)
{
	int i, j, k;
	double norm, **u_prev;

	u_prev = matrix_alloc(n + 1, n + 1);

	for (i = 0; i <= n; i++) {
		for (j = 0; j <= n; j++) {
			u_prev[i][j] = u[i][j];
		}
	}

	for (k = 0; k < MAX_ITER; k++) {
		for (i = 1; i < n; i++) {
			for (j = 1; j < n; j++) {
				u[i][j] = 0.25 * (u_prev[i-1][j] + u_prev[i+1][j] +
								  u_prev[i][j-1] + u_prev[i][j+1] +
								  h*h * f_eval(x[i], y[j]));
			}
		}
		norm = norm_h_eval(h, n, u, u_prev);
		if (norm < TOL * h) { // método convergiu
			 free(u_prev);
			 return (k); 
		}
		for (i = 1; i < n; i++) {
			for (j = 1; j < n; j++) {
				u_prev[i][j] = u[i][j];
			}
		}
	}
	printf("Número máximo de iterações atingido.\n");
	free(u_prev);
	exit(1);
}

/*
 * Função: gauss_seidel
 * Uso: iter = gauss_seidel(h, n, x, y, u);
 * ----------------------------------------
 * Implementa o algoritmo de Gauss Seidel
 * para o resolver o sistema linear
 * oriundo do método numérico para
 * aproximar a solução da equação
 * de Poisson.
 */

static int gauss_seidel(double h, int n, double *x, double *y, double **u)
{
	int i, j, k;
	double norm, **u_prev;

	u_prev = matrix_alloc(n + 1, n + 1);

	for (i = 0; i <= n; i++) {
		for (j = 0; j <= n; j++) {
			u_prev[i][j] = u[i][j];
		}
	}

	for (k = 0; k < MAX_ITER; k++) {
		for (i = 1; i < n; i++) {
			for (j = 1; j < n; j++) {
				u[i][j] = 0.25 * (u[i-1][j] + u_prev[i+1][j] +
								  u[i][j-1] + u_prev[i][j+1] +
								  h*h * f_eval(x[i], y[j]));
			}
		}
		norm = norm_h_eval(h, n, u, u_prev);
		if (norm < TOL * h) { // método convergiu
			 free(u_prev);
			 return (k); 
		}
		for (i = 1; i < n; i++) {
			for (j = 1; j < n; j++) {
				u_prev[i][j] = u[i][j];
			}
		}
	}
	printf("Número máximo de iterações atingido.\n");
	free(u_prev);
	exit(1);
}

/*
 * Função: sor 
 * Uso: iter = sor(h, n, x, y, u);
 * -------------------------------
 * Implementa o algoritmo SOR
 * para o resolver o sistema linear
 * oriundo do método numérico para
 * aproximar a solução da equação
 * de Poisson.
 */

static int sor(double h, int n, double *x, double *y, double **u)
{
	int i, j, k;
	double norm, omega, **u_prev;

	u_prev = matrix_alloc(n + 1, n + 1);

	omega = 2 / (1 + sin(M_PI * h));

	for (i = 0; i <= n; i++) {
		for (j = 0; j <= n; j++) {
			u_prev[i][j] = u[i][j];
		}
	}

	for (k = 0; k < MAX_ITER; k++) {
		for (i = 1; i < n; i++) {
			for (j = 1; j < n; j++) {
				u[i][j] = (1 - omega)*u_prev[i][j] + 0.25*omega *
						  ( u[i-1][j] + u_prev[i+1][j] +
						    u[i][j-1] + u_prev[i][j+1] +
						    h*h * f_eval(x[i], y[j]) );
			}
		}
		norm = norm_h_eval(h, n, u, u_prev);
		if (norm < TOL * h) { // método convergiu
			 free(u_prev);
			 return (k); 
		}
		for (i = 1; i < n; i++) {
			for (j = 1; j < n; j++) {
				u_prev[i][j] = u[i][j];
			}
		}
	}
	printf("Número máximo de iterações atingido.\n");
	free(u_prev);
	exit(1);
}

/*
 * Função: set_init_cond
 * Uso: set_init_cond(u);
 * ----------------------
 * Define a condição inicial para a
 * solução aproximada.
 */

static void set_init_cond(int n, double **u)
{
	int i, j;

	for (i = 1; i < n; i++) {
		for (j = 1; j < n; j++) {
			u[i][j] = 0;
		}
	}
}

/*
 * Função: norm_h_eval
 * Uso: norm_h_eval(u, u_prev, n, h);
 * ----------------------------------
 * Esta função calcula a norma da diferença
 * entre os vetores u e u_prev, onde u e u_prev são
 * as aproximações calculadas na iteração k e
 * k - 1, respectivamente. A variável h é
 * a distância - a princípio, uniforme - entre
 * os pontos da malha.
 */

static double norm_h_eval(double h, int n, double **u, double **u_prev)
{
    int i, j;
    double result;

    result = 0;
    for (i = 1; i < n; i++) {
		for (j = 1; j < n; j++) {
			result += (u[i][j] - u_prev[i][j]) * (u[i][j] - u_prev[i][j]);
		}
	}
	result = sqrt(result) * h;
	return (result);
}

/*
 * Function: norm_inf_eval
 * Uso: norm_inf = norm_inf_eval(.); //TODO
 * ---------------------------------
 * Calcula a norma infinito das aproximações nas
 * sucessivas malhas, em relação à solução exata.
 */

static double norm_inf_eval(int n, int m, double **u, double **ue)
{
	int i, j;
	double norm_inf;

	norm_inf = 0;
	for (i = 1; i < n - 1; i++) {
		for (j = 1; j < m - 1; j++) {
			if (fabs(u[i][j] - ue[i][j]) > norm_inf ) {
				norm_inf = fabs(u[i][j] - ue[i][j]);
			}
		}
	}
	return (norm_inf);
}

/*
 * Function: vector_alloc
 * Usage: arr = vector_alloc(n);
 * ---------------------------
 * This function allocates memmory
 * for array arr with n positions.
 * The implementation handles
 * invalid negative n, and out of
 * memmory.
 */

static double *vector_alloc(int n)
{
	double *vec; /* ponteiro para o vetor */

	if (n < 1) {
		printf("Erro: Parâmetro inválido\n");
		return(NULL);
	}
	vec = malloc(n * sizeof(double));
	if (vec == NULL) {
		printf("Erro: Memória Insuficiente");
		return(NULL);
	}
	return vec;
}


/*
 * Function: matrix_alloc 
 * Usage: d_arr = matrix_alloc(m, n);
 * ---------------------------------
 * This function allocates memmory
 * for double array d_ar with m lines
 * and n columns.
 * The implementation handles
 * invalid negative n, and out of
 * memmory.
 */

static double **matrix_alloc(int m, int n)
{
	double **v;  /* ponteiro para a matriz */
	int   i;    /* variavel auxiliar      */
	if (m < 1 || n < 1) { /* verifica parametros recebidos */
		printf ("** Erro: Parametro invalido **\n");
		return (NULL);
	}
	/* aloca as linhas da matriz */
	v = (double **) calloc (m, sizeof(double *));
	if (v == NULL) {
		printf ("** Erro: Memoria Insuficiente **");
		return (NULL);
	}
	/* aloca as colunas da matriz */
	for (i = 0; i < m; i++) {
		v[i] = (double*) calloc (n, sizeof(double));
		if (v[i] == NULL) {
			printf ("** Erro: Memoria Insuficiente **");
			return (NULL);
		}
	}
	return (v); /* retorna o ponteiro para a matriz */
}

/*
 * Funçao auxilar para o gdb - debug
 * ---------------------------------
 */

static void print_matrix(int n, int m, double **matrix)
{
	int i, j;

	for (i = 0; i <= n; i++) {
		for (j = 0; j <= m; j++) {
			printf("%f ", matrix[i][j]);
		}
		printf("\n");
	}
}
