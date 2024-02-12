/*
 * matriz.c
 *
 *  Created on: 2 de out de 2020
 *      Author: frederico
 */
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <omp.h>
# include "matriz.h"


struct matrix {
	int lin;
	int col;
	double *v;
};


Matriz* matrixKronProd(Matriz *A, Matriz *B) {
	Matriz *C;
	int m, n, p, q;
	int i, j, k, l;
	int inicioLinha, inicioColuna;
	m = A->lin;
	n = A->col;
	p = B->lin;
	q = B->col;

	C = matrixCreate(m*p, n*q);

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			inicioLinha = i*p;
			inicioColuna = j*q;
			for (k = 0; k < p; k++) {
				for (l = 0; l < q; l++) {
					matrixSetElem (C, inicioLinha + k  , inicioColuna + l , (double) (matrixGetElem (A,i,j) * matrixGetElem(B,k,l)));
				}
			}
		}
	}
	return (C);
}


Matriz* matrixMult(Matriz *A, Matriz *B) {
	Matriz *C;
	int m, n, p, q;
	int i, j, k;
	double result;
	m = A->lin;
	n = A->col;
	p = B->lin;
	q = B->col;
	if (n != p) {
		printf("Incompatible types: A is %d x %d and B is %d x %d\n",m,n,p,q);
		exit(1);
	}
	C = matrixCreate(m, q);

	for (i = 0; i < m; i++) {
		for (j = 0; j < p; j++) {
			result = 0;
			for (k = 0; k < n; k++) {
				result = result + matrixGetElem(A, i, k) * matrixGetElem(B, k, j);
			}
			matrixSetElem(C, i, j, result);
		}
	}

	return (C);
}

Matriz* matrixCreate(int m, int n) {
	int i;
	Matriz *mat = (Matriz*) malloc(sizeof(Matriz));
	if (mat == NULL) {
		printf("Insufficient Memory!\n");
		exit(1);
	}
	mat->lin = m;
	mat->col = n;
	//mat->v = (double**) malloc(m*n* sizeof(double*));

	mat->v = (double*) malloc(m * n * sizeof(double));

	/*
	for (i = 0; i < m; i++) {
		mat->v[i] = (double*) malloc(n * sizeof(double));
	}
	*/
	return mat;
}

Matriz* matrixAdd(Matriz *A, Matriz *B) {
	int m, n;
	double plus;
	Matriz *C;
	m = A->lin;
	n = A->col;
	C = matrixCreate(m, n);

	int start, end, n_threads;

    double start_time = omp_get_wtime();
    #pragma omp parallel private(start, end, n_threads)
    {
		int i,j;
        n_threads = omp_get_num_threads();
        int thread_id = omp_get_thread_num();
        int rows_per_thread = m / n_threads;
        int remainder = m % n_threads;

        start = rows_per_thread*thread_id;
        end = start + rows_per_thread -1;

        if(thread_id+1 <= remainder){
            start = start+thread_id;
            end = end+thread_id+1;
        }
        else{
            start = start+remainder;
            end = end+remainder;
        }

		//printf ("%d --> %d , %d\n", thread_id, start, end);

		for (i = start; i <= end; i++) {
			for (j = 0; j < n; j++) {				
				matrixSetElem(C, i, j, (matrixGetElem(A, i, j)) + (matrixGetElem(B, i, j)));
			}
		}
	}
	double stop_time = omp_get_wtime();
	printf ("Tempo de execução: %f segundos\n\n", (stop_time - start_time));

	return (C);
}

Matriz* matrixSub(Matriz *A, Matriz *B) {
	int m, n, i, j;
	double plus;
	Matriz *C;
	m = A->lin;
	n = A->col;
	C = matrixCreate(m, n);
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			plus = (matrixGetElem(A, i, j)) - (matrixGetElem(B, i, j));
			matrixSetElem(C, i, j, plus);
		}
	}
	return (C);
}

Matriz* matrixAdjoint(Matriz *A) {
	int m, n, i, j;
	double result;
	Matriz *C;
	m = A->lin;
	n = A->col;
	C = matrixCreate(m, n);
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			result = (matrixGetElem(A, i, j));
			//result = creal(result) - cimag(result)*I;
			matrixSetElem(C, j, i, result);
		}
	}
	return (C);
}


void matrixFree(Matriz *mat) {
	int i;
	free(mat->v);
	//free(mat);
}

double matrixGetElem(Matriz *mat, int i, int j) {
    if (i >= 0 && i < mat->lin && j >= 0 && j < mat->col) {
        return mat->v[i * mat->col + j]; // Acesso ajustado para array unidimensional
    } else {
        printf("Index out of bounds\n");
        exit(1);
    }
}


void matrixSetElem(Matriz *mat, int i, int j, double value) {
    if (i >= 0 && i < mat->lin && j >= 0 && j < mat->col) {
        mat->v[i * mat->col + j] = value; // Atribuição ajustada para array unidimensional
    } else {
        printf("Index out of bounds\n");
        exit(1);
    }
}

int matrixGetNumOfLines(Matriz *mat) {
	return (mat->lin);
}

int matrixGetNumOfColumns(Matriz *mat) {
	return (mat->col);

}

double matrixTrace (Matriz *A){
	int i;
	double trace = 0;
	if (matrixGetNumOfLines(A)!= matrixGetNumOfColumns(A)){
		printf ("Impossible to compute trace: not a square matrix");
		exit(1);
	}

	for (i=0;i<matrixGetNumOfLines(A);i++){
		trace = trace + matrixGetElem(A,i,i);
	}
	return(trace);
}



double matrixHilbertSchmidtDistance(Matriz *A, Matriz *B){
	double d,e;
	Matriz *M, *N, *AUX;
	M = matrixSub(A,B);
	AUX = matrixAdjoint(M);
	N = matrixMult(M,AUX);
	d = (double) matrixTrace(N);
	e = sqrt (d);
	matrixFree(AUX);
	matrixFree(M);
	matrixFree(N);
	return (e);
}



