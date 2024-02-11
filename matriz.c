/*
 * matriz.c
 *
 *  Created on: 2 de out de 2020
 *      Author: frederico
 */
# include <stdio.h>
# include <stdlib.h>
# include <complex.h>
# include <math.h>
# include "matriz.h"

struct matrix {
	int lin;
	int col;
	complex **v;
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
					matrixSetElem (C, inicioLinha + k  , inicioColuna + l , (complex) (matrixGetElem (A,i,j) * matrixGetElem(B,k,l)));
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
	complex result;
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
	mat->v = (complex**) malloc(m * sizeof(complex*));
	for (i = 0; i < m; i++) {
		mat->v[i] = (complex*) malloc(n * sizeof(complex));
	}
	return mat;
}

Matriz* matrixAdd(Matriz *A, Matriz *B) {
	int m, n, i, j;
	complex plus;
	Matriz *C;
	m = A->lin;
	n = A->col;
	C = matrixCreate(m, n);
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			plus = (matrixGetElem(A, i, j)) + (matrixGetElem(B, i, j));
			matrixSetElem(C, i, j, plus);
		}
	}
	return (C);
}

Matriz* matrixSub(Matriz *A, Matriz *B) {
	int m, n, i, j;
	complex plus;
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
	complex result;
	Matriz *C;
	m = A->lin;
	n = A->col;
	C = matrixCreate(m, n);
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			result = (matrixGetElem(A, i, j));
			result = creal(result) - cimag(result)*I;
			matrixSetElem(C, j, i, result);
		}
	}
	return (C);
}


void matrixFree(Matriz *mat) {
	int i;
	for (i = 0; i < mat->lin; i++) {
		free(mat->v[i]);
	}
	free(mat->v);
	free(mat);
}

complex matrixGetElem(Matriz *mat, int i, int j) {
	if (i < 0 || i >= mat->lin || j < 0 || j >= mat->col) {
		printf("Invalid Access!\n");
		exit(1);
	}

	return (mat->v[i][j]);
}

void matrixSetElem(Matriz *mat, int i, int j, complex v) {
	if (i < 0 || i >= mat->lin || j < 0 || j >= mat->col) {
		printf("Impossible to assign the value!\n");
		exit(1);
	}
	mat->v[i][j] = v;
}

int matrixGetNumOfLines(Matriz *mat) {
	return (mat->lin);
}

int matrixGetNumOfColumns(Matriz *mat) {
	return (mat->col);

}

complex matrixTrace (Matriz *A){
	int i;
	complex trace = 0;
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



