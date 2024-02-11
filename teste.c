/*
 * teste.c
 *
 *  Created on: 3 de jan de 2021
 *      Author: frederico
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "matriz.h"

int main (){
	Matriz *A , *B , *C;
	int i,j,l,m,p,q;
	int passou;

	int N = 5000;

	printf ("Iniciando os testes da TAD matriz\n\n");

	printf ("Testando a criação da matriz e definição de elementos:\n");

	A = matrixCreate(N,N);
	l = matrixGetNumOfLines(A);
	m = matrixGetNumOfColumns(A);

	B = matrixCreate(N,N);
	p = matrixGetNumOfLines(B);
	q = matrixGetNumOfColumns(B);

	for (i=0;i<l;i++){
		for (j=0;j<m;j++){
			matrixSetElem(A,i,j,1);
		}
	}

	for (i=0;i<p;i++){
		for (j=0;j<q;j++){
			matrixSetElem(B,i,j,2);
		}
	}

	passou = 1;

	for (i=0;i<l;i++){
		for (j=0;j<m;j++){
			if (matrixGetElem(A,i,j)!=1){
				passou = 0;
			}
		}
	}

	for (i=0;i<l;i++){
		for (j=0;j<m;j++){
			if (matrixGetElem(B,i,j)!=2){
				passou = 0;
			}
		}
	}
	
	if (passou == 0){
		printf ("Falhou ...\n\n");
	} else {
		printf ("Passou ...\n\n");
	}

	printf ("Testando a adição de duas matrizes:\n");

	C = matrixAdd(A,B);

	passou = 1;

	for (i=0;i<l;i++){
		for (j=0;j<m;j++){
			//printf ("%d , %d = %f\n", i, j,creal(matrixGetElem(C,i,j)));
			if (matrixGetElem(C,i,j)!= 3){
				printf ("%d , %d = %f\n", i, j,matrixGetElem(C,i,j));
				passou = 0;
			}
		}
	}

	if (passou == 0){
		printf ("Falhou ...\n\n");
	} else {
		printf ("Passou ...\n\n");
	}

/*

	printf ("Testando a subtração de duas matrizes:\n");

	C = matrixSub(A,B);

	passou = 1;

	for (i=0;i<l;i++){
		for (j=0;j<m;j++){
			if (matrixGetElem(C,i,j)!= (i+j)-(i*j)){
				passou = 0;
			}
		}
	}

	if (passou == 0){
		printf ("Falhou ...\n\n");
	} else {
		printf ("Passou ...\n\n");
	}

	matrixSetElem(A, 0, 0, 2);
	matrixSetElem(A, 0, 1, 1);
	matrixSetElem(A, 1, 0, 5);
	matrixSetElem(A, 1, 1, 3);

	matrixSetElem(B, 0, 0, 3);
	matrixSetElem(B, 0, 1, -1);
	matrixSetElem(B, 1, 0, -5);
	matrixSetElem(B, 1, 1, 2);

	C = matrixMult(A, B);

	passou = 1;

	for (i=0;i<l;i++){
		for (j=0;j<m;j++){
			if (i == j){
				if (matrixGetElem(C,i,j)!=1)
					passou = 0;
			} else {
				if (matrixGetElem(C,i,j)!=0)
					passou = 0;
			}
		}
	}

	printf ("Testando a multiplicação de duas matrizes:\n");
	if (passou == 0){
		printf ("Falhou ...\n\n");
	} else {
		printf ("Passou ...\n\n");
	}


	matrixFree(C);

	matrixCreate(4,4);

	for (i=0;i<2;i++){
		for (j=0;j<2;j++){
			matrixSetElem(A,i,j,1);
			matrixSetElem(B,i,j,1);
		}
	}

	C = matrixKronProd(A,B);
	p = matrixGetNumOfLines(C);
	q = matrixGetNumOfColumns(C);

	passou = 1;

	for (i=0;i<p;i++){
		for (j=0;j<q;j++){
			if (matrixGetElem(C,i,j)!= 1){
				passou = 0;
			}
		}
	}

	printf ("Testando o produto de kronecker:\n");
	if (passou == 0){
		printf ("Falhou ...\n\n");
	} else {
		printf ("Passou ...\n\n");
	}

	printf ("Testando o cálculo do traço:\n");

	if (matrixTrace(C)==4){
		printf ("Passou ...\n\n");
	} else {
		printf ("Falhou ...\n\n");
	}

	for (i=0;i<2;i++){
		for (j=0;j<2;j++){
			matrixSetElem(A,i,j,i+j*I);
		}
	}

	B = matrixAdjoint(A);

	passou = 1;

	for (i=0;i<2;i++){
		for (j=0;j<2;j++){
			if (matrixGetElem(B,i,j)!= j - i*I){
				passou = 0;
			}
		}
	}

	printf ("Testando o cálculo da matriz adjunta:\n");
	if (passou == 0){
		printf ("Falhou ...\n\n");
	} else {
		printf ("Passou ...\n\n");
	}

	B = matrixAdjoint(B);

	printf ("Testando o cálculo da distância de Hilbert Schmidt:\n");
	if (matrixHilbertSchmidtDistance(A,B)==0){
		printf ("Passou ...\n\n");
	} else {
		printf ("Falhou ...\n\n");
	}

	*/

	printf ("Fim dos testes !!!\n\n");

	return (0);
}




