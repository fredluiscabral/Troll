/*
 * decomp_v2.c
 *
 *  Created on: 14 de out de 2020
 *      Author: frederico
 */

/*
 * decomp.c
 *
 *  Created on: 7 de out de 2020
 *      Author: frederico
 */
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "matriz.h"

struct matrixName {
	char name[100];
	Matriz * matrix;
};

typedef struct columns {
	char sequenceOfMatrices[1000];
	Matriz * matrix;
} columns;



/*
int main() {


	unsigned long int i;
	int j,l;
	double N = 4; // U port dimension
	double k = (int) (log10(N) / log10(2)); // Column size
	int r=7; // Size of the Universal Set os Gates
	int nC; // Number of columns
	complex result;
	int base;
	Matriz *U, *Iden, *solutionAnt, *solution; //Target Matrix U and identity N x N
	int maxIter = 3;
	double hilSchDist, prec = 0.01;
	int * sizesOfUnivGatesSet;
	int numOfPartitions , extraSize;


	// Fourier Matrix 4x4

	U = matrixCreate(4,4);
	matrixSetElem(U,0,0,1);
	matrixSetElem(U,0,1,1);
	matrixSetElem(U,0,2,1);
	matrixSetElem(U,0,3,1);
	matrixSetElem(U,1,0,1);
	matrixSetElem(U,1,1,I);
	matrixSetElem(U,1,2,-1);
	matrixSetElem(U,1,3,-I);
	matrixSetElem(U,2,0,1);
	matrixSetElem(U,2,1,-1);
	matrixSetElem(U,2,2,1);
	matrixSetElem(U,2,3,-1);
	matrixSetElem(U,3,0,1);
	matrixSetElem(U,3,1,-I);
	matrixSetElem(U,3,2,-1);
	matrixSetElem(U,3,3,I);



	// Fourier Matrix 4x4
	U = matrixCreate(4,4);
	matrixSetElem(U,0,0,0);
	matrixSetElem(U,0,1,sqrt(2)/2);
	matrixSetElem(U,0,2,0);
	matrixSetElem(U,0,3,sqrt(2)/2);

	matrixSetElem(U,1,0,(1 + I)/2);
	matrixSetElem(U,1,1,0);
	matrixSetElem(U,1,2,(1 + I)/2);
	matrixSetElem(U,1,3,0);

	matrixSetElem(U,2,0,0);
	matrixSetElem(U,2,1,(1 + I)/2);
	matrixSetElem(U,2,2,0);
	matrixSetElem(U,2,3,-(1 + I)/2);

	matrixSetElem(U,3,0,(I*sqrt(2))/2);
	matrixSetElem(U,3,1,0);
	matrixSetElem(U,3,2,-(I*sqrt(2))/2);
	matrixSetElem(U,3,3,0);


	Iden = matrixCreate(4,4);
	matrixSetElem(Iden,0,0,1);
	matrixSetElem(Iden,0,1,0);
	matrixSetElem(Iden,0,2,0);
	matrixSetElem(Iden,0,3,0);
	matrixSetElem(Iden,1,0,0);
	matrixSetElem(Iden,1,1,1);
	matrixSetElem(Iden,1,2,0);
	matrixSetElem(Iden,1,3,0);
	matrixSetElem(Iden,2,0,0);
	matrixSetElem(Iden,2,1,0);
	matrixSetElem(Iden,2,2,1);
	matrixSetElem(Iden,2,3,0);
	matrixSetElem(Iden,3,0,0);
	matrixSetElem(Iden,3,1,0);
	matrixSetElem(Iden,3,2,0);
	matrixSetElem(Iden,3,3,1);


	// ============================== Universal Set ==============================

	sizesOfUnivGatesSet = malloc (3*sizeof(int));

	sizesOfUnivGatesSet[0] = 5;
	sizesOfUnivGatesSet[1] = 2;
	sizesOfUnivGatesSet[2] = 0;

	numOfPartitions = sizesOfUnivGatesSet[0] + sizesOfUnivGatesSet[1];

	struct matrixName *universalGates = malloc((numOfPartitions)*sizeof(struct matrixName));

	universalGates[0].matrix = matrixCreate(2, 2);
	strcpy (universalGates[0].name,"Id");
	matrixSetElem(universalGates[0].matrix, 0, 0, 1);
	matrixSetElem(universalGates[0].matrix, 0, 1, 0);
	matrixSetElem(universalGates[0].matrix, 1, 0, 0);
	matrixSetElem(universalGates[0].matrix, 1, 1, 1);

	universalGates[1].matrix = matrixCreate(2, 2);
	strcpy (universalGates[1].name,"H");
	matrixSetElem(universalGates[1].matrix, 0, 0, (1 / sqrt(2)) * (+1));
	matrixSetElem(universalGates[1].matrix, 0, 1, (1 / sqrt(2)) * (+1));
	matrixSetElem(universalGates[1].matrix, 1, 0, (1 / sqrt(2)) * (+1));
	matrixSetElem(universalGates[1].matrix, 1, 1, (1 / sqrt(2)) * (-1));

	universalGates[2].matrix = matrixCreate(2, 2);
	strcpy (universalGates[2].name,"T");
	matrixSetElem(universalGates[2].matrix, 0, 0, 1);
	matrixSetElem(universalGates[2].matrix, 0, 1, 0);
	matrixSetElem(universalGates[2].matrix, 1, 0, 0);
	matrixSetElem(universalGates[2].matrix, 1, 1, sqrt(2) / 2 + I * sqrt(2) / 2);

	universalGates[3].matrix = matrixCreate(2, 2);
	strcpy (universalGates[3].name,"Tg");
	matrixSetElem(universalGates[3].matrix, 0, 0, 1);
	matrixSetElem(universalGates[3].matrix, 0, 1, 0);
	matrixSetElem(universalGates[3].matrix, 1, 0, 0);
	matrixSetElem(universalGates[3].matrix, 1, 1, sqrt(2) / 2 - I * sqrt(2) / 2);

	universalGates[4].matrix = matrixCreate(2, 2);
	strcpy (universalGates[4].name,"X");
	matrixSetElem(universalGates[4].matrix, 0, 0, 0);
	matrixSetElem(universalGates[4].matrix, 0, 1, 1);
	matrixSetElem(universalGates[4].matrix, 1, 0, 1);
	matrixSetElem(universalGates[4].matrix, 1, 1, 0);

	// Define gates (4x4 columns) CNOT12 and CNOT21

	universalGates[5].matrix = matrixCreate(4, 4);
	strcpy (universalGates[5].name,"CNOT12");
	matrixSetElem(universalGates[5].matrix, 0, 0, 1);
	matrixSetElem(universalGates[5].matrix, 0, 1, 0);
	matrixSetElem(universalGates[5].matrix, 0, 2, 0);
	matrixSetElem(universalGates[5].matrix, 0, 3, 0);

	matrixSetElem(universalGates[5].matrix, 1, 0, 0);
	matrixSetElem(universalGates[5].matrix, 1, 1, 1);
	matrixSetElem(universalGates[5].matrix, 1, 2, 0);
	matrixSetElem(universalGates[5].matrix, 1, 3, 0);

	matrixSetElem(universalGates[5].matrix, 2, 0, 0);
	matrixSetElem(universalGates[5].matrix, 2, 1, 0);
	matrixSetElem(universalGates[5].matrix, 2, 2, 0);
	matrixSetElem(universalGates[5].matrix, 2, 3, 1);

	matrixSetElem(universalGates[5].matrix, 3, 0, 0);
	matrixSetElem(universalGates[5].matrix, 3, 1, 0);
	matrixSetElem(universalGates[5].matrix, 3, 2, 1);
	matrixSetElem(universalGates[5].matrix, 3, 3, 0);

	universalGates[6].matrix = matrixCreate(4, 4);
	strcpy (universalGates[6].name,"CNOT21");
	matrixSetElem(universalGates[6].matrix, 0, 0, 1);
	matrixSetElem(universalGates[6].matrix, 0, 1, 0);
	matrixSetElem(universalGates[6].matrix, 0, 2, 0);
	matrixSetElem(universalGates[6].matrix, 0, 3, 0);

	matrixSetElem(universalGates[6].matrix, 1, 0, 0);
	matrixSetElem(universalGates[6].matrix, 1, 1, 0);
	matrixSetElem(universalGates[6].matrix, 1, 2, 0);
	matrixSetElem(universalGates[6].matrix, 1, 3, 1);

	matrixSetElem(universalGates[6].matrix, 2, 0, 0);
	matrixSetElem(universalGates[6].matrix, 2, 1, 0);
	matrixSetElem(universalGates[6].matrix, 2, 2, 1);
	matrixSetElem(universalGates[6].matrix, 2, 3, 0);

	matrixSetElem(universalGates[6].matrix, 3, 0, 0);
	matrixSetElem(universalGates[6].matrix, 3, 1, 1);
	matrixSetElem(universalGates[6].matrix, 3, 2, 0);
	matrixSetElem(universalGates[6].matrix, 3, 3, 0);


	// ================================== Generate Columns =================================

	nC = (int) pow (r,k); // Total number of columns

	// Stores in 'allColumns_ant' all 2x2 matrices (Universal Gates) and their names
	// This is necessary for the first iteration of all possible columns
	r = 7;

	struct columns *allColums_ant = malloc(r*sizeof(struct columns));
	for (i=0;i<r;i++){
		strcpy (allColums_ant[i].sequenceOfMatrices, universalGates[i].name);
		if (matrixGetNumOfColumns(universalGates[i].matrix)==2) {
			allColums_ant[i].matrix = matrixMult(universalGates[i].matrix,universalGates[0].matrix);
		} else {
			allColums_ant[i].matrix = matrixMult(universalGates[i].matrix,Iden);
		}
	}

 	// Prints all 2x2 matrices (Universal Set of size = r)

	for (l=0;l<r;l++){
		printf ("\n\nMatriz %s\n" , allColums_aux[l].sequenceOfMatrices);
		for (i=0;i<2;i++){
			for (j=0;j<2;j++){
				result = matrixGetElem(allColums_aux[l].matrix, i, j);
				printf("%f + %fi ", creal(result), cimag(result));
			}
			printf ("\n");
		}
	}



	//printf ("Distancia = %f" , matrixHilbertSchmidtDistance(universalGates[0].matrix, universalGates[1].matrix));

	struct columns *allColums;
	// Compute all possible (r^k) columns and stores them in allColumns
	char nome[1000];
	r = 7;
	base = r;
	for (i=1;i<=k-1;i++){
		//printf ("Base = %d e r = %d\n" , base, r);
		//nC = base*r;
		allColums = malloc(base*r*sizeof(struct columns));
		//printf ("Aloquei %d colunas\n" , nC);
		for (j=0;j<base;j++){
			for (l=0;l<r;l++){
				//printf ("antes j = %d e l = %d ---> posicao = %d\n", j , l , j*r+l);
				allColums[j*r+l].matrix = matrixKronProd(allColums_ant[j].matrix,universalGates[l].matrix);
				//printf("depois\n\n");
				strcpy (nome, allColums_ant[j].sequenceOfMatrices);
				strcat (nome, ",");
				strcat (nome , universalGates[l].name);
				strcpy (allColums[j*r+l].sequenceOfMatrices , nome);
			}
		}

		// Releases the memory used in previous iteration
		for (l=0;l<base;l++){
			matrixFree(allColums_ant[l].matrix);
		}
		free (allColums_ant);

		//Prepares for the next iteration
		allColums_ant = allColums;
		base = base*r;
		//r = r + sizesOfUnivGatesSet[i];

	}

	printf ("Na saida: Nc = %d, i = %d, base = %d e r = %d\n" , nC, i , base, r);

	//nC = 25;

	extraSize = i - 1;

	nC = nC + sizesOfUnivGatesSet[extraSize];

	allColums = realloc(allColums, nC*sizeof(struct columns));

	for (i = 0; i < sizesOfUnivGatesSet[extraSize]; i++){

		//printf ("Posição %d\n\n",nC+i-2);

		allColums[nC+i-2].matrix = universalGates[sizesOfUnivGatesSet[0]+i].matrix;
		strcpy (allColums[nC+i-2].sequenceOfMatrices , universalGates[sizesOfUnivGatesSet[0]+i].name);
	}

	printf ("Na saida: Nc = %d\n" , nC);



	// Prints all columns and their names
	printf ("%d colunas\n", nC);
	for (l=0;l<nC;l++){
		printf ("\n\nMatriz %s   %d x %d\n\n" , allColums[l].sequenceOfMatrices , matrixGetNumOfLines(allColums[l].matrix) , matrixGetNumOfColumns (allColums[l].matrix));
		for (i=0;i<N;i++){
			for (j=0;j<N;j++){
				result = matrixGetElem(allColums[l].matrix, i, j);
				printf("%.2f + %.2fi ; ", creal(result), cimag(result));
			}
			printf ("\n");
		}
	}

	printf("\n\n");

	exit (0);

	// ================================== Creates the Data Structure with all possibilities ================================

	int * combinations;
	combinations = malloc (maxIter*sizeof (int));
	for (i=0;i<maxIter;i++){
		combinations[i] = 0;
	}

	// ================================== Starts Looking for Good Solutions =================================

	printf ("Starting looking for good solutions ...\n\n");

	for (l=1;l<=pow (nC,maxIter);l++){

		//solutionAnt = matrixCreate(N,N);
		solutionAnt = matrixMult(Iden,Iden);
		//solution = matrixCreate(N,N);

		//printf ("Tentativa número %d\n", l);

		for (i=0;i<maxIter;i++){
			solution = matrixMult(solutionAnt , allColums[combinations[i]].matrix);
			matrixFree(solutionAnt);
			solutionAnt = solution;

			hilSchDist = matrixHilbertSchmidtDistance(U , solution);

			if (hilSchDist <= prec){
				//maxIter = i;
				//printf ("Achei com i = %d ", i);
				printf ("\n\nGood Solution with %lu columns:", i+1);
				for (j=0;j<=i;j++){
					if (j != i){
						printf ("%s x " , allColums[combinations[j]].sequenceOfMatrices);
					} else {
						printf ("%s" , allColums[combinations[j]].sequenceOfMatrices);
					}
				}

			//	printf ("\n");
			//	for (j=0;j<=i;j++){
			//		printf ("%d " , combinations[j]);
			//	}
			//	printf ("\n");

			//	int w;
			//	for (w=0;w<N;w++){
			//		for (j=0;j<N;j++){
			//			result = matrixGetElem(solution, w, j);
			//			printf("%.2f + %.2fi ; ", creal(result), cimag(result));
			//		}
			//		printf ("\n");
			//	}

			}
		}

		matrixFree(solution);

		//printf ("\n\n");

		//for (j=0;j<maxIter;j++){
			//printf ("%d " , combinations[j]);
		//}
		//printf ("\n");


		for (j=(maxIter-1);j>=0;j--){
			if (combinations[j] == (nC-1)){
				combinations[j] = 0;
			} else {
				combinations[j] = combinations[j] + 1;
				j = -1;
			}
		}
	}

	printf("\n\nSolutions found ...\n\n");


	// ================================== Finishes the Program =================================

	// Releases memory allocated to Universal Gate matrices
	for (i = 0; i < r; i++){
		matrixFree(universalGates[i].matrix);
	}

	printf ("Desalocou matrizes universais ... Nc = %d\n\n", nC);

	// Releases memory allocated for all columns

	for (i=0; i<nC; i++){
		matrixFree(allColums[i].matrix);
	}
	free (allColums);

	printf ("Desalocou colunas\n\n");

	printf ("Program finished normally\n");

	return 0;
}

*/


