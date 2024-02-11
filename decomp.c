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

struct columns {
	char sequenceOfMatrices[1000];
	Matriz * matrix;
};



struct table {
	char sequenceOfColumns[1000];
	Matriz * matrix;
};




int main() {

	int i,j,l;
	double N = 4; // U port dimension
	double k = (int) (log10(N) / log10(2)); // Column size
	int r = 5; // Size of the Universal Set os Gates
	int nC; // Number of columns
	struct matrixName *universalGates = malloc((r+2)*sizeof(struct matrixName));
	complex result;
	int base;
	Matriz *U, *Iden; //Target Matrix U and identity N x N
	int maxIter = 6;
	double hilSchDist, prec = 0.01;


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

/*
	// Fourier Matrix 4x4
	U = matrixCreate(4,4);
	matrixSetElem(U,0,0,sqrt(2)/2);
	matrixSetElem(U,0,1,sqrt(2)/2);
	matrixSetElem(U,0,2,0);
	matrixSetElem(U,0,3,0);

	matrixSetElem(U,1,0,(1 + I)/2);
	matrixSetElem(U,1,1,-(1 + I)/2);
	matrixSetElem(U,1,2,0);
	matrixSetElem(U,1,3,0);

	matrixSetElem(U,2,0,0);
	matrixSetElem(U,2,1,0);
	matrixSetElem(U,2,2,(1 + I)/2);
	matrixSetElem(U,2,3,(1 + I)/2);

	matrixSetElem(U,3,0,0);
	matrixSetElem(U,3,1,0);
	matrixSetElem(U,3,2,(I*sqrt(2))/2);
	matrixSetElem(U,3,3,-(I*sqrt(2))/2);
*/

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
	strcpy (universalGates[4].name,"S");
	matrixSetElem(universalGates[4].matrix, 0, 0, 1);
	matrixSetElem(universalGates[4].matrix, 0, 1, 0);
	matrixSetElem(universalGates[4].matrix, 1, 0, 0);
	matrixSetElem(universalGates[4].matrix, 1, 1, I);


	// Define extra gates (4x4 columns) CNOT12 and CNOT21

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
	struct columns *allColums_ant = malloc(r*sizeof(struct columns));
	for (i=0;i<r;i++){
		strcpy (allColums_ant[i].sequenceOfMatrices, universalGates[i].name);
		allColums_ant[i].matrix = matrixMult(universalGates[i].matrix,universalGates[0].matrix);
	}


 	// Prints all 2x2 matrices (Universal Set of size = r)
	for (l=0;l<r;l++){
		printf ("\n\nMatriz %s\n" , allColums_ant[l].sequenceOfMatrices);
		for (i=0;i<2;i++){
			for (j=0;j<2;j++){
				result = matrixGetElem(allColums_ant[l].matrix, i, j);
				printf("%f + %fi ", creal(result), cimag(result));
			}
			printf ("\n");
		}
	}




	//printf ("Distancia = %f" , matrixHilbertSchmidtDistance(universalGates[0].matrix, universalGates[1].matrix));

	struct columns *allColums;
	// Compute all possible (r^k) columns and stores them in allColumns
	char nome[1000];
	base = r;
	for (i=1;i<=k-1;i++){
		allColums = malloc(base*r*sizeof(struct columns));
		for (j=0;j<base;j++){
			for (l=0;l<r;l++){
				allColums[j*r+l].matrix = matrixKronProd(allColums_ant[j].matrix,universalGates[l].matrix);
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

	}

	allColums = realloc(allColums, (nC + 2)*sizeof(struct columns));

	allColums[nC].matrix = universalGates[r].matrix;
	strcpy (allColums[nC].sequenceOfMatrices , universalGates[r].name);

	allColums[nC+1].matrix = universalGates[r+1].matrix;
	strcpy (allColums[nC+1].sequenceOfMatrices , universalGates[r+1].name);



	nC = nC + 2;

	// Prints all (nC = r^k) columns and their names
	printf ("\n\n%d colunas\n", nC);
	for (l=0;l<nC;l++){
		printf ("\n\nMatriz %s   %d x %d\n" , allColums[l].sequenceOfMatrices , matrixGetNumOfLines(allColums[l].matrix) , matrixGetNumOfColumns (allColums[l].matrix));
		for (i=0;i<N;i++){
			for (j=0;j<N;j++){
				result = matrixGetElem(allColums[l].matrix, i, j);
				printf(" %f + %fi ; ", creal(result), cimag(result));
			}
			printf ("\n");
		}
	}


	// ================================== Starts Looking for Good Solutions =================================

	// Prepares for the first iteration



	struct table *allMatrices;
	struct table *allMatrices_ant = malloc(nC*sizeof(struct table));
	for (i=0;i<nC;i++){
		strcpy (allMatrices_ant[i].sequenceOfColumns, allColums[i].sequenceOfMatrices);
		allMatrices_ant[i].matrix = matrixMult(allColums[i].matrix , Iden);
	}

	i = 1;
	base = nC;
	int foundSolution = 0;
	int shouldStop = 0;
	do {

		printf ("Looking for solutions with %d columns ...\n",i);

		allMatrices = malloc(nC*base*sizeof(struct table));

		for (j=0;j<base;j++){
			for (l=0;l<nC;l++){
				allMatrices[j*nC+l].matrix = matrixMult(allMatrices_ant[j].matrix,allColums[l].matrix);
				strcpy (nome, allMatrices_ant[j].sequenceOfColumns);
				strcat (nome, "x");
				strcat (nome , allColums[l].sequenceOfMatrices);
				strcpy (allMatrices[j*nC+l].sequenceOfColumns , nome);
			}
		}

		for (j=0;j<base;j++){
			hilSchDist = matrixHilbertSchmidtDistance(U , allMatrices_ant[j].matrix);
			if (hilSchDist <= prec){
				printf ("Good Solution %s\n\n" , allMatrices_ant[j].sequenceOfColumns);
				foundSolution = 1;
			}
		}

		if (i == maxIter){
			shouldStop = 1;
		}

		i++;

		// Releases memory used at previous iteration
		for (l=0;l<base;l++){
			matrixFree(allMatrices_ant[l].matrix);
		}
		free (allMatrices_ant);

		allMatrices_ant = allMatrices;
		base = base*nC;
	} while (!foundSolution && !shouldStop);

	if (!foundSolution){
		printf ("Could not find a good solution with %d columns or less\n" , maxIter);
	}

	// ================================== Finishes the Program =================================

	// Releases memory allocated to Universal Gate matrices
	for (i = 0; i < r; i++){
		matrixFree(universalGates[i].matrix);
	}

	// Releases memory allocated for all columns 
	for (i = 0; i < nC; i++){
		matrixFree(allColums[i].matrix);
	}
	free (allColums);

	// Releases memory allocated for all matrices
	for (i = 0; i < nC; i++){
		matrixFree(allMatrices[i].matrix);
	}
	free (allMatrices);

	printf ("Program finished normally\n");

	return 0;



}



