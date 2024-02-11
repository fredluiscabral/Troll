/*
 * matriz.h
 *
 *  Created on: 4 de out de 2020
 *      Author: frederico
 */

#ifndef MATRIZ_H_
#define MATRIZ_H_

typedef struct matrix Matriz;

Matriz* matrixKronProd(Matriz *A, Matriz *B);
Matriz* matrixMult (Matriz* A, Matriz* B);
Matriz* matrixCreate (int m, int n);
Matriz* matrixAdd (Matriz* A, Matriz* B);
Matriz* matrixSub (Matriz* A, Matriz* B);
Matriz* matrixAdjoint(Matriz *A);
void matrixFree (Matriz* mat);
complex matrixGetElem (Matriz* mat, int i, int j);
void matrixSetElem (Matriz* mat, int i, int j, complex v);
int matrixGetNumOfLines (Matriz* mat);
int matrixGetNumOfColumns (Matriz* mat);
complex matrixTrace (Matriz *A);
double matrixHilbertSchmidtDistance(Matriz *A, Matriz *B);

#endif /* MATRIZ_H_ */
