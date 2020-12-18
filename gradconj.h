#ifndef GRADCONJ_H
#define GRADCONJ_H
#include"sparse.h"

int GradConj (int n, Linha * matrizEsparsa, double* b, double* x, double tol);
int GradConjPreCond(int n, Linha * matrizEsparsa, double* b, double* x, double tol, double w);
#endif