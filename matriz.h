#ifndef MATRIZ_H
#define MATRIZ_H

// cria (aloca) um vetor de dimensão n, retornando seu ponteiro
double* criavet (int n);

// libera (a memória) de um vetor previamente criado
void liberavet (double* v);

// calcula e retorna o produto escalar entre dois vetores de dimensão n
double prodescalar (int n, double* v, double* w);

// calcula e retorna a norma-2 de um vetor de dimensão n
double norma2 (int n, double* v);

// calcula a produto de um vetor v por um escalar s;
// o resultado deve ser armazenado no vetor w, previamente criado
void multvs (int n, double* v, double s, double *w);


// cria (aloca) uma matriz de dimensão m x n, retornando seu ponteiro;
// a matriz é representado por vetor de vetores linha
double** criamat (int m, int n);

// cria (aloca) uma matriz triangular inferior de dimensão n x n;
// a matriz é representado por vetor de vetores linha:
// o primeiro vetor linha tem dimensão 1, o segundo 2, e assim por diante
double** criamattri (int n);

// libera (a memória) de uma matriz previamente criada
void liberamat (int m, double** A);

// preenche a matriz transposta de A em T, previamente criada;
// A tem dimensão m x n; T tem dimensão n x m
void transposta (int m, int n, double** A, double** T);

// calcula o produto de uma matriz A (m x n) por um vetor v (m),
// resultando no vetor w (m), previamente criado
void multmv (int m, int n, double** A, double* v, double* w);

// calcula o produto de uma matriz A (m x n) por uma matriz B (n x q),
// armazenando o resultado na matriz C (m x q), previamente criada
void multmm (int m, int n, int q, double** A, double** B, double** C);

#endif
