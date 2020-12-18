#include"matriz.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// cria (aloca) um vetor de dimensão n, retornando seu ponteiro
double* criavet (int n){
    double *v = (double*) malloc(n*sizeof(double));
    return v;
}

// libera (a memória) de um vetor previamente criado
void liberavet (double* v){
    free(v);
}

// calcula e retorna o produto escalar entre dois vetores de dimensão n
double prodescalar (int n, double* v, double* w){
    double produtoEscalar = 0;
    for(int i = 0; i < n; i++) {
        produtoEscalar+= *(v+i) * *(w+i);
    }
    return produtoEscalar;
}

// calcula e retorna a norma-2 de um vetor de dimensão n
double norma2 (int n, double* v){
    double somaDosQuadrados = 0;
    for(int i = 0; i < n; i++) {
        somaDosQuadrados += ((*(v+i)) * (*(v+i)));
    }
    double norma = sqrt(somaDosQuadrados); 
    return norma;
}

// calcula a produto de um vetor v por um escalar s;
// o resultado deve ser armazenado no vetor w, previamente criado
void multvs (int n, double* v, double s, double *w){
    for(int i=0; i< n;i++){
        *(w+i) = *(v+i) * s;
    }
}


// cria (aloca) uma matriz de dimensão m x n, retornando seu ponteiro;
// a matriz é representado por vetor de vetores linha
double** criamat (int m, int n){
    double ** matriz = (double**) malloc(m*sizeof(double*));
    for(int i =0;i<m;i++){
        *(matriz+i) = (double *) malloc(n*sizeof(double));
    }
    return matriz;
}

// cria (aloca) uma matriz triangular inferior de dimensão n x n;
// a matriz é representado por vetor de vetores linha:
// o primeiro vetor linha tem dimensão 1, o segundo 2, e assim por diante
double** criamattri (int n){
    double ** matrizTriangular = (double**) malloc(n*sizeof(double*));
    for(int i =0;i<n;i++){
        *(matrizTriangular+i) = (double *) malloc((i+1)*sizeof(double));
    }
    return matrizTriangular;
}


// libera (a memória) de uma matriz previamente criada
void liberamat (int m, double** A){
    for(int i =0; i< m; i++) {
        free(*(A+i));
    }
}

// preenche a matriz transposta de A em T, previamente criada;
// A tem dimensão m x n; T tem dimensão n x m
void transposta (int m, int n, double** A, double** T){
    for (int i=0; i< m; i++){
        for(int j = 0; j<n ; j++ ){
            *(*(T+j)+i) = *(*(A+i)+j); 
        }
    }
}

// calcula o produto de uma matriz A (m x n) por um vetor v (m),
// resultando no vetor w (m), previamente criado
void multmv (int m, int n, double** A, double* v, double* w){
    for(int i=0;i<m;i++){
        *(w+i) = 0;
        for(int j=0; j< n; j++){
            *(w+i) += *(*(A+i)+j) * *(v+j);
        }
    } 
}

// calcula o produto de uma matriz A (m x n) por uma matriz B (n x q),
// armazenando o resultado na matriz C (m x q), previamente criada
void multmm (int m, int n, int q, double** A, double** B, double** C){
    // i -> m
    // j -> n
    // k -> q
    for (int i = 0; i < m ; i++) {
        for (int k = 0 ; k< q ; k++) {
            *(*(C+i)+k) = 0;
            for(int j = 0; j < n; j ++) {
                *(*(C+i)+k) += *(*(A+i)+j) * *(*(B+j)+k); 
            }
        }
    }
}