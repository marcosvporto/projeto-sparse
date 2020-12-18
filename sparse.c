#include"sparse.h"
#include<stdlib.h>
#include<stdio.h>
#include"gradconj.h"


double buscaBinariaElementoLinha (int coluna, int n, Elemento * elementos) {
    int inicio = 0, fim = n-1;
    while ( inicio <= fim ){
        int meio = (inicio + fim)/2;
        if (elementos[meio].coluna == coluna) {
            return elementos[meio].valor;
        } else if ( coluna < elementos[meio].coluna) {
            fim = meio -1;
        } else if (coluna > elementos[meio].coluna) {
            inicio = meio +1;
        }
    }
    return 0.0;
}

Linha * criaMatrizEsparsa(int n){
    Linha * matrizEsparsa = (Linha *) malloc(n*sizeof(Linha));
    if (matrizEsparsa == NULL) {
        printf("Erro ao criar matriz esparsa!!\n");
    }
    return matrizEsparsa;
}

void insereElementoNaMatrizEsparsa(int linha, int coluna, double valor, Linha * matrizEsparsa) {
    Elemento * aux;
    
    int elementosLinha = matrizEsparsa[linha].numeroDeElementos;
    if( elementosLinha == 0) {
        
        matrizEsparsa[linha].elementos = (Elemento*)malloc(sizeof(Elemento));
        matrizEsparsa[linha].elementos[0].valor = valor;
        matrizEsparsa[linha].elementos[0].coluna = coluna;
        matrizEsparsa[linha].numeroDeElementos = 1;
    } else {
        
        elementosLinha = matrizEsparsa[linha].numeroDeElementos;
        aux = matrizEsparsa[linha].elementos;
        matrizEsparsa[linha].elementos = (Elemento*)malloc((1+elementosLinha)*sizeof(Elemento));
        for (int i =0; i<elementosLinha; i++){
            matrizEsparsa[linha].elementos[i] = aux[i];
        }
        matrizEsparsa[linha].elementos[elementosLinha].valor = valor;
        matrizEsparsa[linha].elementos[elementosLinha].coluna = coluna;
        free(aux);
        matrizEsparsa[linha].numeroDeElementos = 1+elementosLinha;
    }
    // printf("Numero de elementos %d\n",matrizEsparsa[linha].numeroDeElementos);
    // printf("Inserido valor %.2f\n", matrizEsparsa[linha].elementos[0].valor);
    // printf("Na coluna %d\n", matrizEsparsa[linha].elementos[0].coluna);
    // printf("Na linha %d\n", linha);
}

void imprimeMatrizEsparsa(int n, Linha * matrizEsparsa){

    for (int i = 0; i < n; i++) {
        printf("|");
        int elementoPorLinha = matrizEsparsa[i].numeroDeElementos;
        for (int j = 0; j < n; j++){
            double valor = buscaBinariaElementoLinha(j, elementoPorLinha, matrizEsparsa[i].elementos);
            if (valor == 0.0 && j!=i) {
                elementoPorLinha = matrizEsparsa[j].numeroDeElementos;
                valor = buscaBinariaElementoLinha(i, elementoPorLinha,matrizEsparsa[j].elementos);
            }
            printf(" %.2f ", valor );
        }
        printf("|\n");
    }
}

void imprimeMatrizEsparsaNaoSimetrica(int n, Linha * matrizEsparsa){

    for (int i = 0; i < n; i++) {
        printf("|");
        int elementoPorLinha = matrizEsparsa[i].numeroDeElementos;
        for (int j = 0; j < n; j++){
            double valor = buscaBinariaElementoLinha(j, elementoPorLinha, matrizEsparsa[i].elementos);
            printf(" %.2f ", valor );
        }
        printf("|\n");
    }
}

Linha * criaPreCond(int n,Linha * matrizEsparsa){
    int elementoPorLinhaI;
    int elementoPorLinhaJ;
    int elementoPorLinhaDwLDinv;
    int elementoPorLinhaDU;
    double valor = 0;
    Linha * matrizDwLDinv = criaMatrizEsparsa(n);
    Linha * M = criaMatrizEsparsa(n);
    for(int i = 0; i<n;i++){
        for(int k=0;k<n;k++){
            valor = 0;
            for(int j = 0; j< n; j++){
                if (j==k){
                    double Aij;
                    double Bjk;
                    elementoPorLinhaI = matrizEsparsa[i].numeroDeElementos;
                    elementoPorLinhaJ = matrizEsparsa[j].numeroDeElementos;
                    Aij = buscaBinariaElementoLinha(i,elementoPorLinhaJ,matrizEsparsa[j].elementos);
                    if(Aij != 0) {
                        Bjk = buscaBinariaElementoLinha(j,elementoPorLinhaI,matrizEsparsa[j].elementos);
                        if(Bjk != 0) {
                            valor += Aij * (1.0/Bjk);
                        }
                    }
                }
            }
            if(valor != 0){
                insereElementoNaMatrizEsparsa(i,k,valor,matrizDwLDinv);
            }
        }
    }
    printf("Gerou o D+wLDinv\n");
    //imprimeMatrizEsparsaNaoSimetrica(n,matrizDwLDinv);
    // printf("*****\n");
    // imprimeMatrizEsparsaNaoSimetrica(n,matrizEsparsa);
    // printf("*****\n");

    for(int i = 0; i<n;i++){
        for(int k=0;k<n;k++){
            valor = 0;
            elementoPorLinhaDwLDinv = matrizDwLDinv[i].numeroDeElementos;
            double Aij;
            double Bjk;
            for(int y = 0; y < elementoPorLinhaDwLDinv; y++){
                int j = matrizDwLDinv[i].elementos[y].coluna;
                elementoPorLinhaDU = matrizEsparsa[j].numeroDeElementos;
                Aij = matrizDwLDinv[i].elementos[y].valor;
                Bjk = buscaBinariaElementoLinha(k,elementoPorLinhaDU,matrizEsparsa[j].elementos);
                valor += Aij*Bjk;
            }
            
            if(valor != 0){
                insereElementoNaMatrizEsparsa(i,k,valor,M);
            }
        }
    }
    return M;
}

void multiplicaMatrizEsparsaPorVetor(int n, Linha * matrizEsparsa, double *v, double * w) {
    int elementoPorLinha;
    for(int i = 0; i < n; i++){
        w[i] = 0.0;
    }
    for(int i = 0; i < n; i++){
        elementoPorLinha = matrizEsparsa[i].numeroDeElementos;
        for(int j = 0; j < elementoPorLinha; j++){
            int coluna = matrizEsparsa[i].elementos[j].coluna;
            double valor = matrizEsparsa[i].elementos[j].valor; // A[i][coluna] = valor
            w[i]      += valor * v[coluna]; // w[i] = w[i] + (A[i][coluna]*v[coluna])
            if(i != coluna){
                w[coluna] += valor * v[i];
            }
        }
    }
}