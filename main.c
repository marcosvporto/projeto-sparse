#include<stdio.h>
#include<stdlib.h>
#include"sparse.h"
#include"gradconj.h"

#define N 1000
#define TOL 0.0000001
#define W 0
void imprimeVetor(double *v, int n){
    for(int i =0; i< n; i++) {
        printf("|%.6f|\n", *(v+i));
    }
    printf("\n");
}
int main(){
    double *x = (double*)malloc(N*sizeof(double));
    double *xe = (double*)malloc(N*sizeof(double));
    double *b = (double*)malloc(N*sizeof(double));
    Linha * matriz = criaMatrizEsparsa(N);

    
    for(int i = 0; i<N; i++) {
        xe[i] = 1.0;
        x[i] = 0.0;
        for(int j = 0; j<N; j++){
            if(j==i){
                insereElementoNaMatrizEsparsa(i,j,1.0*(i+1),matriz);
            } else if(j == i+1 || j == i+2 || j == i*2) {
                insereElementoNaMatrizEsparsa(i,j,0.5,matriz);
            }
        }
    }
    multiplicaMatrizEsparsaPorVetor(N,matriz,xe,b);
    int tentativas = GradConj(N,matriz,b,x,TOL);
    printf("%d tentativas\n", tentativas);
    multiplicaMatrizEsparsaPorVetor(N,matriz,xe,b);
    for(int i = 0; i < N ; i++){
        x[i] = 0.0;
    }
    tentativas = GradConjPreCond(N,matriz,b,x,TOL,W);
    printf("%d tentativas\n", tentativas);
    return 0;
}