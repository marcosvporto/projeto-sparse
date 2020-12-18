#include<stdio.h>
#include<stdlib.h>
#include"sparse.h"
#include"gradconj.h"

#define N 49
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
    //imprimeMatrizEsparsa(N, matriz);
    multiplicaMatrizEsparsaPorVetor(N,matriz,xe,b);
    int tentativas = GradConj(N,matriz,b,x,0.01);
    printf("%d tentativas\n", tentativas);
    for(int i = 0; i<N ; i++){
        x[i] = 0;
    }
    tentativas = GradConjPreCondGaussSeidel(N,matriz,b,x,0.01);
    printf("%d tentativas\n", tentativas);
    imprimeVetor(x,N);
    return 0;
}