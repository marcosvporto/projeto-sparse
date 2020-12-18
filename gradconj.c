#include"gradconj.h"
#include"sparse.h"
#include"matriz.h"
#include<math.h>
#include<stdio.h>


void somavet(int n,double* v, double* w, double* u) {
    for(int i = 0; i < n ; i++){
        u[i] = v[i]+w[i];
    }
}
void subtraivet(int n,double* v, double* w, double* u) {
    for(int i = 0; i < n ; i++){
        u[i] = v[i] - w[i];
    }
}
void clonavet(int n, double *origem, double*destino) {
    for(int i = 0; i< n;i++){
        destino[i] = origem[i];
    }
}
int GradConj(int n, Linha * matrizEsparsa, double* b, double* x, double tol) {
    int k; 
    double alpha;
    double beta;
    double rr;
    double normadois;
    double* betaDk = criavet(n);
    double* alphaDk = criavet(n);
    double* alphaAdk = criavet(n);
    double* Adk = criavet(n);
    double* Ax = criavet(n);
    double* rk0 = criavet(n);
    double* rk1 = criavet(n);
    double* dk0 = criavet(n);
    double* dk1 = criavet(n);
    double* xn = criavet(n);
    
    multiplicaMatrizEsparsaPorVetor(n,matrizEsparsa,x,Ax);
    subtraivet(n,b,Ax,rk0);
    clonavet(n,rk0,dk0);

    for(k = 0; k< n; k++) {
        normadois = norma2(n,rk0);
        if(normadois < tol) {
            return k;
        }
        multiplicaMatrizEsparsaPorVetor(n,matrizEsparsa,dk0,Adk);
        rr = prodescalar(n,rk0,rk0);
        alpha =  rr / prodescalar(n,dk0,Adk); 
        multvs(n,dk0,alpha,alphaDk);
        multvs(n,Adk,alpha,alphaAdk);
        somavet(n,x,alphaDk,xn);
        subtraivet(n,rk0,alphaAdk,rk1);
        beta = prodescalar(n,rk1,rk1) / prodescalar(n,rk0,rk0);
        multvs(n,dk0,beta,betaDk);
        somavet(n,rk1,betaDk,dk1);
        clonavet(n,rk1,rk0);
        clonavet(n,dk1,dk0);
        clonavet(n,xn,x);
    }
    return k;
}
int GradConjPreCond(int n, Linha * matrizEsparsa, double* b, double* x, double tol, double w) {
    int k; 
    double alpha;
    double beta;
    double rz;
    double normadois;
    double* betaDk = criavet(n);
    double* alphaDk = criavet(n);
    double* alphaAdk = criavet(n);
    double* Adk = criavet(n);
    double* Ax = criavet(n);
    double* rk0 = criavet(n);
    double* rk1 = criavet(n);
    double* dk0 = criavet(n);
    double* dk1 = criavet(n);
    double* zk0 = criavet(n);
    double* zk1 = criavet(n);
    double* xn = criavet(n);
    printf("Iniciando geração do Precondicionador\n");
    Linha * M = criaPreCond(n,matrizEsparsa,w);
    printf("Gerou o Precondicionador\n");
    multiplicaMatrizEsparsaPorVetor(n,matrizEsparsa,x,Ax);
    printf("Multiplicou Matriz Esparsa por Vetor\n");

    subtraivet(n,b,Ax,rk0);
    GradConj(n,M,rk0,zk0,tol);
    clonavet(n,zk0,dk0);
    for(k = 0; k< n; k++) {
        normadois = norma2(n,rk0);
        if(normadois < tol) {
            return k;
        }
        multiplicaMatrizEsparsaPorVetor(n,matrizEsparsa,dk0,Adk);
        rz = prodescalar(n,rk0,zk0);
        alpha =  rz / prodescalar(n,dk0,Adk); 
        multvs(n,dk0,alpha,alphaDk);
        multvs(n,Adk,alpha,alphaAdk);
        somavet(n,x,alphaDk,xn);
        subtraivet(n,rk0,alphaAdk,rk1);
        GradConj(n,M,rk1,zk1,tol);
        beta = prodescalar(n,rk1,zk1) / prodescalar(n,rk0,zk0);
        multvs(n,dk0,beta,betaDk);
        somavet(n,zk1,betaDk,dk1);
        clonavet(n,rk1,rk0);
        clonavet(n,dk1,dk0);
        clonavet(n,zk1,zk0);
        clonavet(n,xn,x);
    }
    return k;
}

