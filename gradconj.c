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
        //multmv(n,n,A,d,Adk);
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
int GradConjPreCondGaussSeidel(int n, Linha * matrizEsparsa, double* b, double* x, double tol) {
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
    Linha * M = criaPreCond(n,matrizEsparsa);
    printf("Gerou o Precondicionador\n");
    multiplicaMatrizEsparsaPorVetor(n,matrizEsparsa,x,Ax);
    printf("Multiplicou Matriz Esparsa por Vetor\n");

    subtraivet(n,b,Ax,rk0);
    //imprimeMatrizEsparsa(n,M);  
    GradConj(n,M,rk0,zk0,tol);
    clonavet(n,zk0,dk0);
    for(k = 0; k< n; k++) {
        normadois = norma2(n,rk0);
        if(normadois < tol) {
            return k;
        }
        //multmv(n,n,A,d,Adk);
        multiplicaMatrizEsparsaPorVetor(n,matrizEsparsa,dk0,Adk);
        rz = prodescalar(n,rk0,zk0);
        alpha =  rz / prodescalar(n,dk0,Adk); 
        multvs(n,dk0,alpha,alphaDk);
        multvs(n,Adk,alpha,alphaAdk);
        somavet(n,x,alphaDk,xn);
        subtraivet(n,rk0,alphaAdk,rk1);
        //multiplicaPreCondPorVetor(n,matrizEsparsa, 1.0, rk1, zk1);
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
// int GradConjPreCond (int n, double** A, double* b, double* x, double tol) {
//     int k; 
//     double alpha;
//     double beta;
//     double rz;
//     double normadois;
//     double* betaDk = criavet(n);
//     double* alphaDk = criavet(n);
//     double* alphaAdk = criavet(n);
//     double* Adk = criavet(n);
//     double* Ax = criavet(n);
//     double* r = criavet(n);
//     double* d = criavet(n);
//     double* z = criavet(n);
//     multiplicaMatrizEsparsaPorVetor(n,4,matrizEsparsa,x,Ax);
//     //multmv(n,n,A,x,Ax);
//     subtraivet(n,b,Ax,r);
//     multmv(n,n,Minversa,r,d);
//     z = d;
//     for(k = 0; k< n; k++) {
//         normadois = norma2(n,r);
//         if(normadois < tol) {
//             return k;
//         }
//         multmv(n,n,A,d,Adk);
//         rz = prodescalar(n,r,z);
//         alpha =  rz / (1+prodescalar(n,d,Adk)); 
//         multvs(n,d,alpha,alphaDk);
//         multvs(n,Adk,alpha,alphaAdk);
//         somavet(n,x,alphaDk,x);
//         subtraivet(n,r,alphaAdk,r);
//         multmv(n,n,Minversa,r,z);
//         beta = prodescalar(n,r,z) / (1+rz);
//         multvs(n,d,beta,betaDk);
//         somavet(n,z,betaDk,d);

//     } 
//    return k;
//}

