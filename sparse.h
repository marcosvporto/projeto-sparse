#ifndef SPARSE_H
#define SPARSE_H

typedef struct elemento{
    double valor;
    int coluna;

} Elemento;

typedef struct linha {
    Elemento * elementos;
    int numeroDeElementos;
} Linha;

Linha * criaMatrizEsparsa(int n);

void imprimeMatrizEsparsa(int n, Linha * matrizEsparsa);
void insereElementoNaMatrizEsparsa(int linha, int coluna, double valor, Linha * matrizEsparsa);
void multiplicaMatrizEsparsaPorVetor(int n, Linha * matrizEsparsa, double *v, double * w);
double buscaBinariaElementoLinha (int coluna, int n, Elemento * elementos);
//void multiplicaPreCondPorVetor(int n,Linha * matrizEsparsa,double w, double *v, double* u);
Linha * criaPreCond(int n,Linha * matrizEsparsa);
#endif