#include "armaMex.hpp"
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <armadillo>

using namespace std;
using namespace arma;

float **A;
int Numb, Bsize;

/*Protótipos das funções */

void allocateA(mat)
mat BRI_fwdA(int, int*, int*);
mat BRI_fwdB(int, int*, int*);
mat BRI_fwdC(int, int*, int*);
mat BRI_fwdD(int, int*, int*);
float get_element(int, int);
mat gM (int, int);

/****************************************** FUNÇÃO main ********************************************/
mat init(mat X, int num_blocos, int tam_bloco) {

	Bsize = tam_bloco;
	Numb = num_blocos;
	allocateA(X);
	int * Brows = new int[Numb];
	int * Bcols = new int[Numb];
    
	for (int i = 0; i < Numb; i++) {
		Brows[i] = i;
	}
	for (int i = 0; i < Numb; i++) {
		Bcols[i] = i;
	}

    mat biM = pinv(BRI_fwdA(Numb,Brows,Bcols));

	biM.print("Bloco Final Invertido");

	delete[] Brows;
	delete[] Bcols;

	free (A);
	A = NULL;
    
    return biM;
}

/****************************************** FUNÇÃO allocateA ********************************************/
void allocateA(mat X) {
    int A_dimensao = Bsize*Numb;
	A = (float**)malloc(sizeof(float*)*A_dimensao);
    
    for(int i=0;i<A_dimensao;i++){
    	A[i] = (float*)malloc(sizeof(float)*A_dimensao);
        for(int j=0;j<A_dimensao;j++){
            A[i][j] = X(j+i*A_dimensao);
        }
    }
 
}

/****************************************** FUNÇÃO get_element ********************************************/

float get_element(int i, int j) {
	return A[i][j];
}

/********************************************* FUNÇÃO gM ********************************************/

mat gM (int row, int col){

	mat bloco (Bsize, Bsize);

	for (int i=0; i < Bsize; i++){
		for (int j=0; j < Bsize; j++){
			int A_i = row*Bsize + i;
			int A_j = col*Bsize + j;
			bloco(i,j) = get_element(A_i,A_j);
		}
	}
	return bloco;
}


/*************************************** FUNÇÃO BRI_fwdA ********************************************/
mat BRI_fwdA(int _numb, int* browsA, int* bcolsA){

    if (_numb>2){

   //     EXPRESSÃO bwdA = A - B*pinv(D)*C

/*************** TRABALHANDO COM B **********************/
	bcolsA+=1;

	mat bwdA = BRI_fwdB(_numb-1,browsA,bcolsA);

	bcolsA-=1;

/*************** TRABALHANDO COM D **********************/
	browsA+=1;
	bcolsA+=1;

    bwdA = bwdA * pinv(BRI_fwdD(_numb-1,browsA,bcolsA));

	browsA-=1;
	bcolsA-=1;

/*************** TRABALHANDO COM C **********************/
	browsA+=1;

    bwdA = bwdA * BRI_fwdC(_numb-1,browsA,bcolsA);

	browsA-=1;

/*************** TRABALHANDO COM A **********************/

   bwdA = BRI_fwdA(_numb-1,browsA,bcolsA) - bwdA;

	/*************** RETORNANDO bwdA **********************/
	return bwdA;
    }
    else{         
        mat bwdA = gM(browsA[0],bcolsA[1]);
        bwdA = bwdA * pinv(gM(browsA[1],bcolsA[1]));
        bwdA = bwdA * gM(browsA[1],bcolsA[0]);
        bwdA = gM(browsA[0],bcolsA[0]) - bwdA;

    	return bwdA;
     }
}  

/*************************************** FUNÇÃO BRI_fwdB ********************************************/
mat BRI_fwdB(int _numb, int* browsB, int* bcolsB){

    if (_numb>2){
	/*************** INVERTENDO bcolsB **********************/
	int temp = bcolsB[1];
	bcolsB[1] = bcolsB[0];
	bcolsB[0] = temp;

   //     EXPRESSÃO  bwdB = B - A*pinv(C)*D;

/*************** TRABALHANDO COM A **********************/

        mat bwdB = BRI_fwdA(_numb-1,browsB,bcolsB);

/*************** TRABALHANDO COM C **********************/
	browsB+=1;

    bwdB = bwdB * pinv(BRI_fwdC(_numb-1,browsB,bcolsB));
        
	browsB-=1;

/*************** TRABALHANDO COM D **********************/
	browsB+=1;
	bcolsB+=1;

    bwdB = bwdB * BRI_fwdD(_numb-1,browsB,bcolsB);

	browsB-=1;
	bcolsB-=1;

/*************** TRABALHANDO COM B **********************/
	bcolsB+=1;

    bwdB = BRI_fwdB(_numb-1,browsB,bcolsB) - bwdB;

	bcolsB-=1;

	/*************** RETORNANDO bwdB **********************/
	return bwdB;
    }
   else{  
        mat bwdB = gM(browsB[0],bcolsB[0]);
    	bwdB = bwdB * pinv(gM(browsB[1],bcolsB[0]));
        bwdB = bwdB * gM(browsB[1],bcolsB[1]);
        bwdB = gM(browsB[0],bcolsB[1]) - bwdB;

        return bwdB;
     }
} 

/*************************************** FUNÇÃO BRI_fwdC ********************************************/
mat BRI_fwdC(int _numb, int* browsC, int* bcolsC){

    if (_numb>2){
	/*************** INVERTENDO browsC **********************/
	int temp = browsC[1];
	browsC[1] = browsC[0];
	browsC[0] = temp;

   //     EXPRESSÃO  bwdC = C - D*pinv(B)*A;

/*************** TRABALHANDO COM D **********************/
	browsC+=1;
	bcolsC+=1;

    mat bwdC = BRI_fwdD(_numb-1,browsC,bcolsC);

	browsC-=1;
	bcolsC-=1;

/*************** TRABALHANDO COM B **********************/
	bcolsC+=1;

    bwdC = bwdC * pinv(BRI_fwdB(_numb-1,browsC,bcolsC));

	bcolsC-=1;

/*************** TRABALHANDO COM A **********************/

    bwdC = bwdC * BRI_fwdA(_numb-1,browsC,bcolsC);


/*************** TRABALHANDO COM C **********************/
	browsC+=1;

    bwdC = BRI_fwdC(_numb-1,browsC,bcolsC) - bwdC;

	browsC-=1;

	/*************** RETORNANDO bwdC **********************/

	return bwdC;
    }
    else{         
        mat bwdC = gM(browsC[1],bcolsC[1]);
        bwdC = bwdC * pinv(gM(browsC[0],bcolsC[1]));
        bwdC = bwdC * gM(browsC[0],bcolsC[0]);
        bwdC = gM(browsC[1],bcolsC[0]) - bwdC;

        return bwdC;
     }
} 

/*************************************** FUNÇÃO BRI_fwdD ********************************************/
mat BRI_fwdD(int _numb, int* browsD, int* bcolsD){

    if (_numb>2){
	/*************** INVERTENDO bcolsD e browsD **********************/
	int temp = browsD[1];
	browsD[1] = browsD[0];
	browsD[0] = temp;

	temp = bcolsD[1];
	bcolsD[1] = bcolsD[0];
	bcolsD[0] = temp;

   //     EXPRESSÃO   bwdD = D - C*pinv(A)*B

/*************** TRABALHANDO COM C **********************/
	browsD+=1;

    mat bwdD = BRI_fwdC(_numb-1,browsD,bcolsD);

	browsD-=1;


/*************** TRABALHANDO COM A **********************/

        bwdD = bwdD * pinv(BRI_fwdA(_numb-1,browsD,bcolsD));

/*************** TRABALHANDO COM B **********************/
	bcolsD+=1;

    bwdD = bwdD * BRI_fwdB(_numb-1,browsD,bcolsD);
        
	bcolsD-=1;

/*************** TRABALHANDO COM D **********************/
	browsD+=1;
	bcolsD+=1;

    bwdD = BRI_fwdD(_numb-1,browsD,bcolsD) - bwdD;

	browsD-=1;
	bcolsD-=1;

	/*************** RETORNANDO bwdD **********************/
	return bwdD;
    }
    else{         
        mat bwdD = gM(browsD[1],bcolsD[0]);
        bwdD = bwdD * pinv(gM(browsD[0],bcolsD[0]));
        bwdD = bwdD * gM(browsD[0],bcolsD[1]);
        bwdD = gM(browsD[1],bcolsD[1]) - bwdD;

        return bwdD;
     }
} 

void mexFunction(int nlhs, mxArray *plhs[],
         int nrhs, const mxArray *prhs[]) 
{   
	//BRI_in(1): numero de blocos BRI_in(2): dimensao dos blocos
    mat X = conv_to<mat>::from(armaGetPr(prhs[0], true));
    mat BRI_in = conv_to<mat>::from(armaGetPr(prhs[1], true));
    plhs[0] = armaCreateMxMatrix(BRI_in(1), BRI_in(1), mxDOUBLE_CLASS, mxREAL);
    armaSetPr(plhs[0], conv_to<mat>::from(init(X,BRI_in(0),BRI_in(1))));
}
