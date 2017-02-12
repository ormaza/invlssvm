/*Algoritmo BRI - versão 1.5: Valores da matrix de amostras X são aleatórios, mas buscam-se realizar testes mantendo o tamanho do bloco constante e variando apenas a qtde de blocos. Usa-se a função "randn" */
#include "armaMex.hpp"
#include <armadillo>

using namespace std;
using namespace arma;

double **X;
unsigned int Numb, Bsize, d;
float Sigma, Gamma;

/*Protótipos das funções */

void allocateX(unsigned int, unsigned int, mat);
mat BRI_fwdA(unsigned int, unsigned int*, unsigned int*, unsigned int*, unsigned int*);
mat BRI_fwdB(unsigned int, unsigned int*, unsigned int*, unsigned int*, unsigned int*);
mat BRI_fwdC(unsigned int, unsigned int*, unsigned int*, unsigned int*, unsigned int*);
mat BRI_fwdD(unsigned int, unsigned int*, unsigned int*, unsigned int*, unsigned int*);
double get_element(unsigned int, unsigned int);
mat gM (unsigned int, unsigned int, unsigned int*, unsigned int*);


/****************************************** FUNÇÃO main ********************************************/
mat init(mat dados, unsigned int block_amount, unsigned int block_size, unsigned int data_dim, float gam, float sig, unsigned int block_select) {

	unsigned int qx=0;
	
	Bsize = block_size;
	Numb = block_amount;
	Gamma = gam;
	Sigma = sig;
	
	qx = (Bsize*Numb)-1;
	//A_dimension = qx+1
	d=data_dim; //block dimension
	
	cout << "\nQtde de amostras (X) = " << qx << ' ' << endl;


	allocateX(qx,d,dados);

	// printf ("\n\n\nDisplaying Matrix X:\n\n");
	// for (unsigned int i=0;i<qx;i++) {
	// 	printf("\n");
	// 	unsigned int j;
	// 	for (j = 0; j < d; j++){ 
 //   	        	printf("%f \t", X[i][j]); 
 //        	} 
	// }
	// printf("\n\n");

	unsigned int * Brows = new unsigned int[Numb];
	unsigned int * Bcols = new unsigned int[Numb];

    for (unsigned int i = 0; i < Numb; i++) {
		Brows[i] = i;
	}
	for (unsigned int i = 0; i < Numb; i++) {
		Bcols[i] = i;
	}


	unsigned int * Brows_index = new unsigned int[Numb];
	unsigned int * Bcols_index = new unsigned int[Numb];

    for (unsigned int i = 0; i < Numb; i++) {
		Brows_index[i] = i;
	}
	for (unsigned int i = 0; i < Numb; i++) {
		Bcols_index[i] = i;
	}

	Brows_index[0] = block_select;
	Bcols_index[0] = block_select;
	Brows_index[block_select] = 0;
	Brows_index[block_select] = 0;

    
    mat biM = pinv(BRI_fwdA(Numb,Brows,Bcols,Brows_index,Bcols_index));

	// biM.print("Bloco Final pinvertido");

	delete[] Brows;
	delete[] Bcols;

	unsigned int i;	
	for (i=0;i<qx;i++) {
		free(X[i]);
	}
	free (X);
	X = NULL; 

	return biM;
}

/********************************************* FUNÇÃO allocateX ********************************************/
void allocateX(unsigned int x_quantity, unsigned int x_dimension, mat dados) {
    // x_dimension==d e x_quantity==sample quantity
    // Allocate memory for Matrix X

    X = (double**)malloc(sizeof(double*)*x_quantity);

    unsigned int i,j;
    for (i=0; i<x_quantity; i++) {
        X[i] = (double*)malloc(sizeof(double)*x_dimension);
        for (j = 0; j < x_dimension; j++){
		X[i][j] = dados(j+i*x_dimension);		
		}
    }

}

/********************************************* FUNÇÃO get_element ********************************************/
double get_element(unsigned int A_i, unsigned int A_j) {
	if (A_i==(Bsize*Numb)-1 && A_j==(Bsize*Numb)-1){
        return 0;
        }
        else if (A_i==(Bsize*Numb)-1 || A_j==(Bsize*Numb)-1) {
            return 1;
        }
	
	double elem = 0;
        unsigned int k;
	for (k=0;k<d;k++) {
        elem = elem - pow(X[A_i][k],2) - pow(X[A_j][k],2) + 2*X[A_i][k]*X[A_j][k];
	}
	elem = expf(elem/(2*(Sigma*Sigma))) + (A_i==A_j)/Gamma; 

	return elem;
}


/********************************************* FUNÇÃO gM ********************************************/

mat gM (unsigned int row, unsigned int col, unsigned int * Brows_index, unsigned int * Bcols_index){

	row = Brows_index[row];
	col = Bcols_index[col];

	mat bloco (Bsize, Bsize);

	for (unsigned int i=0; i < Bsize; i++){
		for (unsigned int j=0; j < Bsize; j++){
			unsigned int A_i = row*Bsize + i;
			unsigned int A_j = col*Bsize + j;
			bloco(i,j) = get_element(A_i,A_j);
		}
	}
	return bloco;
}


/*************************************** FUNÇÃO BRI_fwdA ********************************************/
mat BRI_fwdA(unsigned int _numb, unsigned int* browsA, unsigned int* bcolsA, unsigned int* Brows_index, unsigned int* Bcols_index){

    //printf ("\n\n\nEntrou na funcao:\n\n");
    if (_numb>2){

   //     EXPRESSÃO bwdA = A - B*pinv(D)*C

/*************** TRABALHANDO COM B **********************/
        //printf ("\n\n\nB - bcolsA passa a apontar para o elemento 1 do vetor: \n\n");
	bcolsA+=1;

	mat bwdA = BRI_fwdB(_numb-1,browsA,bcolsA,Brows_index,Bcols_index);

        //printf ("\n\n\nB - bcolsA volta a apontar para o elemento 0 do vetor: \n\n");
	bcolsA-=1;

/*************** TRABALHANDO COM D **********************/
	//printf ("\n\n\nD - browsA e bcolsA passam a apontar para o elemento 1 do vetor: \n\n");
	browsA+=1;
	bcolsA+=1;

        bwdA = bwdA * pinv(BRI_fwdD(_numb-1,browsA,bcolsA,Brows_index,Bcols_index));

	//printf ("\n\n\nD - browsA e bcolsA voltam a apontar para o elemento 0 do vetor: \n\n");
	browsA-=1;
	bcolsA-=1;

/*************** TRABALHANDO COM C **********************/
	//printf ("\n\n\nC - browsA passa a apontar para o elemento 1 do vetor: \n\n");
	browsA+=1;

        bwdA = bwdA * BRI_fwdC(_numb-1,browsA,bcolsA,Brows_index,Bcols_index);

	//printf ("\n\n\nC - browsA volta a apontar para o elemento 0 do vetor: \n\n");
	browsA-=1;

/*************** TRABALHANDO COM A **********************/
        //printf ("\n\n A - Nenhuma mudança é feita nos elementos 0 e 1: \n\n");

        bwdA = BRI_fwdA(_numb-1,browsA,bcolsA,Brows_index,Bcols_index) - bwdA;

        //printf ("\n\n A - Continua sem mudança... \n\n");

	/*************** RETORNANDO bwdA **********************/
	return bwdA;
    }
   else{         
	//printf ("\n\n\nChegou a Numb igual a 2 em A: \n\n");
        mat bwdA = gM(browsA[0],bcolsA[1],Brows_index,Bcols_index);
	bwdA = bwdA * pinv(gM(browsA[1],bcolsA[1],Brows_index,Bcols_index));
	bwdA = bwdA * gM(browsA[1],bcolsA[0],Brows_index,Bcols_index);
	bwdA = gM(browsA[0],bcolsA[0],Brows_index,Bcols_index) - bwdA;

    	return bwdA;
     }
}  

/*************************************** FUNÇÃO BRI_fwdB ********************************************/
mat BRI_fwdB(unsigned int _numb, unsigned int* browsB, unsigned int* bcolsB, unsigned int* Brows_index, unsigned int* Bcols_index){

    //printf ("\n\n\nEntrou na funcao:\n\n");
    if (_numb>2){
	/*************** pinvERTENDO bcolsB **********************/
	unsigned int temp = bcolsB[1];
	bcolsB[1] = bcolsB[0];
	bcolsB[0] = temp;

   //     EXPRESSÃO  bwdB = B - A*pinv(C)*D;

/*************** TRABALHANDO COM A **********************/
        //printf ("\n\n A - Nenhuma mudança é feita nos elementos 0 e 1: \n\n");

        mat bwdB = BRI_fwdA(_numb-1,browsB,bcolsB,Brows_index,Bcols_index);

        //printf ("\n\n A - Continua sem mudança... \n\n");

/*************** TRABALHANDO COM C **********************/
	//printf ("\n\n\nC - browsB passa a apontar para o elemento 1 do vetor: \n\n");
	browsB+=1;

        bwdB = bwdB * pinv(BRI_fwdC(_numb-1,browsB,bcolsB,Brows_index,Bcols_index));

	//printf ("\n\n\nC - browsB volta a apontar para o elemento 0 do vetor: \n\n");
	browsB-=1;

/*************** TRABALHANDO COM D **********************/
	//printf ("\n\n\nD - browsB e bcolsB passam a apontar para o elemento 1 do vetor: \n\n");
	browsB+=1;
	bcolsB+=1;

        bwdB = bwdB * BRI_fwdD(_numb-1,browsB,bcolsB,Brows_index,Bcols_index);

	//printf ("\n\n\nD - browsB e bcolsB voltam a apontar para o elemento 0 do vetor: \n\n");
	browsB-=1;
	bcolsB-=1;

/*************** TRABALHANDO COM B **********************/
        //printf ("\n\n\nB - bcolsB passa a apontar para o elemento 1 do vetor: \n\n");
	bcolsB+=1;

        bwdB = BRI_fwdB(_numb-1,browsB,bcolsB,Brows_index,Bcols_index) - bwdB;

        //printf ("\n\n\nB - bcolsB volta a apontar para o elemento 0 do vetor: \n\n");
	bcolsB-=1;


	/*************** pinvERTENDO NOVAMENTE bcolsB **********************/
	temp = bcolsB[1];
	bcolsB[1] = bcolsB[0];
	bcolsB[0] = temp;

	/*************** RETORNANDO bwdB **********************/
	return bwdB;
    }
   else{  
	//printf ("\n\n\nChegou a Numb igual a 2 em B: \n\n");
        mat bwdB = gM(browsB[0],bcolsB[0],Brows_index,Bcols_index);
	bwdB = bwdB * pinv(gM(browsB[1],bcolsB[0],Brows_index,Bcols_index));
	bwdB = bwdB * gM(browsB[1],bcolsB[1],Brows_index,Bcols_index);
	bwdB = gM(browsB[0],bcolsB[1],Brows_index,Bcols_index) - bwdB;

	return bwdB;
     }
} 

/*************************************** FUNÇÃO BRI_fwdC ********************************************/
mat BRI_fwdC(unsigned int _numb, unsigned int* browsC, unsigned int* bcolsC, unsigned int* Brows_index, unsigned int* Bcols_index){

    //printf ("\n\n\nEntrou na funcao:\n\n");
    if (_numb>2){
	/*************** pinvERTENDO browsC **********************/
	unsigned int temp = browsC[1];
	browsC[1] = browsC[0];
	browsC[0] = temp;

   //     EXPRESSÃO  bwdC = C - D*pinv(B)*A;

/*************** TRABALHANDO COM D **********************/
	//printf ("\n\n\nD - browsC e bcolsC passam a apontar para o elemento 1 do vetor: \n\n");
	browsC+=1;
	bcolsC+=1;

        mat bwdC = BRI_fwdD(_numb-1,browsC,bcolsC,Brows_index,Bcols_index);

	//printf ("\n\n\nD - browsC e bcolsC voltam a apontar para o elemento 0 do vetor: \n\n");
	browsC-=1;
	bcolsC-=1;

/*************** TRABALHANDO COM B **********************/
        //printf ("\n\n\nB - bcolsC passa a apontar para o elemento 1 do vetor: \n\n");
	bcolsC+=1;

        bwdC = bwdC * pinv(BRI_fwdB(_numb-1,browsC,bcolsC,Brows_index,Bcols_index));

        //printf ("\n\n\nB - bcolsC volta a apontar para o elemento 0 do vetor: \n\n");
	bcolsC-=1;

/*************** TRABALHANDO COM A **********************/
     	//printf ("\n\n A - Nenhuma mudança é feita nos elementos 0 e 1: \n\n");

        bwdC = bwdC * BRI_fwdA(_numb-1,browsC,bcolsC,Brows_index,Bcols_index);

        //printf ("\n\n A - Continua sem mudança... \n\n");

/*************** TRABALHANDO COM C **********************/
	//printf ("\n\n\nC - browsC passa a apontar para o elemento 1 do vetor: \n\n");
	browsC+=1;

        bwdC = BRI_fwdC(_numb-1,browsC,bcolsC,Brows_index,Bcols_index) - bwdC;

	//printf ("\n\n\nC - browsC volta a apontar para o elemento 0 do vetor: \n\n");
	browsC-=1;

	/*************** pinvERTENDO NOVAMENTE browsC **********************/
	temp = browsC[1];
	browsC[1] = browsC[0];
	browsC[0] = temp;

	/*************** RETORNANDO bwdC **********************/

	return bwdC;
    }
   else{         
	//printf ("\n\n\nChegou a Numb igual a 2 em C: \n\n");
        mat bwdC = gM(browsC[1],bcolsC[1],Brows_index,Bcols_index);
	bwdC = bwdC * pinv(gM(browsC[0],bcolsC[1],Brows_index,Bcols_index));
	bwdC = bwdC * gM(browsC[0],bcolsC[0],Brows_index,Bcols_index);
	bwdC = gM(browsC[1],bcolsC[0],Brows_index,Bcols_index) - bwdC;

	return bwdC;
     }
} 

/*************************************** FUNÇÃO BRI_fwdD ********************************************/
mat BRI_fwdD(unsigned int _numb, unsigned int* browsD, unsigned int* bcolsD, unsigned int* Brows_index, unsigned int* Bcols_index){

    //printf ("\n\n\nEntrou na funcao:\n\n");
    if (_numb>2){
	/*************** pinvERTENDO bcolsD e browsD **********************/
	//printf("ANTES: Imprimindo browsD[0]= %i e browsD[1]= %i \n",browsD[0], browsD[1] ); 	
	//printf("ANTES: Imprimindo bcolsD[0]= %i e bcolsD[1]= %i \n",bcolsD[0], bcolsD[1] ); 
	unsigned int temp = browsD[1];
	browsD[1] = browsD[0];
	browsD[0] = temp;
	

	//printf("DEPOIS: Imprimindo browsD[0]= %i e browsD[1]= %i \n",browsD[0], browsD[1] ); 	temp=0;
	temp = bcolsD[1];
	bcolsD[1] = bcolsD[0];
	bcolsD[0] = temp;
	//printf("DEPOIS: Imprimindo bcolsD[0]= %i e bcolsD[1]= %i \n",bcolsD[0], bcolsD[1] ); 

   //     EXPRESSÃO   bwdD = D - C*pinv(A)*B

/*************** TRABALHANDO COM C **********************/
	//printf ("\n\n\nC - browsD passa a apontar para o elemento 1 do vetor: \n\n");
	browsD+=1;

        mat bwdD = BRI_fwdC(_numb-1,browsD,bcolsD,Brows_index,Bcols_index);

	//printf ("\n\n\nC - browsD volta a apontar para o elemento 0 do vetor: \n\n");
	browsD-=1;


/*************** TRABALHANDO COM A **********************/
     	//printf ("\n\n A - Nenhuma mudança é feita nos elementos 0 e 1: \n\n");

        bwdD = bwdD * pinv(BRI_fwdA(_numb-1,browsD,bcolsD,Brows_index,Bcols_index));

        //printf ("\n\n A - Continua sem mudança... \n\n");

/*************** TRABALHANDO COM B **********************/
        //printf ("\n\n\nB - bcolsD passa a apontar para o elemento 1 do vetor: \n\n");
	bcolsD+=1;

        bwdD = bwdD * BRI_fwdB(_numb-1,browsD,bcolsD,Brows_index,Bcols_index);

        //printf ("\n\n\nB - bcolsD volta a apontar para o elemento 0 do vetor: \n\n");
	bcolsD-=1;


/*************** TRABALHANDO COM D **********************/
	//printf ("\n\n\nD - browsD e bcolsD passam a apontar para o elemento 1 do vetor: \n\n");
	browsD+=1;
	bcolsD+=1;

        bwdD = BRI_fwdD(_numb-1,browsD,bcolsD,Brows_index,Bcols_index) - bwdD;

	//printf ("\n\n\nD - browsD e bcolsD voltam a apontar para o elemento 0 do vetor: \n\n");
	browsD-=1;
	bcolsD-=1;

	/*************** pinvERTENDO NOVAMENTE bcolsD e browsD **********************/
	temp = browsD[1];
	browsD[1] = browsD[0];
	browsD[0] = temp;
	
	temp = bcolsD[1];
	bcolsD[1] = bcolsD[0];
	bcolsD[0] = temp;

	/*************** RETORNANDO bwdD **********************/
	return bwdD;
    }
   else{         
	//printf ("\n\n\nChegou a Numb igual a 2 em D: \n\n");
        mat bwdD = gM(browsD[1],bcolsD[0],Brows_index,Bcols_index);
	bwdD = bwdD * pinv(gM(browsD[0],bcolsD[0],Brows_index,Bcols_index));
	bwdD = bwdD * gM(browsD[0],bcolsD[1],Brows_index,Bcols_index);
	bwdD = gM(browsD[1],bcolsD[1],Brows_index,Bcols_index) - bwdD;

	return bwdD;
     }
}


void mexFunction(int nlhs, mxArray *plhs[],
         int nrhs, const mxArray *prhs[]) 
{   
    mat data = conv_to<mat>::from(armaGetPr(prhs[0], true));
    mat BRI_in = conv_to<mat>::from(armaGetPr(prhs[1], true));
    plhs[0] = armaCreateMxMatrix(BRI_in[1], BRI_in[1], mxDOUBLE_CLASS, mxREAL);
    armaSetPr(plhs[0], conv_to<mat>::from(init(data,BRI_in(0),BRI_in(1),BRI_in(2),BRI_in(3),BRI_in(4),BRI_in(5))));
}
