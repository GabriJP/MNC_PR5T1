#include <cstdio>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "mkl_dss.h"
#include "mkl_types.h"

#define NNZ 9
#define N 5
#define NRHS 1

using namespace std;

int main(int argc, char *argv[]){
	int nnz = NNZ;
	int nrow = N;
	int ncol = N;
	double value[NNZ] = {1,2,3,4,5,6,7,8,9};
	int column[NNZ] = {1,1,2,1,3,2,4,3,5};
	int rowIndex[N + 1] = {1,2,4,6,8,10};

	double B[N] = { 1, 2, 3, 4, 5 };
	int nrhs = NRHS;

	double X[N*NRHS];

	_MKL_DSS_HANDLE_t handle;
	int opt = MKL_DSS_DEFAULTS;
	int sym = MKL_DSS_NON_SYMMETRIC;
	int type = MKL_DSS_INDEFINITE;
	int error;

	error = dss_create(handle, opt);
	if (error != MKL_DSS_SUCCESS){
		fprintf(stderr, "DSS Error. Crear el handle\n");
		exit(1);
	}

	error = dss_define_structure(handle, sym, rowIndex, nrow, ncol, column, nnz);
	if (error != MKL_DSS_SUCCESS){
		fprintf(stderr, "DSS Error. Definir la estructura\n");
		exit(1);
	}

	error = dss_reorder(handle, opt, 0);
	if (error != MKL_DSS_SUCCESS){
		fprintf(stderr, "DSS Error. Reordenar\n");
		exit(1);
	}


	error = dss_factor_real(handle, type, value);
	if (error != MKL_DSS_SUCCESS){
		fprintf(stderr, "DSS Error. Factorizacion\n");
		exit(1);
	}

	error = dss_solve_real(handle, opt, B, nrhs, X);
	if (error != MKL_DSS_SUCCESS){
		fprintf(stderr, "DSS Error. Resolucion\n");
		exit(1);
	}

	for (int i = 0; i < N; i++){
		printf("%.4f\n", X[i]);
	}

	_DOUBLE_PRECISION_t values[9];
	memset(values, 0, 9 * sizeof(_DOUBLE_PRECISION_t));
	error = dss_statistics(handle, opt, "ReorderTime,FactorTime,SolveTime,Determinant,Flops,Peakmem,Factormem", values);
	//ReorderTime
	printf("Tiempo utilizado en reordenar: %lf s\n", values[0]);
	//FactorTime
	printf("Tiempo utilizado en factorizar: %lf s\n", values[1]);
	//SolveTime
	printf("Tiempo utilizado en la resolucion: %lf s\n", values[2]);
	//Determinant
	printf("Determinante de A: %lfe%lf, es decir: %lf\n", values[4], values[3], values[4] * pow(10, values[3]));
	//Flops
	printf("Numero de FLOPS durante la factorizacion: %lf\n", values[6]);
	//Peakmem
	printf("Peakmem: %lf KiB\n", values[7]);
	//Factormem
	printf("Factormem: %lf KiB\n", values[8]);


	error = dss_delete(handle, opt);
	if (error != MKL_DSS_SUCCESS){
		fprintf(stderr, "DSS Error. Delete\n");
		exit(1);
	}


	char c = getc(stdin);
	exit(0);
}