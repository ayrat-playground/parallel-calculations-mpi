#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <io.h>
#include "mpi.h"

void matrixVectorMultiply(int argc, char* argv[]);
void hello_world(int argc, char* argv[]);
void vectorComponentsSum(int argc, char* argv[]);
void vectorMax(int argc, char* argv[]);
void maxMatrix(int argc, char* argv[]);
void dotProduct(int argc, char* argv[]);
void monteCarloPi(int argc, char* argv[]);
void averagePositiveCompon(int argc, char* argv[]);
void vectorInvert(int argc, char* argv[]);
void matrixVectorMultiplyColumn(int argc, char* argv[]);
int MatrixNormal(int argc, char **argv);
int MatrixInverse(int argc, char **argv);

int main(int argc, char* argv[]){

	//matrixVectorMultiply(argc, argv);
	//hello_world(argc, argv);
	//vectorComponentsSum(argc, argv);
	//vectorMax(argc, argv);
	//maxMatrix(argc, argv);
	//dotProduct(argc, argv);
	//monteCarloPi(argc, argv);
	//averagePositiveCompon(argc, argv);
	//matrixVectorMultiply(argc, argv);
	 MatrixNormal(argc, argv);

	
}

/*Hello World 2�����*/ 
void hello_world(int argc, char* argv[]){
	int ProcNum, ProcRank, RecvRank;
	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	if (ProcRank == 0){
		// ��������, ����������� ������ ��������� � ������ 0
		printf("\n Hello from process %3d", ProcRank);
		for (int i = 1; i < ProcNum; i++) {
			MPI_Recv(&RecvRank, 1, MPI_INT, MPI_ANY_SOURCE,
				MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
			printf("\n Hello from process %3d", RecvRank);
		}
	}
	else
		// ���������, ������������ ����� ����������,
		// ����� �������� � ������ 0
		MPI_Send(&ProcRank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	MPI_Finalize();
}

/*����� ������. �������*/


void vectorComponentsSum(int argc, char* argv[]){

	double x[100], TotalSum, ProcSum = 0.0;
	int ProcRank, ProcNum, N = 100;
	MPI_Status Status;
	// �������������
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	// ���������� ������
	if (ProcRank == 0) {
		for (int i = 0; i < N; i++){
			x[i] = 1;
		}
	}
	// �������� ������ �� ��� ��������
	MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// ���������� ��������� ����� �� ������ �� ���������
	// �� ������ �������� ����������� �������� ������� x �� i1 �� i2
	int k = N / ProcNum;
	int i1 = k * ProcRank;
	int i2 = k * (ProcRank + 1);
	if (ProcRank == ProcNum - 1) i2 = N;
	for (int i = i1; i < i2; i++)
		ProcSum = ProcSum + x[i];
	// ������ ��������� ���� �� �������� � ������ 0
	if (ProcRank == 0) {
		TotalSum = ProcSum;
		for (int i = 1; i < ProcNum; i++) {
			MPI_Recv(&ProcSum, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
				&Status);
			TotalSum = TotalSum + ProcSum;
		}
	}
	else // ��� �������� �������� ���� ��������� �����
		MPI_Send(&ProcSum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	// ����� ����������
	if (ProcRank == 0)
		printf("\nTotal Sum = %10.2f", TotalSum);
	MPI_Finalize();
	system("pause");

}

/*2.	Max �������(3 �����)*/
void vectorMax(int argc, char* argv[]){

	double x[100], TotalMax, ProcMax = 0.0;
	int ProcRank, ProcNum, N = 100;
	MPI_Status Status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0) {
		for (int i = 0; i < N; i++){
			int max = 1000;
			int min = 0;
			x[i] = min + (rand() % (int)(max - min + 1));
		}
	}

	MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	int k = N / ProcNum;
	int i1 = k * ProcRank;
	int i2 = k * (ProcRank + 1);
	if (ProcRank == ProcNum - 1) i2 = N;
	for (int i = i1; i < i2; i++)
		if (x[i]>ProcMax){
			ProcMax = x[i];
		}

	MPI_Reduce(&ProcMax, &TotalMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if (ProcRank == 0)
		printf("\nTotal Max = %10.2f", TotalMax);
	MPI_Finalize();
	system("pause");

}

/*2.	Max �������(5 �����)*/
void maxMatrix(int argc, char* argv[]){

	double x[100][100], TotalMax, ProcMax = 0.0;
	int ProcRank, ProcNum, N = 100;
	MPI_Status Status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0) {
		int max = 1000;
		int min = 0;
		for (int i = 0; i < N; i++){
			for (int j = 0; j < N; j++){

				x[i][j] = min + (rand() % (int)(max - min + 1));
			}
		}
	}
	MPI_Bcast(x, N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	int k = N / ProcNum;
	int i1 = k * ProcRank;
	int i2 = k * (ProcRank + 1);
	if (ProcRank == ProcNum - 1) i2 = N;
	for (int i = i1; i < i2; i++)
		for (int j = 0; j<N; j++){
			if (x[i][j]>ProcMax){
				ProcMax = x[i][j];
			}
		}

	MPI_Reduce(&ProcMax, &TotalMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if (ProcRank == 0)
		printf("\nTotal Max = %10.2f", TotalMax);
	MPI_Finalize();
	system("pause");

}


/*4.	��������� ������������(3 �����)*/
void dotProduct(int argc, char* argv[]){

	double x[100],y[100], TotalSum, ProcSum = 0.0;
	int ProcRank, ProcNum, N = 100;
	MPI_Status Status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0) {
		int max = 1000;
		int min = 0;
		for (int i = 0; i < N; i++){
			

				x[i] = 1;
				y[i] = 1;
			
		}
	}
	MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(y, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	int k = N / ProcNum;
	int i1 = k * ProcRank;
	int i2 = k * (ProcRank + 1);
	if (ProcRank == ProcNum - 1) i2 = N;
	for (int i = i1; i < i2; i++){
		ProcSum = ProcSum + x[i] * y[i];
		}

	MPI_Reduce(&ProcSum, &TotalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (ProcRank == 0)
		printf("\nTotal Sum = %10.2f", TotalSum);
	MPI_Finalize();
	system("pause");

}

/*16.	���������� ����� �� ������� �����-����� ( 5 ������)     */
void monteCarloPi(int argc, char* argv[]){

	double x, y, pi;
	int ProcRank, ProcNum, N = 1000000, in = 0,allIn;
	MPI_Status Status;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);


	

	int k = N / ProcNum;
	srand(time(NULL));   // change seed value for random number generation
	for (int i = 1; i <= k; i++){
		x = (double) rand() / RAND_MAX;;
		y = (double)rand() / RAND_MAX;;
		if (x*x + y*y <= 1)
			in++;
	}

	MPI_Reduce(&in, &allIn, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (ProcRank == 0){
		pi =(double) 4*allIn/N;
		printf("\nPI = %10.2f", pi);
	}
	MPI_Finalize();

	system("pause");
	

}

/*14.	������� �������������� ����� ������������� ����� ������� */
void averagePositiveCompon(int argc, char* argv[]){

	double x[100], TotalAver, ProcAver = 0.0;
	int ProcRank, ProcNum, N = 100, TotalCount, ProcCount = 0;
	MPI_Status Status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0) {
		for (int i = 0; i < N; i++){
			int max = 1000;
			int min = -1000;
			x[i] = min + (rand() % (int)(max - min + 1));
		}
	}

	MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	int k = N / ProcNum;
	int i1 = k * ProcRank;
	int i2 = k * (ProcRank + 1);
	if (ProcRank == ProcNum - 1) i2 = N;

	for (int i = i1; i < i2; i++)
		if (x[i]>0){
			ProcAver=+x[i];
			ProcCount++;
			
		}


	MPI_Reduce(&ProcAver, &TotalAver, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&ProcCount, &TotalCount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Finalize();
	double average = (double)TotalAver / TotalCount;
	printf("\nTotal Average = %10.2f", average);
	system("pause");

}

/* ��������� ������� �� ������ �� �������  */
void matrixVectorMultiply(int argc, char* argv[]){

	double x[100][100],y[100],z[100], TotalProd[100];
	int ProcRank, ProcNum, N = 100;
	MPI_Status Status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0) {
		int max = 1000;
		int min = 0;
		for (int i = 0; i < N; i++){
			z[i] = 0;
			y[i] = 1;
			for (int j = 0; j < N; j++){

				x[i][j] = 1;
			}
		}
	}

	MPI_Bcast(x, N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(y, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	int k = N / ProcNum;
	int i1 = k * ProcRank;
	int i2 = k * (ProcRank + 1);
	if (ProcRank == ProcNum - 1) i2 = N;
	for (int i = i1; i < i2; i++)
		for (int j = 0; j<N;j++){
			z[i] =+ y[j] * x[i][j];
		}

	MPI_Reduce(&z, &TotalProd, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	

	MPI_Finalize();
	for (int i = 0; i < N; i++){
		printf("\nProds[%i] = %10.2f", i, TotalProd[i]);
	}
	system("pause");


}

/*15.	������������� ������(7 ������) */
void vectorInvert(int argc, char* argv[]){

	double x[100], z[100], y[100];
	int ProcRank, ProcNum, N = 100;
	MPI_Status Status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0) {
		for (int i = 0; i < N; i++){
			int max = 1000;
			int min = 0;
			y[i] = 0;
			x[i] = i;
		}
	}

	MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(y, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	int k = N / ProcNum;
	int i1 = k * ProcRank;
	int i2 = k * (ProcRank + 1);
	if (ProcRank == ProcNum - 1) i2 = N;
	for (int i = i1; i < i2; i++)
		y[N - 1 - i] = x[i];

	MPI_Reduce(&y, &z, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Finalize();

	for (int i = 0; i < N; i++){
		printf("\nProds[%i] = %10.2f", i, z[i]);
	}
	system("pause");

}
/* ��������� ������� �� ������ �� ��������*/
void matrixVectorMultiplyColumn(int argc, char* argv[]){

	double x[100][100], y[100], z[100],d[100], TotalRes[100], res[100];
	int ProcRank, ProcNum, N = 100;
	MPI_Status Status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0) {
		int max = 1000;
		int min = 0;
		for (int i = 0; i < N; i++){
			z[i] = 0;
			y[i] = 1;
			for (int j = 0; j < N; j++){

				x[i][j]=1;
			}
		}
	}

	MPI_Bcast(x, N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(y, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(z, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	int k = N / ProcNum;
	int i1 = k * ProcRank;
	int i2 = k * (ProcRank + 1);
	if (ProcRank == ProcNum - 1) i2 = N;
	for (int i = i1; i < i2; i++)
		for (int j = 0; j<N; j++){
			z[j] = +y[j] * x[j][i];

		}
	MPI_Reduce(&z, &TotalRes, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


	MPI_Finalize();
	for (int i = 0; i < N; i++){
		printf("\nProds[%i] = %10.2f", i, TotalRes[i]);
	}
	system("pause");


}

int MatrixNormal(int argc, char **argv) {
	int myrank, total, N_PER_PROC=2;

	double *A, *B, *C;	// ������������ ������ � root

	double *a, *b, *c;	// ����� ������� [mxn], ������ [n], ���-� [m]
	int n, m;
	int i, j;
	int intBuf[2];

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &total);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	printf("Total=%d, rank=%d\n", total, myrank);

	if (!myrank) {	// ���������� ���. ������ (������ root)
		n = N_PER_PROC * total;
		A = (double *)malloc(sizeof(double)*n*n);
		B = (double *)malloc(sizeof(double)*n);
		C = (double *)malloc(sizeof(double)*n);
		for (i = 0; i<n; i++) {
			B[i] = (double)i;
			for (j = 0; j<n; j++)
				A[i*n + j] = (double)(i + j);
		};
	};

	if (!myrank) {
		intBuf[0] = n;
		intBuf[1] = N_PER_PROC;
	};
	MPI_Bcast((void *)intBuf, 2, MPI_INT, 0, MPI_COMM_WORLD);
	n = intBuf[0];
	m = intBuf[1];

	a = (double *)malloc(sizeof(double)*n*m);
	b = (double *)malloc(sizeof(double)*n);
	c = (double *)malloc(sizeof(double)*m);

	if (!myrank) {	// ������ ��������, B �� �����
		for (int i = 0; i < n; i++){
			b[i] = B[i];
		}
	};
	MPI_Bcast((void *)b, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Scatter((void *)A, n*m, MPI_DOUBLE,
		(void *)a, n*m, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (i = 0; i<m; i++) {
		c[i] = 0;
		for (j = 0; j<n; j++)
			c[i] += a[n*i + j] * b[j];
	};

	MPI_Gather((void *)c, m, MPI_DOUBLE,
		(void *)C, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (!myrank)
		for (i = 0; i<n; i++)
			printf("%g\n", C[i]);

	MPI_Finalize();
	system("pause");
	exit(0);
}

int MatrixInverse(int argc, char **argv) {
	int myrank, total, N_PER_PROC = 2;

	double *A, *B, *C;	// ������������ ������ � root

	double *a, *b, *c;	// ����� ������� [mxn], ������ [n], ���-� [m]
	int n, m;
	int i, j;
	int intBuf[2];

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &total);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	printf("Total=%d, rank=%d\n", total, myrank);

	if (!myrank) {	// ���������� ���. ������ (������ root)
		n = N_PER_PROC * total;
		A = (double *)malloc(sizeof(double)*n*n);
		B = (double *)malloc(sizeof(double)*n);
		C = (double *)malloc(sizeof(double)*n);
		for (i = 0; i<n; i++) {
			B[i] = (double)i;
			for (j = 0; j<n; j++)
				A[i*n + j] = (double)(i + j);
		};
		for (i = 0; i<n; i++) {
			B[i] = (double)i;
			for (j = 0; j<n; j++)
				A[i*n + j] = (double)(i + j);
		};
		for ( i = 0; i < n; i++)
		{
			for (j = 0; j < i; j++)
			{
				double tmp = A[i*n + j];
				A[i*n + j] = A[j*n + i];
				A[j*n + i] = tmp;
			}
		}
	};

	if (!myrank) {
		intBuf[0] = n;
		intBuf[1] = N_PER_PROC;
	};
	MPI_Bcast((void *)intBuf, 2, MPI_INT, 0, MPI_COMM_WORLD);
	n = intBuf[0];
	m = intBuf[1];

	a = (double *)malloc(sizeof(double)*n*m);
	b = (double *)malloc(sizeof(double)*n);
	c = (double *)malloc(sizeof(double)*m);

	if (!myrank) {	// ������ ��������, B �� �����
		for (int i = 0; i < n; i++){
			b[i] = B[i];
		}
	};
	MPI_Bcast((void *)b, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Scatter((void *)A, n*m, MPI_DOUBLE,
		(void *)a, n*m, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (i = 0; i<m; i++) {
		for (j = 0; j<n; j++)
			c[j] = a[n*i + j] * b[j];
	};

	MPI_Gather((void *)c, m, MPI_DOUBLE,
		(void *)C, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (!myrank)
		for (i = 0; i<n; i++)
			printf("%g\n", C[i]);

	MPI_Finalize();
	system("pause");
	exit(0);
}