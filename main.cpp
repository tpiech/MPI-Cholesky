// MPITest.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h" - for Visual Studio 
#include <iostream>
#include <cmath>
#include "mpi.h"
using namespace std;

double** create_d2_matrix_nxn(int n) {
	double **new_matrix = new double*[n];
	for (int i = 0; i < n; i++) {
		new_matrix[i] = new double[n];
	}
	return new_matrix;
}
double* create_matrix_nxn(int n) {
	double *new_matrix = new double[n*n];
	return new_matrix;
}
void print_matrix(double* matrix, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << matrix[i*n + j] << " ";
		}
		cout << endl;
	}
}


double determinant(double **matrixb, int n) {
	double det = 0;
	double **matrix = create_d2_matrix_nxn(n);
	double **submatrix = create_d2_matrix_nxn(n);
	for (int i = 0; i<n; i++) {
		for (int j = 0; j<n; j++) {
			matrix[i][j] = matrixb[i][j];
		}

	}

	if (n == 2)
		return ((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]));
	else {
		for (int x = 0; x < n; x++) {
			int subi = 0;
			for (int i = 1; i < n; i++) {
				int subj = 0;
				for (int j = 0; j < n; j++) {
					if (j == x)
						continue;
					submatrix[subi][subj] = matrix[i][j];
					subj++;
				}
				subi++;
			}
			det = det + (pow(-1, x) * matrix[0][x] * determinant(submatrix, n - 1));
		}
	}
	return det;
}




double determinant_flat(double *matrix_a, int n) {

	double det = 0;
	double **matrix = create_d2_matrix_nxn(n);
	for (int i = 0; i<n; i++) {
		for (int j = 0; j<n; j++) {
			matrix[i][j] = matrix_a[i*n + j];
		}

	}

	double **submatrix = create_d2_matrix_nxn(n);
	if (n == 2)
		return ((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]));
	else {
		for (int x = 0; x < n; x++) {
			int subi = 0;
			for (int i = 1; i < n; i++) {
				int subj = 0;
				for (int j = 0; j < n; j++) {
					if (j == x)
						continue;
					submatrix[subi][subj] = matrix[i][j];
					subj++;
				}
				subi++;
			}
			det = det + (pow(-1, x) * matrix[0][x] * determinant(submatrix, n - 1));
		}
	}
	return det;
}



bool check_if_matrix_is_symetric(double* matrix, int n) {
	for (int i = 0; i<n; i++) {
		for (int j = i; j<n; j++) {
			if (i != j) {
				if (matrix[i*n + j] != matrix[j * n + i]) {
					return false;
				}
			}
		}
	}
	return true;
}
void remove_matrix(double* matrix, int n) {
	delete[] matrix;
	cout << "matrix removed" << endl;
}

void Choleski_up(double *a, double *l, int n)
{
	double sum = 0.0;

	for (int s = 0; s<n; s++)
	{
		for (int i = s; i<n; i++)
		{
			sum = 0.0;
			for (int j = 0; j<s; j++)
				sum += l[j*n + i] * l[j*n + s];
			l[s *n + i] = a[s*n + i] - sum;
			if (s == i) l[s*n + i] = sqrt((l[i*n + s]));
			else l[s*n + i] = l[s*n + i] / l[s*n + s];
		}
	}
}


void Choleski_down(double *a, double *l, int n)
{
	double sum = 0.0;

	for (int s = 0; s<n; s++)
	{
		for (int i = s; i<n; i++)
		{
			sum = 0.0;
			for (int j = 0; j<s; j++)
				sum += l[i* n + j] * l[s* n + j];
			l[i*n + s] = a[i*n + s] - sum;
			if (s == i) l[i*n + s] = sqrt((l[i*n + s]));
			else l[i*n + s] = l[i*n + s] / l[s*n + s];
		}
	}
}


int main(int argc, char *argv[])
{

	MPI_Init(&argc, &argv);
	//MPI_Init(NULL, NULL);
	int rank;
	int size;
	double matrix_det;
	int is_matrix_symetric;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int n;
	if (rank == 0) {
		cout << "Please input size of matrix, you want to decompose (n x n)" << endl;
		cin >> n;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	double *a = new double[n*n];
	double *l = new double[n*n];
	double *l_2 = new double[n*n];
	if (rank == 0)
	{

		cout << "Please input matrix, you want to decompose" << endl;
		for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
		cout << "a[" << i << "][" << j << "] = ";
		//cin >> a[i][j];
		cin >> a[i*n + j];
		cout << endl;
		}
		} 

		print_matrix(a, n);
		cout << endl;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				l[i*n + j] = 0;
			}
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				l_2[i* n + j] = 0;
			}
		}
		MPI_Send(l, n*n, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD);
		MPI_Send(l_2, n*n, MPI_DOUBLE, 4, 0, MPI_COMM_WORLD);
	}

	
	MPI_Bcast(a, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	if (rank == 1)
	{ // proces 1 - liczy wyznacznik macierzy - otrzymuje macierz A i wymiar n z 0, wysyla wynik wyznacznika do procesu 3
	  //recive n
	  //recive a

		matrix_det = determinant_flat(a, n);
		cout << "Matrix determinant = " << matrix_det << endl;
		//send det to 3 i 4
		MPI_Send(&matrix_det, 1, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD);
		MPI_Send(&matrix_det, 1, MPI_DOUBLE, 4, 0, MPI_COMM_WORLD);
	}

	if (rank == 2) // proces 2 - sprawdza czy macierz jest symetryczna - otrzymuje macierz A i wymiar n z 0, wysyla wynik funkcji do procesu 3
	{
		is_matrix_symetric = check_if_matrix_is_symetric(a, n);
		if (is_matrix_symetric == 1) {
			cout << "Matrix is symetric";
		}
		else {
			cout << "Matrix is not symetric";
		}
		MPI_Send(&is_matrix_symetric, 1, MPI_INT, 3, 0, MPI_COMM_WORLD);
		MPI_Send(&is_matrix_symetric, 1, MPI_INT, 4, 0, MPI_COMM_WORLD);
	}
	if (rank == 3) // dokonuje rozkladu do m_dolnej, otrzymuje A, n, is_m_sym i det, wypisuje wynik
	{	// recive l
		// recive is_sym
		// recive det
		MPI_Recv(&is_matrix_symetric, 1, MPI_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&matrix_det, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(l, n*n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (is_matrix_symetric && matrix_det > 0) {
			cout << "Matrix is ok, let's create Cholesky lower matrix" << endl;
			Choleski_down(a, l, n);
			print_matrix(l, n);
			cout << endl;
		}
		else {
			cout << "Matrix is not symetric or not positive oriented :(" << endl;
		}


	}

	if (rank == 4) {
		// dokonuje rozkladu do m_gornej, otrzymuje A, n, is_m_sym i det, 
		// recive l_2
		//recive A
		// recive n
		// recive is_sym
		// recive det
		MPI_Recv(&is_matrix_symetric, 1, MPI_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&matrix_det, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(l_2, n*n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (is_matrix_symetric && matrix_det > 0) {
			cout << "Matrix is ok, let's create Cholesky upper matrix!" << endl;
			Choleski_up(a, l_2, n);
			cout << endl;
			print_matrix(l_2, n);
			cout << endl;

		}
		else {
			cout << "Matrix is not symetric or not positive oriented :(" << endl;
		}

	}

	MPI_Finalize();
	return 0;
}
