#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#pragma comment(lib, "C:\\Program Files (x86)\\Microsoft SDKs\\MPI\\Lib\\x64\\msmpi.lib")

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>


int** alocareMatrice(int n, int m);
void free_mat(int*** mat);
void initMatrix(int n, int m, int** mat);
void initMatrix_eye(int n, int m, int** mat);
void printf_matrix(int n, int m, int** mat);
int** transp(int n, int m, int** mat);
int** prod_matrix(int n, int m, int p, int** a, int** b);
int** prod_matrix_pseudo(int n, int m, int p, int** a, int** b);

int MPI_Prod_matrix(int n, int m, int p, int** a, int** b, int** c, int root, MPI_Comm comm);
int MPI_Prod_matrix_pseudo(int n, int m, int p, int** a, int** b, int** c, int root, MPI_Comm comm);
int MPI_Prod_matrix_pseudo_row(int n, int m, int p, int** a, int** b, int** c, int root, MPI_Comm comm);
int main(int argc, char** argv)
{
	//declarare variabile
	int rank, size, m = 1500, n = 1500, p = 1500;
	int** a, ** b, ** c;
	srand(time(0));
	//get rank si size
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//alocare memorie matrici
	a = alocareMatrice(n, m);
	b = alocareMatrice(m, p);
	c = alocareMatrice(n, p);
	//initializare matrici
	if (rank == 0)
	{
		initMatrix_eye(n, m, a);
		//printf("Mat a:\n");
		//printf_matrix(n, m, a);
		initMatrix(m, p, b);
		//printf("Mat b:\n");
		//printf_matrix(m, p, b);
	}

	//apelare functie MPI_Prod_matrix
	double t = MPI_Wtime(), alltime;
	//MPI_Prod_matrix(n, m, p, a, b, c, 0, MPI_COMM_WORLD);
	//MPI_Prod_matrix_pseudo(n, m, p, a, b, c, 0, MPI_COMM_WORLD);
	MPI_Prod_matrix_pseudo_row(n, m, p, a, b, c, 0, MPI_COMM_WORLD);
	t = MPI_Wtime() - t;
	MPI_Reduce(&t, &alltime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank == 0)
	{
		printf("Time = %lf\n", alltime);
		//printf("Mat c:\n");
		//printf_matrix(n, p, c);
	}
	//eliberare memorie
	free_mat(&a);
	free_mat(&b);
	free_mat(&c);
	//finalizare context paralel
	MPI_Finalize();
	return 0;
}

int** alocareMatrice(int n, int m)
{
	int** mat = (int**)calloc(n, sizeof(int*));
	mat[0] = (int*)calloc(n * m, sizeof(int));
	for (int i = 1; i < n; i++)
		mat[i] = mat[i - 1] + m;
	return mat;
}

void free_mat(int*** mat)
{
	free((*mat)[0]);
	free(*mat);
}
void initMatrix(int n, int m, int** mat)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			mat[i][j] = rand() % 10;
}
void initMatrix_eye(int n, int m, int** mat)
{
	n = n < m ? n : m;
	for (int i = 0; i < n; i++)
		mat[i][i] = 1;
}

void printf_matrix(int n, int m, int** mat)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
			printf("%d ", mat[i][j]);
		puts("");
	}
}

int** transp(int n, int m, int** mat)
{
	int** mat_t = alocareMatrice(m, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			mat_t[j][i] = mat[i][j];
	return mat_t;
}

int** prod_matrix(int n, int p, int m, int** a, int** b)
{
	int** mat = alocareMatrice(n, m);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
		{
			for (int k = 0; k < p; k++)
				mat[i][j] += a[i][k] * b[k][j];
		}
	return mat;
}

int** prod_matrix_pseudo(int n, int m, int p, int** a, int** b)
{
	int** mat = alocareMatrice(n, m);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < p; j++)
		{
			for (int k = 0; k < m; k++)
				mat[i][j] += a[i][k] * b[j][k];
		}
	return mat;
}

int MPI_Prod_matrix(int n, int m, int p, int** a, int** b, int** c, int root, MPI_Comm comm)
{
	//declarare variabile
	int rank, size, ** local_a, ** local_c, local_rows;
	//get rank si size

	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	local_rows = n / size;
	//alocare memorie pentru matrici locale
	local_a = alocareMatrice(local_rows, m);

	//scatter pentru matricea a
	int error = MPI_Scatter(a[0], local_rows * m, MPI_INT, local_a[0], local_rows * m, MPI_INT, root, comm);
	if (error != MPI_SUCCESS)
	{
		free_mat(&local_a);
		return error;
	}
	//bcast pentru matricea b
	error = MPI_Bcast(b[0], m * p, MPI_INT, root, comm);
	if (error != MPI_SUCCESS)
	{
		free_mat(&local_a);
		return error;
	}
	//inmultire local_a cu b
	local_c = prod_matrix(local_rows, m, p, local_a, b);

	//gather local_res in c
	error = MPI_Gather(local_c[0], local_rows * p, MPI_INT, c[0], local_rows * p, MPI_INT, root, comm);
	//eliberare memorie
	free_mat(&local_a);
	free_mat(&local_c);
	return error;

}

int MPI_Prod_matrix_pseudo(int n, int m, int p, int** a, int** b, int** c, int root, MPI_Comm comm)
{
	//declarare variabile
	int rank, size, ** local_a, ** local_c, local_rows;
	//get rank si size

	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	local_rows = n / size;
	//alocare memorie pentru matrici locale
	local_a = alocareMatrice(local_rows, m);

	//scatter pentru matricea a
	int error = MPI_Scatter(a[0], local_rows * m, MPI_INT, local_a[0], local_rows * m, MPI_INT, root, comm);
	if (error != MPI_SUCCESS)
	{
		free_mat(&local_a);
		return error;
	}
	//bcast pentru matricea b
	error = MPI_Bcast(b[0], m * p, MPI_INT, root, comm);
	if (error != MPI_SUCCESS)
	{
		free_mat(&local_a);
		return error;
	}
	//inmultire local_a cu b
	int** bT = transp(m, p, b);
	local_c = prod_matrix_pseudo(local_rows, m, p, local_a, bT);
	free(bT);
	//gather local_res in c
	error = MPI_Gather(local_c[0], local_rows * p, MPI_INT, c[0], local_rows * p, MPI_INT, root, comm);
	//eliberare memorie
	free_mat(&local_a);
	free_mat(&local_c);
	return error;
}


int MPI_Prod_matrix_pseudo_row(int n, int m, int p, int** a, int** b, int** c, int root, MPI_Comm comm)
{
	//declarare variabile
	int rank, size, ** local_a, ** local_c, local_rows;
	//get rank si size

	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	local_rows = n / size;

	MPI_Datatype row;
	MPI_Type_contiguous(m, MPI_INT, &row);
	MPI_Type_commit(&row);

	//alocare memorie pentru matrici locale
	local_a = alocareMatrice(local_rows, m);

	//scatter pentru matricea a
	int error = MPI_Scatter(a[0], local_rows, row, local_a[0], local_rows, row, root, comm);
	if (error != MPI_SUCCESS)
	{
		free_mat(&local_a);
		return error;
	}
	//bcast pentru matricea b
	error = MPI_Bcast(b[0], p, row, root, comm);
	if (error != MPI_SUCCESS)
	{
		free_mat(&local_a);
		return error;
	}
	//inmultire local_a cu b
	int** bT = transp(m, p, b);
	local_c = prod_matrix_pseudo(local_rows, m, p, local_a, bT);
	free(bT);
	//gather local_res in c
	error = MPI_Gather(local_c[0], local_rows * p, MPI_INT, c[0], local_rows * p, MPI_INT, root, comm);
	//eliberare memorie
	MPI_Type_free(&row);
	free_mat(&local_a);
	free_mat(&local_c);
	return error;
}