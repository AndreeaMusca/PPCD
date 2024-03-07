#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#pragma comment(lib, "C:\\Program Files (x86)\\Microsoft SDKs\\MPI\\Lib\\x64\\msmpi.lib")

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double sum_array(int n, double* v);
int MPI_sum_array(int n, double* v, double* sum, int root, MPI_Comm comm);

int main(int argc, char** argv)
{
	//declarare variabile

	int rank, size, n = 10000000;
	double* array;
	double sum, t, allTime;

	//initializare context paralel
	MPI_Init(&argc, &argv);

	//get rank si size
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//alocare vector
	array = (double*)calloc(n, sizeof(double));

	//initializare vector pe procesorul radacina
	if (rank == 0)
	{
		for (int i = 0; i < n; i++) {
			array[i] = 1.;
		}
	}

	//apelare functie MPI_sum_array
	t = MPI_Wtime();
	int error = MPI_sum_array(n, array, &sum, 0, MPI_COMM_WORLD);
	if (error != MPI_SUCCESS)
	{
		free(array);
		MPI_Finalize();
		return 1;
	}

	t = MPI_Wtime() - t;

	//afisare timp de executite
	error = MPI_Reduce(&t, &allTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank == 0)
	{
		printf("Time = %lf\t sum = %.2lf\n", allTime, sum);
	}

	//eliberare memorie
	free(array);

	//finalizare context paralel
	MPI_Finalize();
}

double sum_array(int n, double* v)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
	{
		sum += v[i];
	}
	return sum;
}

int MPI_sum_array(int n, double* v, double* sum, int root, MPI_Comm comm)
{
	//declarare variabile
	int rank, size, local_n;
	double* local_array, local_sum;

	//get rank, size
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	//alocare local_array
	local_n = n / size;
	local_array = (double*)calloc(local_n, sizeof(double));

	//scatter array in  local_array
	int error = MPI_Scatter(v, local_n, MPI_DOUBLE, local_array, local_n, MPI_DOUBLE, root, comm);
	if (error != MPI_SUCCESS)
	{
		free(local_array);
		return error;
	}

	//calculare local_sum
	local_sum = sum_array(local_n, local_array);

	//reduce local_sum in sum
	error = MPI_Reduce(&local_sum, sum, 1, MPI_DOUBLE, MPI_SUM, root, comm);

	//eliberare memorie
	free(local_array);

	//return
	return error;
}


//La examen da o functie MPI_ceva la care
// sa descriem parametri si o bucatica de cod in care aratam cum se utilizeaza functia