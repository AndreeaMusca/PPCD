#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#pragma comment(lib, "C:\\Program Files (x86)\\Microsoft SDKs\\MPI\\Lib\\x64\\msmpi.lib")

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

double sum_array(int n, double* v);
void sort_array(int n, double* v);
double* merge_array(int n, double* a, int m, double* b);
int MPI_sum_array(int n, double* v, double* sum, int root, MPI_Comm comm);
int MPI_sort_array(int n, double* v, int root, MPI_Comm comm);

int main(int argc, char** argv)
{
	//declarare variabile

	int rank, size, n = 100000;
	double* array, sum, t, allTime;

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
		srand(time(0));
		for (int i = 0; i < n; i++) {
			//array[i] = 1.;
			array[i] = rand() / 1000.;
		}
		/*for (int i = 0; i < n; i++) {
			printf("%.3lf ", array[i]);
		}
		printf("\n");*/
	}

	//apelare functie MPI_sum_array
	t = MPI_Wtime();
	//int error = MPI_sum_array(n, array, &sum, 0, MPI_COMM_WORLD);
	int error = MPI_sort_array(n, array, 0, MPI_COMM_WORLD);
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
		/*for (int i = 0; i < n; i++) {
			printf("%.3lf ", array[i]);
		}*/
		//puts("");
		//printf("Time = %lf\t sum = %.2lf\n", allTime, sum);
		printf("Time = %lf\n", allTime);
	}

	//eliberare memorie
	free(array);

	//finalizare context paralel
	MPI_Finalize();
	return 0;
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

void sort_array(int n, double* v)
{
	for (int i = 0; i < n-1; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			if (v[i] > v[j])
			{
				double aux = v[i];
				v[i] = v[j];
				v[j] = aux;
			}
		}
	}
}

double* merge_array(int n, double* a, int m, double* b)
{
	double* tmp = (double*)calloc(n + m, sizeof(double));
	int i=0, j=0, cnt=0;
	while (i < n && j < m)
	{
		if (a[i] < b[j])
			tmp[cnt++] = a[i++];
		else
			tmp[cnt++] = b[j++];
		while (i < n)
			tmp[cnt++] = a[i++];
		while (j < m)
			tmp[cnt++] = b[j++];
	}
	return tmp;
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

int MPI_sort_array(int n, double* v, int root, MPI_Comm comm)
{
	//declarare variabile
	int rank, size, local_n;
	double* local_array;

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

	//sort local_array
	sort_array(local_n, local_array);

	//gatter all local_array into array
	error = MPI_Gather(local_array, local_n, MPI_DOUBLE, v, local_n, MPI_DOUBLE, root, comm);

	if (error != MPI_SUCCESS)
	{
		free(local_array);
		return error;
	}

	if (rank == root)
	{
		for (int i = 1; i < size; i++)
		{
			double* tmp = merge_array(i * local_n, v, local_n, v + i * local_n);
			for (int j = 0; j < (i + 1) * local_n; j++)
			{
				v[j] = tmp[j];
			}
			free(tmp);
		}
	}

	free(local_array);
	return MPI_SUCCESS;
}


//La examen da o functie MPI_ceva la care
// sa descriem parametri si o bucatica de cod in care aratam cum se utilizeaza functia