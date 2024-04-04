#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#pragma comment(lib, "C:\\Program Files (x86)\\Microsoft SDKs\\MPI\\Lib\\x64\\msmpi.lib")

#include <stdio.h>
#include<malloc.h>
#include <stdlib.h>
#include <time.h>

//MPI methods
int MPI_Exchange(int n, double* array, int rank1, int rank2, MPI_Comm comm);
int MPI_Sort(int n, double* array, int root, MPI_Comm comm);
int MPI_Is_Sorted(int n, double* array, int* isSorted, int root, MPI_Comm comm);

//all in previous labs
double* merge(int n, double* array, int m, double* b);
void merge_sort(int n, double* array);
void swap(double* array, double* b);

// function definitions
int main(int argc, char** argv) {
    //setup
    int size, rank, result, i, * answer;
    int n = 5000000;
    double m = 10.0;
    double* array, processorTime;

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //allocate space for an array of doubles, size n
    array = (double*)calloc(n, sizeof(double));

    //fills array with random values on root proc
    if (rank == 0)
    {
        //get random values for array & output for testing
        srand(((unsigned)time(NULL) + rank));
        for (i = 0; i < n; i++)
        {
            array[i] = ((double)rand() / RAND_MAX) * m;
            //printf("Initial: %f\n", array[i]);
        }
    }

    //get start time for each processor
    processorTime = MPI_Wtime();

    //MPI_Sort does all the heavy work
    result = MPI_Sort(n, array, 0, MPI_COMM_WORLD);
    if (result != MPI_SUCCESS)
    {
        return result;
    }

    //get end time for each processor
    processorTime = MPI_Wtime() - processorTime;
    //TODO COMMENT OUT
    //printf("Processor %d takes %lf sec\n", rank, processorTime);

    double allTime;
    MPI_Reduce(&processorTime, &allTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    //output ordered list for testing
    if (rank == 0)
    {
        printf("Time: %lf\n", allTime);
        //for (i = 0; i < n; i++)
        //{
        //    //TODO COMMENT OUT
        //    printf("Output : %f\n", array[i]);
        //}
    }
    free(array);
    MPI_Finalize();
}

int MPI_Sort(int n, double* array, int root, MPI_Comm comm) {
    int rank, size;
    double* local_array;
    // get rank and size of comm
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    //allocate space for numElements/numProcessors amount of doubles
    local_array = (double*)calloc(n / size, sizeof(double));

    //scatter a to local_a
    MPI_Scatter(array, n / size, MPI_DOUBLE, local_array, n / size, MPI_DOUBLE, root, comm);

    //sort local_a using mergeSort
    merge_sort(n / size, local_array);

    //odd-even iterations

    int isSorted;

    for (int i = 0; i < size; i++)
    {
        if ((i + rank) % 2 == 0)
        {
            if (rank + 1 < size)
                MPI_Exchange(n / size, local_array, rank, rank + 1, comm);
        }
        else
        {
            if (rank > 0)
                MPI_Exchange(n / size, local_array, rank - 1, rank, comm);
        }
        MPI_Barrier(comm);

        MPI_Is_Sorted(n / size, local_array, &isSorted, root, comm);
        if (isSorted)
        {
            if (rank == root)
            {
                printf("Iteratii: %d\n", i);
            }
            break;
        }
    }

    //gather local_a
    MPI_Gather(local_array, n / size, MPI_DOUBLE, array, n / size, MPI_DOUBLE, root, comm);
    free(local_array);
    return MPI_SUCCESS;

}

int MPI_Exchange(int n, double* array, int rank1, int rank2, MPI_Comm comm) {
    int rank, size, result, i, tag1 = 0, tag2 = 1;
    double* b = (double*)calloc(n, sizeof(double));
    double* c = NULL;

    MPI_Status status;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    //L8.6
    if (rank == rank1)
    {
        result = MPI_Send(&array[0], n, MPI_DOUBLE, rank2, tag1, comm);
        result = MPI_Recv(&b[0], n, MPI_DOUBLE, rank2, tag2, comm, &status);
        c = merge(n, array, n, b);
        for (i = 0; i < n; i++)
        {
            array[i] = c[i];
        }
    }
    else if (rank == rank2)
    {
        result = MPI_Recv(&b[0], n, MPI_DOUBLE, rank1, tag1, comm, &status);
        result = MPI_Send(&array[0], n, MPI_DOUBLE, rank1, tag2, comm);
        c = merge(n, array, n, b);
        for (i = 0; i < n; i++)
        {
            array[i] = c[i + n];
        }
    }
    free(b);
    if (c)
        free(c);
    return MPI_SUCCESS;
}




int MPI_Is_Sorted(int n, double* array, int* isSorted, int root, MPI_Comm comm)
{
    int rank, size;
    double* firsts, *lasts;

    // get rank and size of comm
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    //alocate memory
    firsts = (double*)calloc(size, sizeof(double));
    lasts = (double*)calloc(size, sizeof(double));

    MPI_Gather(&array[0], 1, MPI_DOUBLE, firsts, 1, MPI_DOUBLE, root, comm);
    MPI_Gather(&array[n-1], 1, MPI_DOUBLE, lasts, 1, MPI_DOUBLE, root, comm);

    *isSorted = 1;

    for(int i = 0; i < size - 1; i++)
    {
        if (lasts[i] > firsts[i + 1])
            *isSorted = 0;
    }

    MPI_Bcast(isSorted, 1, MPI_INT, root, comm);

    free(firsts);
    free(lasts);

    return MPI_SUCCESS;
}
//notes
double* merge(int n, double* a, int m, double* b) {
    int i, j, k;
    double* c = (double*)calloc(n + m, sizeof(double));

    for (i = j = k = 0; (i < n) && (j < m); )
    {
        if (a[i] <= b[j])
        {
            c[k++] = a[i++];
        }
        else
        {
            c[k++] = b[j++];
        }
    }
    if (i == n)
    {
        for (; j < m; )
        {
            c[k++] = b[j++];
        }
    }
    else
    {
        for (; i < n; )
        {
            c[k++] = a[i++];
        }
    }
    return c;
}

//notes
void merge_sort(int n, double* a) {
    double* c;
    int i;

    if (n <= 1)
    {
        return;
    }
    if (n == 2)
    {
        if (a[0] > a[1])
        {
            swap(&a[0], &a[1]);
        }
        return;
    }

    merge_sort(n / 2, a);
    merge_sort(n - n / 2, a + n / 2);
    c = merge(n / 2, a, n - n / 2, a + n / 2);
    for (i = 0; i < n; i++)
    {
        a[i] = c[i];
    }
    free(c);
}

//notes
void swap(double* a, double* b) {
    double temp;
    temp = *a;
    *a = *b;
    *b = temp;
}



//todo de pus in excel si varianta noua si cei doi speed up