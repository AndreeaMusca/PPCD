#include <stdio.h>
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#pragma comment(lib, "C:\\Program Files (x86)\\Microsoft SDKs\\MPI\\Lib\\x64\\msmpi.lib")


#include <stdio.h>
#include <math.h>

int isprime(long n);

int main(int argc, char** argv) {

    int n = 10000000, i, size, rank, count = 0;
    double time;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    time = MPI_Wtime();
    for (i = rank * n / size + 1; i < (rank + 1) * n / size; i++)
        if (isprime(i)) count++;
    time = MPI_Wtime() - time;

    printf("Processor %d finds %d primes in %lf\n", rank, count, time);



    MPI_Finalize();
}


int isprime(long n) {
    int d, ans = 1;
    if (n == 1)return 0;
    if (n == 2 || n == 3)return 1;
    if (n % 2 == 0)return 0;
    for (d = 3; d <= sqrt(n); d += 2)
        if (n % d == 0) { ans = 0; break; }
    return ans;
}


//int main(int  argc, char* argv[])
//{
//    int rank, size;
//    int i, namelen;
//    char processor_name[MPI_MAX_PROCESSOR_NAME];
//
//
//    MPI_Init(&argc, &argv);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//    //MPI_Get_processor_name(processor_name,&namelen);
//
//    //printf("called on %s\n",processor_name);
//
//    printf("Hello world from process %d of %d\n", rank, size);
//
//    MPI_Finalize();
//    return 0;
//}
 