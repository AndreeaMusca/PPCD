#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#pragma comment(lib, "C:\\Program Files (x86)\\Microsoft SDKs\\MPI\\Lib\\x64\\msmpi.lib")
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int testWinner(char player1, char player2);
const char* rsp = "rsp";

int main(int argc, char** argv)
{
	int currentPlay, otherPlay;
	int rank, size, winner, tag1 = 1, tag2 = 2;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	srand(rank + time(0));

	currentPlay = rand() % 3;
	if (rank == 0)
	{
		MPI_Send(&currentPlay, 1, MPI_INT, 1, tag1, MPI_COMM_WORLD);//ce trimit aici
		MPI_Recv(&otherPlay, 1, MPI_INT, 1, tag2, MPI_COMM_WORLD, &status);
		winner = testWinner(rsp[currentPlay], rsp[otherPlay]);
	}
	else if (rank == 1) //proesor curent 1
	{
		MPI_Recv(&otherPlay, 1, MPI_INT, 0, tag1, MPI_COMM_WORLD, &status);//primesc aici
		MPI_Send(&currentPlay, 1, MPI_INT, 0, tag2, MPI_COMM_WORLD);
		winner = testWinner(rsp[otherPlay], rsp[currentPlay]);
	}

	if (winner == -1)
	{
		printf("Draw\n");
	}
	else {
		printf("the winner from processor %d is %d\n",rank, winner);
	}
	

	MPI_Finalize();
	return 0;
}

int testWinner(char player1, char player2)
{
	if ((player1 == 'r' && player2 == 's') ||
		(player1 == 's' && player2 == 'p') ||
		(player1 == 'p' && player2 == 'r'))
		return 0;
	if (player1 == player2)
		return -1;
	return 1;
}
