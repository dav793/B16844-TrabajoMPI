/**
 * 
 * 	Inputs:		
 * 		1: number of tosses
 * 
 * 	Prep:       export TMPDIR=/tmp
 *  Compile:    $HOME/opt/usr/local/bin/mpic++ -std=c++11 -o ej1 ./ej1.cpp
 *  Exec:       $HOME/opt/usr/local/bin/mpiexec -np 10 ./ej1 1000000
 * 
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <math.h>

void getArgs(char* argv[], int& num_tosses);
float getRandomNum();

std::random_device rd; // obtain a random number from hardware
std::mt19937 eng(rd()); // seed the generator
std::uniform_int_distribution<> distr(0, 100); // define the range

int main(int argc, char* argv[]) {

    int rank;
    int proc_qty;
    MPI_Status mpi_status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_qty);

	int num_tosses_per_proc;
	int num_tosses;
	getArgs(argv, num_tosses);

	if (num_tosses % proc_qty != 0)
		printf("ERROR: Numero de tiros debe ser divisible por la cantidad de procesos.\n");
	else {
		int local_num_in_circle = 0;
		num_tosses_per_proc = num_tosses / proc_qty;

		MPI_Barrier(MPI_COMM_WORLD);
		double local_start = MPI_Wtime();

		for (int i = 0; i < num_tosses_per_proc; ++i) {
			float x = getRandomNum();
			float y = getRandomNum();
			if ( x*x + y*y <= 1 )
				local_num_in_circle++;
		}

		int global_num_in_circle;
    	MPI_Reduce(&local_num_in_circle, &global_num_in_circle, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

		double local_finish = MPI_Wtime();
		double local_elapsed = local_finish - local_start;
		double global_elapsed;
		MPI_Reduce(&local_elapsed, &global_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

		if (rank == 0) {
			float pi_estimate = (4*global_num_in_circle) / ((float) num_tosses);
			printf("PI estimate: %f\n", pi_estimate);
			printf("Precision error: %f\n", abs(4.0*atan(1.0) - pi_estimate));
			printf("Execution time: %f\n", global_elapsed);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}

void getArgs(char* argv[], int& num_tosses) {
	num_tosses = strtol(argv[1], NULL, 10);
}

// generate random float between -1 y 1
float getRandomNum() {
    return ((float) distr(eng)) / 50 - 1;
}