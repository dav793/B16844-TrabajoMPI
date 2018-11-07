/**
 * 
 *  Inputs:		
 *      1: number of intervals
 *      2: min value 
 *      3: max value 
 *      4: number of generated integers
 * 
 *  Prep:       export TMPDIR=/tmp
 *  Compile:    $HOME/opt/usr/local/bin/mpic++ -std=c++11 -o ej2 ./ej2.cpp
 *  Exec:       $HOME/opt/usr/local/bin/mpiexec -np 10 ./ej2 10 0 100000 10000
 * 
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <string>

void getArgs(char* argv[], int& num_intervals, int& min_value, int& max_value, int& num_ints);
void putInBin(int& num, int& num_bins, int& min_value, int& max_value, int* bins);
void printHistogram(int& num_bins, int& min_value, int& max_value, int& num_ints, int* bins);
void printHistogramNumbers(int& num_bins, int& min_value, int& max_value, int* bins);

int main(int argc, char* argv[]) {

    int rank;
    int proc_qty;
    MPI_Status mpi_status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_qty);

	int num_intervals;
    int min_value;
    int max_value;
    int num_ints;
	getArgs(argv, num_intervals, min_value, max_value, num_ints);

    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<> distr(min_value, max_value); // define the range

    /* ejecuci�n del proceso principal */
    if (num_ints % proc_qty != 0)
		printf("ERROR: Numero de enteros generados debe ser divisible por la cantidad de procesos.\n");
	else {
        int ints_per_proc = num_ints / proc_qty;
        int local_bins[num_intervals];
        for (int i = 0; i < num_intervals; ++i) {
            local_bins[i] = 0;
        }

        MPI_Barrier(MPI_COMM_WORLD);
		double local_start = MPI_Wtime();

        for (int i = 0; i < ints_per_proc; ++i) {
            int num = distr(eng);
            putInBin(num, num_intervals, min_value, max_value, local_bins);
        }

        // printf("proc %d local bins\n", rank);
        // printHistogram(num_intervals, min_value, max_value, num_ints, local_bins);
        // printf("\n");

        int* all_bins = NULL;
        if (rank == 0) {
            all_bins = new int[num_intervals * proc_qty];
        }
        MPI_Gather(local_bins, num_intervals, MPI_INT, all_bins, num_intervals, MPI_INT, 0, MPI_COMM_WORLD);

        double local_finish = MPI_Wtime();
		double local_elapsed = local_finish - local_start;
		double global_elapsed;
		MPI_Reduce(&local_elapsed, &global_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            int global_bins[num_intervals];
            for (int i = 0; i < num_intervals; ++i) {
                global_bins[i] = 0;
                for (int j = 0; j < proc_qty; ++j) {
                    global_bins[i] += all_bins[i + j*num_intervals];
                }
            }; 
            
            printHistogram(num_intervals, min_value, max_value, num_ints, global_bins);
            printf("\n");

            // printHistogramNumbers(num_intervals, min_value, max_value, global_bins);
            // printf("\n");

            printf("Execution time: %f\n", global_elapsed);
        }

    }
	
	MPI_Barrier(MPI_COMM_WORLD); // para sincronizar la finalizaci�n de los procesos

    MPI_Finalize();
    return 0;
}

void getArgs(char* argv[], int& num_intervals, int& min_value, int& max_value, int& num_ints) {
	num_intervals = strtol(argv[1], NULL, 10);
    min_value = strtol(argv[2], NULL, 10);
    max_value = strtol(argv[3], NULL, 10);
    num_ints = strtol(argv[4], NULL, 10);
}

void putInBin(int& num, int& num_bins, int& min_value, int& max_value, int* bins) {
    int interval_size = (max_value - min_value) / num_bins;

    for (int i = 0; i < num_bins; ++i) {
        if ( num >= min_value + (interval_size*i) && num <= min_value + (interval_size*(i+1)) ) {
            bins[i]++;
            break;
        }
    }
}

void printHistogram(int& num_bins, int& min_value, int& max_value, int& num_ints, int* bins) {
    int interval_size = (max_value - min_value) / num_bins;

    for (int i = 0; i < num_bins; ++i) {
        std::string label = "" + std::to_string(min_value + (interval_size*i)) + " - " + std::to_string(min_value + (interval_size*(i+1))) + ": ";
        printf("%-20s", label.c_str());
        int num_xs = bins[i] / (num_ints/1000);
        for (int i = 0; i < num_xs; ++i) {
            printf("X");
        }
        printf("\n");
    }    
}

void printHistogramNumbers(int& num_bins, int& min_value, int& max_value, int* bins) {
    int interval_size = (max_value - min_value) / num_bins;

    for (int i = 0; i < num_bins; ++i) {
        std::string label = "" + std::to_string(min_value + (interval_size*i)) + " - " + std::to_string(min_value + (interval_size*(i+1))) + ": ";
        printf("%-20s%d\n", label.c_str(), bins[i]);
    }    
}