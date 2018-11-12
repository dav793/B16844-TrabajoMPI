/**
 * 
 *  Compile:    $HOME/opt/usr/local/bin/mpic++ -std=c++11 -o ej3 ./ej3.cpp
 *  Exec:       $HOME/opt/usr/local/bin/mpiexec -np 2 ./ej3 10
 *
 *  Contraints:
 *      n should be divisible by number of processes
 * 
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

void getArgs(char* argv[], int& num_tosses);
void findPrimesInInterval(int start, int interval_size, int* primes);
void findTerms(int n, int* primes, int& term_a, int& term_b, int& term_c);

int main(int argc, char* argv[]) {

    int rank;
    int proc_qty;
    MPI_Status mpi_status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_qty);

	int n;
	getArgs(argv, n);

    if (n % proc_qty != 0) {
        printf("ERROR: n must be divisible by process count\n");
        return -1;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double local_start = MPI_Wtime();

    // find all primes up to n
    int local_interval_size = n/proc_qty;
    int local_primes_start = rank * local_interval_size;
    int* local_primes = new int[local_interval_size];
    findPrimesInInterval(local_primes_start, local_interval_size, local_primes);

    // gather local primes into global primes
    int* global_primes = NULL;
    if (rank == 0)
        global_primes = new int[n];
    MPI_Gather(local_primes, local_interval_size, MPI_INT, global_primes, local_interval_size, MPI_INT, 0, MPI_COMM_WORLD);

    // merge global primes into array of primes
    int* primes = new int[n];
    if (rank == 0) {
        int currentPrime = 0;
        for (int i = 0; i < n; ++i) {
            if (global_primes[i] > 0) {
                primes[currentPrime] = global_primes[i];
                currentPrime++;
            }
        }
        for (int i = currentPrime; i < n; ++i) {
            primes[i] = -1;
        }
    }

    // broadcast primes to every process
    MPI_Bcast(primes, n, MPI_INT, 0, MPI_COMM_WORLD);

    // find sums for local interval
    int local_min_n = rank * local_interval_size;
    int local_max_n = (rank+1) * local_interval_size - 1;
    // printf("process %d gets from %d to %d\n", rank, local_min_n+1, local_max_n+1);

    int* local_interval = new int[local_interval_size*3];

    for (int i = 0; i < local_interval_size*3; i+=3) {
        int local_n = rank * local_interval_size + i/3 + 1;
        // printf("proc %d: find terms of %d\n", rank, local_n);

        int a, b, c;
        findTerms(local_n, primes, a, b, c);
        local_interval[i] = a;
        local_interval[i+1] = b;
        local_interval[i+2] = c;
    }

    // gather local results into global result
    int* global_result = NULL;
    if (rank == 0)
        global_result = new int[n*3];
    MPI_Gather(local_interval, local_interval_size*3, MPI_INT, global_result, local_interval_size*3, MPI_INT, 0, MPI_COMM_WORLD);

    double local_finish = MPI_Wtime();
    double local_elapsed = local_finish - local_start;
    double global_elapsed;
    MPI_Reduce(&local_elapsed, &global_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // print results
    if (rank == 0) {
        for (int i = 0; i < n*3; i+=3) {

            int num = i/3 + 1;
            int term_a = global_result[i];
            int term_b = global_result[i+1];
            int term_c = global_result[i+2];

            if (term_a < 0 && term_b < 0 && term_c < 0) {
                printf("%d = N/A\n", num);
            }
            else {
                if (term_c < 0) {
                    printf("%d = %d + %d\n", num, term_a, term_b);
                }
                else {
                    printf("%d = %d + %d + %d\n", num, term_a, term_b, term_c);
                }
            }
        }

        printf("\nExecution time: %f\n", global_elapsed);

        // clean up
        delete [] global_result;
        delete [] global_primes;
    }

    // clean up
    delete [] local_primes;
    delete [] local_interval;
    delete [] primes;

	MPI_Barrier(MPI_COMM_WORLD); // para sincronizar la finalizaci�n de los procesos

    MPI_Finalize();
    return 0;
}

void getArgs(char* argv[], int& n) {
	n = strtol(argv[1], NULL, 10);
}

// fill <primes> with all primes from <start> up to but excluding <finish>
// fill the rest of <primes> with -1's
void findPrimesInInterval(int start, int interval_size, int* primes) {
    int end = start + interval_size;
    int length = 0;

    for (int i = start; i < end; ++i) {
        if (i == 0) {
            primes[0] = 1;
            length++;
            continue;
        }

        int num_factors = 0;
        for (int j = i; j > 0; j--) {
            if (i%j == 0) {
                num_factors++;
                if (num_factors > 2)
                    break;
            }
        }

        if (num_factors == 2) {
            primes[length] = i;
            length++;
        }
    }

    for (int i = length; i < interval_size; ++i) {
        primes[length] = -1;
        length++;
    }
}

/**
 * find the prime terms that add up to n
 * 
 * @param n : input - number to find the terms of
 * @param primes : input - array of all primes (up to at least n)
 * @param term_a : output - first term (or -1)
 * @param term_b : output - second term (or -1)
 * @param term_c : output - third term (or -1)
 */ 
void findTerms(int n, int* primes, int& term_a, int& term_b, int& term_c) {

    bool stop = false;
    term_a = -1;
    term_b = -1;
    term_c = -1;

    for (int a = 1; a <= n; ++a) {

        int prime_a = primes[a];
        if (prime_a >= n || prime_a <= 0 || stop)
            break;

        for (int b = 1; b <= n; ++b) {

            int prime_b = primes[b];
            if (prime_b >= n || prime_b <= 0 || stop)
                break;

            for (int c = 1; c <= n; ++c) {

                int prime_c = primes[c];
                if (prime_c >= n || prime_c <= 0 || stop)
                    break;
                
                if ( prime_a + prime_b == n) {
                    stop = true;
                    term_a = prime_a;
                    term_b = prime_b;
                    a = n+1;
                    b = n+1;
                    c = n+1;
                    break;
                }
                else if ( prime_a + prime_b + prime_c == n ) {
                    stop = true;
                    term_a = prime_a;
                    term_b = prime_b;
                    term_c = prime_c;
                    a = n+1;
                    b = n+1;
                    c = n+1;
                    break;
                }

            }
        }    
    }

    // printf("%d + %d + %d\n", term_a, term_b, term_c);
}