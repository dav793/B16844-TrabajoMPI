/**
 * 
 *  Compile:    $HOME/opt/usr/local/bin/mpic++ -std=c++11 -o ej4 ./ej4.cpp
 *  Exec:       $HOME/opt/usr/local/bin/mpiexec -np 4 ./ej4 1024
 * 
 *  Contraints:
 *      n must be a power of 2
 *      number of processes must be a power of 2
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <math.h>
#include <algorithm>
#include <vector>

using namespace std;

void getArgs(char* argv[], int& n);
int getRandomNum();
void mergeSortStep(vector<int>& s1, vector<int>& s2, vector<int>& result);
void mergeSort(int* unsorted, int size, int* sorted);

std::random_device rd; // obtain a random number from hardware
std::mt19937 eng(rd()); // seed the generator
std::uniform_int_distribution<> distr(0, 1000); // define the range

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

    // generate array of unsorted numbers
    int* arr = NULL;
    if (rank == 0) {
        arr = new int[n];
        for (int i = 0; i < n; ++i) {
            arr[i] = getRandomNum();
        }

        // print unsorted array
        // printf("Unsorted:\n");
        // for (int i = 0; i < n; ++i) {
        //     printf("%d ", arr[i]);
        // }
        // printf("\n\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double local_start = MPI_Wtime();

    // give every process a chunk of the array
    int local_arr_size = n/proc_qty;
    int* local_arr = new int[local_arr_size];
    MPI_Scatter(arr, local_arr_size, MPI_INT, local_arr, local_arr_size, MPI_INT, 0, MPI_COMM_WORLD);
    
    // print local unsorted array
    // printf("proc %d unsorted array:\n\t", rank);
    // for (int i = 0; i < local_arr_size; ++i) {
    //     printf("%d ", local_arr[i]);
    // }
    // printf("\n");

    // sort local chunk
    int* local_arr_sorted = new int[local_arr_size];
    mergeSort(local_arr, local_arr_size, local_arr_sorted);

    int iterations = (int) log2(proc_qty); 
    for (int i = 0; i < iterations; ++i) {

        int multiplier = pow( 2, i+1 );
        int partner_proc;
        bool is_receiver;

        if (rank % multiplier == 0) {

            is_receiver = true;
            partner_proc = rank + pow( 2, i );

            int partner_arr[local_arr_size];
            MPI_Recv(partner_arr, local_arr_size, MPI_INT, partner_proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 

            vector<int> local_vec(local_arr_sorted, local_arr_sorted + local_arr_size);
            vector<int> partner_vec(partner_arr, partner_arr + local_arr_size);
            vector<int> merged_vec;
            mergeSortStep(local_vec, partner_vec, merged_vec);

            delete [] local_arr_sorted;
            local_arr_sorted = new int[local_arr_size * 2];
            copy(merged_vec.begin(), merged_vec.end(), local_arr_sorted);

        }
        else if (rank % multiplier == pow( 2, i )) {

            is_receiver = false;
            partner_proc = rank - pow( 2, i );

            MPI_Send(local_arr_sorted, local_arr_size, MPI_INT, partner_proc, 0, MPI_COMM_WORLD);

        }     

        local_arr_size *= 2;

    }

    double local_finish = MPI_Wtime();
    double local_elapsed = local_finish - local_start;
    double global_elapsed;
    MPI_Reduce(&local_elapsed, &global_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // print results
    if (rank == 0) {
        printf("Result:\n");
        for (int i = 0; i < n; ++i) {
            printf("%d ", local_arr_sorted[i]);
        }
        printf("\n\n");

        printf("Execution time: %f\n", global_elapsed);
    }

    // clean up
    if (rank == 0) {
        delete [] arr;
    }
    delete [] local_arr;
    delete [] local_arr_sorted;

	MPI_Barrier(MPI_COMM_WORLD); // para sincronizar la finalizaci�n de los procesos

    MPI_Finalize();
    return 0;
}

void getArgs(char* argv[], int& n) {
	n = strtol(argv[1], NULL, 10);
}

// generate random num
int getRandomNum() {
    return distr(eng);
}

// perform single step of merge sort algoritm (produce <result> whose size is <s1> + <s2>)
// s1 and s2 MUST be ordered
void mergeSortStep(vector<int>& s1, vector<int>& s2, vector<int>& result) {
    result.clear();
    result.resize(s1.size() + s2.size());
    merge(s1.begin(), s1.end(), s2.begin(), s2.end(), result.begin());
}

// perform complete merge sort algoritm
// <unsorted> is an unsorted array
// <sorted> is the produced sorted array
// <size> MUST be a power of 2
void mergeSort(int* unsorted, int size, int* sorted) {

    int sub_arr_size = size;
    int sub_arr_compl_size = 1;

    // setup initial vector of single elements
    vector<vector<int>> s;
    s.clear();
    s.resize(sub_arr_size);
    for (int i = 0; i < sub_arr_size; ++i) {
        s[i].clear();
        s[i].resize(sub_arr_compl_size);
        
        for (int j = 0; j < sub_arr_compl_size; ++j) {
            s[i][j] = unsorted[i+j];
        }
    }

    // int k = 0;
    while (sub_arr_size > 1) {  // merge sort iteratively

        sub_arr_size = sub_arr_size/2;
        sub_arr_compl_size = sub_arr_compl_size*2;

        vector<vector<int>> s2;
        s2.clear();
        s2.resize(sub_arr_size);
        for (int i = 0; i < sub_arr_size; ++i) {
            mergeSortStep(s[(i*2)], s[(i*2)+1], s2[i]);
        }

        s.clear();
        s.resize(sub_arr_size);
        copy(s2.begin(), s2.end(), s.begin());

        // k++;
        // printf("proc %d s%d:\n ", rank, k);
        // for (int i = 0; i < sub_arr_size; ++i) {
        //     printf("\t");    
        //     for (int j = 0; j < s[i].size(); ++j) {
        //         printf("%d, ", s[i][j]);    
        //     }
        //     printf("\n");
        // }
        // printf("\n");

    }

    // copy final vector into output array
    for (int i = 0; i < size; ++i) {
        sorted[i] = s[0][i];
    }

}