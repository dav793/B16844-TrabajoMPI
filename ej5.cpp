
/**
 * 
 *  Compile:    $HOME/opt/usr/local/bin/mpic++ -std=c++11 -o ej5 ./ej5.cpp
 *  Exec:       $HOME/opt/usr/local/bin/mpiexec -np 2 ./ej5 10
 * 
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip> // std::setw, setprecision

using namespace std;

void leeAdyacencias(vector< vector<string> > tokens, vector< vector< int > >& ma, int& cntVertices);
void readFile(string filename, vector<string>& lines);
void tokenize(string& line, vector<string>& tokens);
void algoritmoFloydWarshall(const vector< vector< int > >& ma, int cntVertices, vector< vector< int > >& mc, vector< vector<int> >& next);
void obtenerCaminoMasCorto(int v1, int v2, vector< vector<int> >& next, vector<int>& camino);

int main(int argc, char* argv[]) {

    int rank;
    int proc_qty;
    MPI_Status mpi_status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_qty);

	string nombreArchivoEntrada = "gpequenyo.txt"; // formato "*.txt"
	vector< vector< int > > matrizAdyacencias;
	vector< vector< int > > matrizCostos;
    vector< vector< int > > caminosMasCortos;
	int cntVertices;

    vector<string> lines;
    readFile(nombreArchivoEntrada, lines);

    vector< vector<string> > tokens;
    for (int i = 0; i < lines.size(); ++i) {
        vector<string> temp;
        tokenize(lines[i], temp);
        tokens.push_back(temp);
    }

    // se cargan las adyacencias del archivo
    leeAdyacencias(tokens, matrizAdyacencias, cntVertices);

    if (rank == 0) {

        // imprimir matriz de adyacencias
        leeAdyacencias(tokens, matrizAdyacencias, cntVertices);
        // cout << "MATRIZ DE ADYACENCIAS" << endl;
        // for (int i = 0; i < cntVertices; i++) {
        //     for (int j = 0; j < cntVertices; j++)
        //         cout << (matrizAdyacencias[i][j] == INT_MAX ? 0 : matrizAdyacencias[i][j]) << ",";
        //     cout << endl;
        // }

        // generar matriz de costos y la de caminos
	    algoritmoFloydWarshall(matrizAdyacencias, cntVertices, matrizCostos, caminosMasCortos);

        // imprimir matriz de costos:
        // cout << endl << "MATRIZ DE COSTOS" << endl;
        // for (int i = 0; i < cntVertices; i++) {
        //     for (int j = 0; j < cntVertices; j++)
        //         cout << matrizCostos[i][j] << ',';
        //     cout << endl;
        // }

        // // imprimir matriz de caminos mas cortos:
        // cout << endl << "MATRIZ DE CAMINOS MAS CORTOS" << endl;
        // for (int i = 0; i < cntVertices; i++) {
        //     for (int j = 0; j < cntVertices; j++)
        //         cout << caminosMasCortos[i][j] << ',';
        //     cout << endl;
        // }

        // imprimir caminos mas cortos
        cout << endl << "CAMINOS MAS CORTOS" << endl;
        // cout << setw(10) << "origen " << setw(10) << "dest " << setw(10) << "dist " << setw(20) << "camino " << endl;
        for (int i = 0; i < cntVertices; i++) {
            for (int j = 0; j < cntVertices; j++) {
                if (i != j) {

                    // cout << setw(5) << i << setw(5) << j << setw(5) << matrizCostos[i][j] << setw(10);
                    cout << i << " -> " << j << "\t:\tdist: " << matrizCostos[i][j] << "\t\t ruta: ";

                    vector<int> camino;
                    obtenerCaminoMasCorto(i, j, caminosMasCortos, camino);
                    for (int k = 0; k < camino.size(); ++k) {
                        if (k > 0)
                            cout << "-> ";
                        cout << camino[k] << " ";
                    }
                    cout << endl;

                }
            }
        }

    }



	MPI_Barrier(MPI_COMM_WORLD); // para sincronizar la finalizaci�n de los procesos

    MPI_Finalize();
    return 0;
}

void leeAdyacencias(vector< vector<string> > tokens, vector< vector< int > >& ma, int& cntVertices) {

    cntVertices = stoi(tokens[0][0]);  // el primer n�mero del archivo es la cantidad de v�rtices

    vector< int > v;
	v.resize(cntVertices, INT_MAX);
	ma.resize(cntVertices, v);

    for (int i = 1; i < tokens.size(); ++i) {
        for (int j = 0; j < tokens[i].size(); ++j) {
            int ady = stoi(tokens[i][j]);
            ma[i-1][ady] = 1;
        }
    }

}

void tokenize(string& line, vector<string>& tokens) {
    std::stringstream ss(line);
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> vstrings(begin, end);
    tokens.resize(vstrings.size());
    std::copy(vstrings.begin(), vstrings.end(), tokens.begin());
}

void readFile(string filename, vector<string>& lines) {
    std::ifstream f (filename);
    string line;

    while(getline(f, line)) {
        lines.push_back(line);
    }
}

void algoritmoFloydWarshall(const vector< vector< int > >& ma, int cntVertices, vector< vector< int > >& mc, vector< vector<int> >& next) {
	
    vector< int > v;
	v.resize(cntVertices, INT_MAX);
	mc.resize(cntVertices, v);

    vector< int > w;
    w.resize(cntVertices, -1);
    next.resize(cntVertices, w);

    for (int i = 0; i < cntVertices; ++i) {
        for (int j = 0; j < cntVertices; ++j) {
            mc[i][j] = ma[i][j];
            next[i][j] = j;
        }
    }

    for (int k = 0; k < cntVertices; ++k) {
        for (int i = 0; i < cntVertices; ++i) {
            for (int j = 0; j < cntVertices; ++j) {
                if (mc[i][j] > mc[i][k] + mc[k][j] && mc[i][k] < INT_MAX && mc[k][j] < INT_MAX) {
                    mc[i][j] = mc[i][k] + mc[k][j];
                    next[i][j] = next[i][k];
                }
            }
        }  
    }

}

void obtenerCaminoMasCorto(int v1, int v2, vector< vector<int> >& next, vector<int>& camino) {
    camino.clear();
    if (next[v1][v2] == -1)
        return;

    camino.push_back(v1);
    while (v1 != v2) {
        v1 = next[v1][v2];
        camino.push_back(v1);
    }
}