#include <iostream>
#include <omp.h>

#define N 1500

using namespace std;


int main() {
	int **matrix1 = new int* [N], **matrix2 = new int* [N], **result_matrix = new int* [N];
	for (int i = 0; i < N; i++) {
		matrix1[i] = new int[N];
		matrix2[i] = new int[N];
		result_matrix[i] = new int[N];
		for (int j = 0; j < N; j++) {
			matrix1[i][j] = rand() % N;
			matrix2[i][j] = rand() % N;
		}
	}

	double start_time;
	int all_threads = omp_get_max_threads();
	for (int i = 1; i <= all_threads; i++) {
		start_time = omp_get_wtime();
		omp_set_num_threads(i);

		#pragma omp parallel for
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				result_matrix[i][j] = 0;
				for (int k = 0; k < N; k++) {
					result_matrix[i][j] += matrix1[i][k] * matrix2[k][j];
				}
			}
		}

		cout << "Execution time in " << i << " threads:\n" << omp_get_wtime() - start_time << endl;
	}

	// по красоте чистим память для всех выделенных массивов
	for (int i = 0; i < N; i++) {
		delete[] matrix1[i], matrix2[i], result_matrix[i];
	}
	delete[] matrix1, matrix2, result_matrix;
}
