#include <iostream>
#include <mpi.h>

#define N 1500

const int MAIN_THREAD_NUM = 0, MATRIX1_TAG = 1, MATRIX2_TAG = 2, RESULT_ROW_TAG = 3;

using namespace std;


int main(int argc, char* argv[]) {
	int *matrix1 = new int[N * N], *matrix2 = new int[N * N], *result_matrix = new int[N * N];
	int *result_thread_row = new int[N];

	MPI_Init(&argc, &argv);

	int cur_thread, all_threads;
	MPI_Comm_rank(MPI_COMM_WORLD, &cur_thread);
	MPI_Comm_size(MPI_COMM_WORLD, &all_threads);

	double start_time = MPI_Wtime(), end_time;
	if (cur_thread == MAIN_THREAD_NUM) {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				matrix1[i * N + j] = rand() % N;
				matrix2[i * N + j] = rand() % N;
			}
		}

		for (int thread = 1; thread < all_threads; thread++) {
			MPI_Send(matrix1, N * N, MPI_INTEGER, thread, MATRIX1_TAG, MPI_COMM_WORLD);
			MPI_Send(matrix2, N * N, MPI_INTEGER, thread, MATRIX2_TAG, MPI_COMM_WORLD);
		}

		for (int row = 0; row < N; row += all_threads) {
			for (int col = 0; col < N; col++) {
				result_matrix[row * N + col] = 0;
				for (int k = 0; k < N; k++) {
					result_matrix[row * N + col] += matrix1[row * N + k] * matrix2[k * N + col];
				}
			}
		}

		for (int thread = 1; thread < all_threads; thread++) {
			for (int row = thread; row < N; row += all_threads) {
				MPI_Recv(result_thread_row, N, MPI_INTEGER, thread, row, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				for (int j = 0; j < N; j++) {
					result_matrix[row * N + j] = result_thread_row[j];
				}
			}
		}
	}
	else {
		MPI_Recv(matrix1, N * N, MPI_INTEGER, MAIN_THREAD_NUM, MATRIX1_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(matrix2, N * N, MPI_INTEGER, MAIN_THREAD_NUM, MATRIX2_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		for (int row = cur_thread; row < N; row += all_threads) {
			for (int col = 0; col < N; col++) {
				result_thread_row[col] = 0;
				for (int k = 0; k < N; k++) {
					result_thread_row[col] += matrix1[row * N + k] * matrix2[k * N + col];
				}
			}
			MPI_Send(result_thread_row, N, MPI_INTEGER, MAIN_THREAD_NUM, row, MPI_COMM_WORLD);
		}
	}

	end_time = MPI_Wtime();
	MPI_Finalize();

	if (cur_thread == 0) {

		cout << "Execution time in " << all_threads << " threads:\n" << end_time - start_time << endl;

		// по красоте чистим память для всех выделенных массивов
		delete[] matrix1, matrix2, result_matrix, result_thread_row;
	}
}
