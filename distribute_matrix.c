#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define global_height 15
#define global_width 15

void allocate_matrix(int ***matrix_ptr, int height, int width)
{
    int **matrix = (int **)malloc(sizeof(int *) * height);
    int *matrix_data = (int *)malloc(sizeof(int) * height * width);
    for (int i = 0; i < height; i++)
    {
        matrix[i] = &(matrix_data[i * width]);
    }
    *matrix_ptr = matrix;
}

void fill_matrix(int **matrix, int height, int width)
{
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            matrix[i][j] = j + i;
        }
    }
}

void calc_offset(int *send_from, int *size_to_send, int batch, int larger_batch, int nprocess)
{
    for (int i = 0; i < nprocess; i++)
    {
        send_from[i] = (i == 0) ? 0 : send_from[i - 1] + size_to_send[i - 1];
        size_to_send[i] = (i < larger_batch) ? batch + 1 : batch;
    }
    send_from[nprocess] = global_height;
}

int main(int argc, char *argv[])
{
    // Parámetros MPI
    int myrank;
    int nprocess;
    int root = 0;
    int tag = 0;

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocess);

    // Matriz de datos
    int **data_matrix;

    // Datos par alos procesos
    int batch = global_height / nprocess;
    int larger_batch = global_height % nprocess;
    int send_from[nprocess + 1];
    int size_to_send[nprocess];

    calc_offset(send_from, size_to_send, batch, larger_batch, nprocess);

    // Reserva de memoria
    if (myrank == root)
    {
        allocate_matrix(&data_matrix, global_height, global_width);
        fill_matrix(data_matrix, global_height, global_width);
        for (int i = 1; i < nprocess; i++)
        {
            for (int j = send_from[i]; j < send_from[i + 1]; j++)
            {
                MPI_Send(&data_matrix[j][0], global_width, MPI_INT, i, tag, MPI_COMM_WORLD);
            }
        }

        printf("Myrank: %d\n", myrank);
        for (int i = 0; i < global_height; i++)
        {
            for (int j = 0; j < global_width; j++)
            {
                printf("%d ", data_matrix[i][j]);
            }
            printf("\n");
        }
    }
    else
    {
        allocate_matrix(&data_matrix, size_to_send[myrank], global_width);
        for (int i = 0; i < size_to_send[myrank]; i++)
        {
            MPI_Recv(&data_matrix[i][0], global_width, MPI_INT, root, tag, MPI_COMM_WORLD, &status);
        }

        printf("Myrank: %d\n", myrank);
        for (int i = 0; i < size_to_send[myrank]; i++)
        {
            for (int j = 0; j < global_width; j++)
            {
                printf("%d ", data_matrix[i][j]);
            }
            printf("\n");
        }
    }

    MPI_Finalize();

    return 0;
}
