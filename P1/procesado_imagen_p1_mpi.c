#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

// Reserva de memoria
unsigned char **allocate_matrix(int, int);

// Lectura de fichero binario
void read_binary_data(unsigned char **, const char *, int, int);

// Volcado de información en archivo binario
void write_binary_data(unsigned char **, const char *, int, int);

// Añade padding a cualquier matrix NxM o NxN
void padding(unsigned char **, int, int);

// Calcula el offset para las matrices que almacenan la imagen con padding
void calc_pad_offset(int *, int *, int, int, int);

// Ajusta las filas necesarias para cada matriz teniendo en cuenta el padding
void adjust_rows(int *, int *, int, int);

// Calcula el offset para las matrices resultantes
void calc_result_offset(int *, int *, int *, int);

// Algoritmo para ordenar un vector
void bubble_sort(int[], int);

// Utilidad para imprimir matrices
void print_matrix(unsigned char **, int, int);

// Filtros
void mean_filter(unsigned char **, unsigned char **, int, int);

void median_filter(unsigned char **, unsigned char **, int, int);

void sobel_filter(unsigned char **, unsigned char **, int, int);

// TODO

// witdth, height, [1,2,3 {media, mediana, sobel}], input file
int main(int argc, char *argv[])
{
    int selected_filter = atoi(argv[3]);

    // Archivos input/output
    const char *input_file = argv[4], *output_file = "output.raw";

    // Matriz de datos y imagen resultado
    unsigned char **matrix, **result;

    // Atributos de las imágenes
    int width = atoi(argv[1]), height = atoi(argv[2]);

    // Parámetros MPI
    int nprocess;
    int myrank;
    int tag = 0;
    int root = 0;
    double start, end;

    // Inicio MPI
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocess);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // Parámetros para las matrices con padding
    int pad_batch;         // Filas por procesos
    int pad_bigger_batch;  // Cantidad de Procesos que tendran una fila más
    int *pad_send_from;    // Índices de filas iniciales
    int *pad_size_to_send; // Desplazamiento desde filas iniciales (offset)

    // Parámetros para las matrices resultantes
    int *res_send_from;
    int *res_size_to_send;

    // Tipos de datos personalizados para utilizar con scatterv y gatherv
    MPI_Datatype pad_row_type; // Tipo de dato para matrices con padding
    MPI_Type_contiguous(width + 2, MPI_UNSIGNED_CHAR, &pad_row_type);
    MPI_Type_commit(&pad_row_type);

    MPI_Datatype row_type; // Tipo de dato para matrices resultantes
    MPI_Type_contiguous(width, MPI_UNSIGNED_CHAR, &row_type);
    MPI_Type_commit(&row_type);

    // Carga de datos en root
    if (myrank == root)
    {
        matrix = allocate_matrix(width + 2, height + 2);
        result = allocate_matrix(width, height);
        read_binary_data(matrix, input_file, width, height);
        padding(matrix, width, height);

        // Si el número de procesos es 1, se procesa directamente en root
        if (nprocess == 1)
        {
            start = MPI_Wtime();
            sobel_filter(matrix, result, width, height);
            end = MPI_Wtime();
            printf("Runtime = %f\n", end - start);
            write_binary_data(result, output_file, width, height);
            MPI_Finalize();
            return 0;
        }
    }

    // Inicio distribución de los datos
    pad_batch = (height + 2) / nprocess;
    pad_bigger_batch = (height + 2) % nprocess;

    pad_send_from = malloc((height + 2) * sizeof(int));
    pad_size_to_send = malloc((height + 2) * sizeof(int));
    int *pad_rows = malloc(sizeof(int) * nprocess);

    calc_pad_offset(pad_send_from, pad_size_to_send, pad_batch, pad_bigger_batch, height + 2);
    adjust_rows(pad_size_to_send, pad_rows, myrank, nprocess);

    // Matrices locales para cada proceso
    unsigned char **local_matrix = allocate_matrix(width + 2, pad_rows[myrank]);

    // El condicional lo utilizo exclusivamente porque utilizo mi ordenador personal para probar código y sin este no puedo compilar el programa
    if (myrank == root)
    {
        MPI_Scatterv(&matrix[0][0], pad_size_to_send, pad_send_from, pad_row_type, &local_matrix[0][0], pad_size_to_send[myrank], pad_row_type, 0, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Scatterv(NULL, pad_size_to_send, pad_send_from, pad_row_type, &local_matrix[1][0], pad_size_to_send[myrank], pad_row_type, 0, MPI_COMM_WORLD);
    }

    // Se envían las filas extras a los procesos 2, 3 y 4
    if (nprocess == 4)
    {
        if (myrank == 0)
        {
            MPI_Send(&local_matrix[pad_rows[myrank] - 2][0], 1, pad_row_type, 1, tag, MPI_COMM_WORLD);
            MPI_Recv(&local_matrix[pad_rows[myrank] - 1][0], 1, pad_row_type, 1, tag, MPI_COMM_WORLD, &status);
        }
        else if (myrank == 1 || myrank == 2)
        {
            MPI_Send(&local_matrix[1][0], 1, pad_row_type, myrank - 1, tag, MPI_COMM_WORLD);
            MPI_Recv(&local_matrix[0][0], 1, pad_row_type, myrank - 1, tag, MPI_COMM_WORLD, &status);

            MPI_Send(&local_matrix[pad_rows[myrank] - 2][0], 1, pad_row_type, myrank + 1, tag, MPI_COMM_WORLD);
            MPI_Recv(&local_matrix[pad_rows[myrank] - 1][0], 1, pad_row_type, myrank + 1, tag, MPI_COMM_WORLD, &status);
        }
        else if (myrank == 3)
        {
            MPI_Send(&local_matrix[1][0], 1, pad_row_type, 2, tag, MPI_COMM_WORLD);
            MPI_Recv(&local_matrix[0][0], 1, pad_row_type, 2, tag, MPI_COMM_WORLD, &status);
        }
    }

    if (nprocess == 2)
    {
        if (myrank == 0)
        {
            MPI_Send(&local_matrix[pad_rows[myrank] - 2][0], 1, pad_row_type, 1, tag, MPI_COMM_WORLD);
            MPI_Recv(&local_matrix[pad_rows[myrank] - 1][0], 1, pad_row_type, 1, tag, MPI_COMM_WORLD, &status);
        }
        if (myrank == 1)
        {
            MPI_Send(&local_matrix[1][0], 1, pad_row_type, 0, tag, MPI_COMM_WORLD);
            MPI_Recv(&local_matrix[0][0], 1, pad_row_type, 0, tag, MPI_COMM_WORLD, &status);
        }
    }
    if (nprocess == 3)
    {
        if (myrank == 0)
        {
            MPI_Recv(&local_matrix[pad_rows[myrank] - 1][0], 1, pad_row_type, 1, tag, MPI_COMM_WORLD, &status);
            MPI_Send(&local_matrix[pad_rows[myrank] - 2][0], 1, pad_row_type, 1, tag, MPI_COMM_WORLD);
        }
        else if (myrank == 1)
        {
            MPI_Send(&local_matrix[1][0], 1, pad_row_type, 0, tag, MPI_COMM_WORLD);
            MPI_Send(&local_matrix[pad_rows[myrank] - 2][0], 1, pad_row_type, 2, tag, MPI_COMM_WORLD);

            MPI_Recv(&local_matrix[0][0], 1, pad_row_type, 0, tag, MPI_COMM_WORLD, &status);
            MPI_Recv(&local_matrix[pad_rows[myrank] - 1][0], 1, pad_row_type, 2, tag, MPI_COMM_WORLD, &status);
        }
        else if (myrank == 2)
        {
            MPI_Recv(&local_matrix[0][0], 1, pad_row_type, 1, tag, MPI_COMM_WORLD, &status);
            MPI_Send(&local_matrix[1][0], 1, pad_row_type, 1, tag, MPI_COMM_WORLD);
        }
    }

    // Fin Inicio distribución de los datos
    res_send_from = (int *)malloc(sizeof(int) * nprocess);
    res_size_to_send = (int *)malloc(sizeof(int) * nprocess);

    // Calcula el offset y filas para el scatterv
    calc_result_offset(res_size_to_send, pad_size_to_send, res_send_from, nprocess);

    // Se almacena memoria para los resultados locales
    unsigned char **local_result = allocate_matrix(width, res_size_to_send[myrank]);

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    // Convolute
    switch (selected_filter)
    {
    case 1:
        mean_filter(local_matrix, local_result, width, res_size_to_send[myrank]);
        break;

    case 2:
        median_filter(local_matrix, local_result, width, res_size_to_send[myrank]);
        break;

    case 3:
        sobel_filter(local_matrix, local_result, width, res_size_to_send[myrank]);
        break;

    default:
        break;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();

    //Recolecta la información de todos los resultados locales y los envía al proceso 0, el condicional también depende de en qué equipo se ejecute
    if (myrank == root)
    {
        MPI_Gatherv(&local_result[0][0], res_size_to_send[myrank], row_type, &result[0][0], res_size_to_send, res_send_from, row_type, 0, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Gatherv(&local_result[0][0], res_size_to_send[myrank], row_type, NULL, res_size_to_send, res_send_from, row_type, 0, MPI_COMM_WORLD);
    }

    if (myrank == 0)
    {
        write_binary_data(result, output_file, width, height);
    }

    MPI_Finalize();

    if (myrank == 0)
    {
        printf("Tiempo de ejecucion en s: %f\n", end - start);
    }
    return 0;
}

// Definición de funciones

void print_matrix(unsigned char **matrix, int width, int height)
{
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

void adjust_rows(int *size_to_send, int *new_rows, int myrank, int nprocess)
{
    if (nprocess == 4)
    {
        if (myrank == 0 || myrank == 3)
        {
            new_rows[myrank] = size_to_send[myrank] + 1;
        }
        else if (myrank == 1 || myrank == 2)
        {
            new_rows[myrank] = size_to_send[myrank] + 2;
        }
    }

    if (nprocess == 3)
    {
        if (myrank == 0 || myrank == 2)
        {
            new_rows[myrank] = size_to_send[myrank] + 1;
        }
        else if (myrank == 1)
        {
            new_rows[myrank] = size_to_send[myrank] + 2;
        }
    }

    if (nprocess == 2)
    {
        new_rows[myrank] = size_to_send[myrank] + 1;
    }
}
void calc_pad_offset(int *send_from, int *size_to_send, int batch, int bigger_batch, int height)
{
    for (int i = 0; i < height; i++)
    {
        send_from[i] = (i == 0) ? 0 : send_from[i - 1] + size_to_send[i - 1];
        size_to_send[i] = (i < bigger_batch) ? (batch + 1) : batch;
    }
}

void calc_result_offset(int *a_sts, int *f_sts, int *sf, int n)
{
    if (n == 4)
    {
        for (int i = 0; i < n; i++)
        {
            a_sts[i] = f_sts[i];
            if (i == 0 || i == n - 1)
                a_sts[i] += 1;
            if (i == 1 || i == 2)
                a_sts[i] += 2;
            a_sts[i] -= 2;
            sf[i] = (i == 0) ? 0 : a_sts[i - 1] + sf[i - 1];
        }
    }

    if (n == 3)
    {
        for (int i = 0; i < n; i++)
        {
            a_sts[i] = f_sts[i];
            if (i == 0 || i == 2)
                a_sts[i] += 1;
            if (i == 1)
                a_sts[i] += 2;
            a_sts[i] -= 2;
            sf[i] = (i == 0) ? 0 : a_sts[i - 1] + sf[i - 1];
        }
    }

    if (n == 2)
    {
        for (int i = 0; i < n; i++)
        {
            a_sts[i] = f_sts[i];
            a_sts[i] -= 2;
            sf[i] = (i == 0) ? 0 : a_sts[i - 1] + sf[i - 1];
        }
    }
}

unsigned char **allocate_matrix(int width, int height)
{
    unsigned char **matrix, *data;
    matrix = (unsigned char **)malloc(sizeof(unsigned char *) * height);
    data = (unsigned char *)malloc(sizeof(unsigned char) * width * height);
    for (int i = 0; i < height; i++)
    {
        matrix[i] = &(data[i * width]);
    }
    return matrix;
}

void read_binary_data(unsigned char **matrix, const char *input_file, int width, int height)
{
    FILE *raw_data;

    raw_data = fopen(input_file, "rb");

    for (int i = 1; i < height + 1; i++)
    {
        fread(matrix[i] + 1, sizeof(unsigned char), width, raw_data);
    }

    fclose(raw_data);
}

void write_binary_data(unsigned char **matrix, const char *output_file, int width, int height)
{
    FILE *output;

    output = fopen(output_file, "wb");
    for (int i = 0; i < height; i++)
    {
        fwrite(matrix[i], sizeof(unsigned char), width, output);
    }
    fclose(output);
}

void padding(unsigned char **matrix, int width, int height)
{
    matrix[0][0] = matrix[2][2];

    matrix[0][width + 1] = matrix[2][width - 1];

    matrix[height + 1][0] = matrix[height - 1][2];

    matrix[height + 1][width + 1] = matrix[height - 1][width - 1];

    for (int i = 1; i <= width; i++)
    {
        matrix[0][i] = matrix[2][i];
    }

    for (int i = 1; i <= width; i++)
    {
        matrix[height + 1][i] = matrix[height - 1][i];
    }

    for (int i = 1; i <= height; i++)
    {
        matrix[i][0] = matrix[i][2];
    }

    for (int i = 0; i <= height; i++)
    {
        matrix[i][width + 1] = matrix[i][width - 1];
    }
}

void mean_filter(unsigned char **matrix, unsigned char **result, int width, int height)
{
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            result[i][j] = (matrix[i][j] +
                            matrix[i][j + 1] +
                            matrix[i][j + 2] +
                            matrix[i + 1][j] +
                            matrix[i + 1][j + 1] +
                            matrix[i + 1][j + 2] +
                            matrix[i + 2][j] +
                            matrix[i + 2][j + 1] +
                            matrix[i + 2][j + 2]) /
                           9;
        }
    }
}

void median_filter(unsigned char **matrix, unsigned char **result, int width, int height)
{
    int data[9];

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            data[0] = matrix[i][j];
            data[1] = matrix[i][j + 1];
            data[2] = matrix[i][j + 2];
            data[3] = matrix[i + 1][j];
            data[4] = matrix[i + 1][j + 1];
            data[5] = matrix[i + 1][j + 2];
            data[6] = matrix[i + 2][j];
            data[7] = matrix[i + 2][j + 1];
            data[8] = matrix[i + 2][j + 2];

            // Se ordena el vector con los 9 datos y se escoge el del medio para el bit procesado
            bubble_sort(data, 9);

            result[i][j] = data[4];
        }
    }
}

void bubble_sort(int array[], int size)
{
    for (int i = 0; i < size - 1; i++)
    {
        for (int j = 0; j < size - i - 1; j++)
        {
            if (array[j] > array[i + 1])
            {
                int aux = array[j];
                array[j] = array[i];
                array[i] = aux;
            }
        }
    }
}

void sobel_filter(unsigned char **matrix, unsigned char **result, int width, int height)
{
    int Gx[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1}, Gy[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
    unsigned char data[9];

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int gx_sum = 0, gy_sum = 0;

            data[0] = matrix[i][j];
            data[1] = matrix[i][j + 1];
            data[2] = matrix[i][j + 2];
            data[3] = matrix[i + 1][j];
            data[4] = matrix[i + 1][j + 1];
            data[5] = matrix[i + 1][j + 2];
            data[6] = matrix[i + 2][j];
            data[7] = matrix[i + 2][j + 1];
            data[8] = matrix[i + 2][j + 2];

            for (int k = 0; k < 9; k++)
            {
                gx_sum += data[k] * Gx[k];
                gy_sum += data[k] * Gy[k];
            }

            // Esta parte del código es opcional ya que sirve para reducir la intensidad de los bordes
            gx_sum /= 4;
            gy_sum /= 4;

            int processed_bit = sqrt(gx_sum * gx_sum + gy_sum * gy_sum);

            result[i][j] = processed_bit;
        }
    }
}