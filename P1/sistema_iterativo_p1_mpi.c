#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define N 15000

void read_binary_data(double **matrix, const char *file_name)
{
    FILE *data;
    data = fopen(file_name, "rb");

    for (int i = 0; i < N; i++)
    {
        fread(matrix[i], sizeof(double), N, data);
    }
}

void write_binary_data(double **matrix, const char *file_name)
{
    FILE *output;

    output = fopen(file_name, "wb");

    for (int i = 0; i < N; i++)
    {
        fwrite(matrix[i], sizeof(double), N, output);
    }

    fclose(output);
}

// Convierte un valor con signo a absoluto mediante un operador ternario
double absolute(double n)
{
    return n > 0 ? n : -n;
}

// Encuentra el mayor valor absoluto en un vector y devuelve el índice en el que se encuentra
double find_local_abs_max(double vec[], int size)
{
    double max = absolute(vec[0]), current = 0;
    int position = 0;
    for (int i = 1; i < size; i++)
    {
        current = absolute(vec[i]);
        if (current > max)
        {
            position = i;
            max = current;
        }
    }
    return vec[position];
}

// Divido cada componente de un vector por el parámetro max
void divide(double vec[], double max, int size)
{
    for (int i = 0; i < size; i++)
    {
        vec[i] /= max;
    }
}

// Reserva de memoria en bloque para matriz
double **allocate_matrix(int width, int height)
{
    double **matrix, *data;
    matrix = (double **)malloc(sizeof(double *) * height);
    data = (double *)malloc(sizeof(double) * width * height);
    for (int i = 0; i < height; i++)
    {
        matrix[i] = &(data[i * width]);
    }
    return matrix;
}

// Rellena la matriz
void fill_matrix(double **matrix, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            // Diagonal matriz
            if (j == i)
            {
                matrix[i][j] = 1;
            }
            // Mitad triangular superior
            else if (i > j)
            {
                matrix[i][j] = (double)-50 * (i + 1) * (j + 1) / (n + n);
            }
            // Mitad triangular inferior
            else
            {
                matrix[i][j] = (double)50 * (i + 1) * (j + 1) / (n + n);
            }
        }
    }
}

/*
 * Operación MPI para encontrar el máximo absoluto
 * Se envían dos direcciones de memoria de tipo void porque la función allreduce sólo admite punteros a void
 * Se hace el casting de puntero a double y la desreferencia de los punteros para obtener el valor
 */
void abs_max(void *in, void *inout, int *len, MPI_Datatype *type)
{
    double n = *(double *)in, r = *(double *)inout, m;
    m = (abs(n) > abs(r)) ? n : r;
    *(double *)inout = m;
}

// Los parámetros son: m, fichero de salida de texto, fichero binario que almacena la matriz
int main(int argc, char *argv[])
{
    // Número de iteraciones
    int m = atoi(argv[1]);

    // Fichero de entrada y salida
    FILE *output_text_file = fopen(argv[2], "w");

    // Estructuras de datos globales
    double **global_matrix;

    // Parámetros MPI
    int nprocess;
    int myrank;
    int tag = 0;
    int root = 0;

    // Parámetros para medir el tiempo
    double global_start, global_end;
    double local_start;
    double execution_start;

    // Inicio MPI
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocess);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // Parámetros para los procesos
    int batch;              // Filas por procesos
    double max_value;       // Contendrá el valor máximo de cada iteración
    double local_max_value; // Valor local máximo de producto 
    double **local_matrix;  // Sub matrices locales
    double local_result[N]; // Producto matriz-vector resultante
    double u_vector[N];     // Vector unitario

    // Tipo de dato para las filas de las matrices
    MPI_Datatype row_type;
    MPI_Type_contiguous(N, MPI_DOUBLE, &row_type);
    MPI_Type_commit(&row_type);

    // Operación personalizada para allreduce
    MPI_Op MPI_abs_max;
    MPI_Op_create(abs_max, 1, &MPI_abs_max);

    MPI_Barrier(MPI_COMM_WORLD);
    global_start = MPI_Wtime();

    // Se carga la matriz en el proceso 0
    if (myrank == root)
    {
        global_matrix = allocate_matrix(N, N);

        // Si no se envía el nombre de la matriz como parámetro, se rellena y guarda
        if (argc < 4)
        {
            fill_matrix(global_matrix, N);
            write_binary_data(global_matrix, "matrix.raw");
        }
        else
        {
            read_binary_data(global_matrix, argv[3]);
        }
    }

    local_start = MPI_Wtime();

    batch = N / nprocess; // Se asume el número de filas correspondientes para cada proceso es un número entero positivo

    // Sub matrices locales
    local_matrix = allocate_matrix(N, batch);

    // Se reparte la información de la matriz original en sub matrices locales
    if (myrank == 0)
    {
        MPI_Scatter(&global_matrix[0][0], batch, row_type, &local_matrix[0][0], batch, row_type, 0, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Scatter(NULL, batch, row_type, &local_matrix[0][0], batch, row_type, 0, MPI_COMM_WORLD);
    }

    // Se inicializan los vectores unitarios de cada proceso
    for (int i = 0; i < N; i++)
    {
        u_vector[i] = 1;
    }

    execution_start = MPI_Wtime();
    // Inicio del proceso iterativo
    for (int k = 0; k < m; k++)
    {
        double sum = 0;
        for (int i = 0; i < batch; i++)
        {
            for (int j = 0; j < N; j++)
            {
                sum += local_matrix[i][j] * u_vector[j];
            }
            local_result[i] = sum;
            sum = 0;
        }

        if (k > 0)
        {
            // Se almacena el valor máximo local de cada producto matriz-vectori local
            double local_max_value = find_local_abs_max(local_result, batch);

            // Mediante allreduce, encontramos el valor máximo entre los n posibles y se distribuyen a todos los procesos, equivale a un reduce + bcast
            MPI_Allreduce(&local_max_value, &max_value, 1, MPI_DOUBLE, MPI_abs_max, MPI_COMM_WORLD);

            // El proceso 0 guarda e imprime el valor máximo junto a la iteración
            if (myrank == 0)
            {
                fprintf(output_text_file, "Iteracion: %d, valor max: %.10e\n", k, max_value);
                printf("maxvalue: %.10e\n", max_value);
            }

            // Se divide cada matriz local por el valor máximo global devuelto por el reduce
            divide(local_result, max_value, batch);
        }

        // Se distribuyen los nuevos valores computados de cada producto matriz-vector local en los vectores unitarios de cada proceso
        for (int i = 0; i < nprocess; i++)
        {
            MPI_Gather(&local_result[0], batch, MPI_DOUBLE, &u_vector[0], batch, MPI_DOUBLE, i, MPI_COMM_WORLD);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    global_end = MPI_Wtime();

    MPI_Finalize();

    // Se imprimen los datos en el fichero de salida
    if (myrank == 0)
    {
        fprintf(output_text_file, "Numero de iteraciones: %d\nNumero de procesos: %d\n", m, nprocess);
        fprintf(output_text_file, "Tiempo de global: %.2fs\n", global_end - global_start);
        fprintf(output_text_file, "Tiempo de comunicaciones iniciales y ejecucion: %.2fs\n", global_end - local_start);
        fprintf(output_text_file, "Tiempo de ejecucion: %.2fs\n", global_end - execution_start);
    }

    fclose(output_text_file);

    return 0;
}