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

void write_binary_data(double **matrix)
{
    FILE *output;

    output = fopen("matrix.raw", "wb");

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
double find_max_absolute(double vec[], int size)
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

// Reserva e inicialización el vector unidad
void init_unity_vector(double *x0_ptr[], int size)
{
    double *vector_aux = (double *)malloc(sizeof(double) * size);
    for (int i = 0; i < size; i++)
    {
        vector_aux[i] = 1;
    }
    *x0_ptr = vector_aux;
}

// Operación MPI para encontrar el máximo absoluto
void abs_max(void *in, void *inout, int *len, MPI_Datatype *type)
{
    double n = *(double *)in, r = *(double *)inout, m;
    m = (abs(n) > abs(r)) ? n : r;
    *(double *)inout = m;
}

int main(int argc, char *argv[])
{
    // Número de iteraciones
    int m = atoi(argv[1]);

    // Fichero de entrada
    FILE *result = fopen("informacion_sistema_iterativo.txt", "w");

    // Estructuras de datos globales
    double **global_matrix;

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

    // Parámetros para los procesos
    int batch;              // Filas por procesos
    double max_value;       // Contendrá el valor máximo de cada iteración
    double local_max_value; // Valor local máximo de cada proceso
    double **local_matrix;  // Sub matrices locales
    double local_result[N]; // Producto vectorial resultante
    double u_vector[N];     // Vector unitario

    // Tipo de dato para las filas de las matrices
    MPI_Datatype row_type;
    MPI_Type_contiguous(N, MPI_DOUBLE, &row_type);
    MPI_Type_commit(&row_type);

    // Operación personalizada para allreduce
    MPI_Op calc_abs_max;
    MPI_Op_create(abs_max, 1, &calc_abs_max);

    if (myrank == root)
    {
        global_matrix = allocate_matrix(N, N);

        if (argc < 3)
        {
            fill_matrix(global_matrix, N);
            write_binary_data(global_matrix);
        }
        else
        {
            read_binary_data(global_matrix, argv[2]);
        }
    }

    batch = N / nprocess; // Se asume que el cociente de la división es un número entero positivo

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

    // Se inicializa el vector unitario
    for (int i = 0; i < N; i++)
    {
        u_vector[i] = 1;
    }

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
            // Se almacena el máximo local junto al signo original del valor en un struct
            double local_max_value = find_max_absolute(local_result, batch);

            MPI_Allreduce(&local_max_value, &max_value, 1, MPI_DOUBLE, calc_abs_max, MPI_COMM_WORLD);
            if (myrank == 0)
            {
                fprintf(result, "Iteracion: %d, valor max: %.10e\n", k, max_value);
                printf("maxvalue: %.10e\n", max_value);
            }
            divide(local_result, max_value, batch);
        }

        for (int i = 0; i < nprocess; i++)
        {
            MPI_Gather(&local_result[0], batch, MPI_DOUBLE, &u_vector[0], batch, MPI_DOUBLE, i, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();
    return 0;
}