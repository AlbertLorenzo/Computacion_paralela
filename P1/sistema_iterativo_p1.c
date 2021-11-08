#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 15000

void write_matrix_data(double **matrix)
{
    FILE *output;

    output = fopen("matriz.txt", "w");

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            fprintf(output, "%f ", matrix[i][j]);
        }
        fprintf(output, "\n");
    }

    fclose(output);
}

void read_matrix_from_file(double **matrix, const char *file_name)
{
    FILE *data;
    data = fopen(file_name, "r");

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            fscanf(data, "%lf ", &matrix[i][j]);
        }
    }
}

// Función para copiar valores de un array a otro
void copy_vector_values(double vec_a[], double vec_b[])
{
    for (int i = 0; i < N; i++)
    {
        vec_b[i] = vec_a[i];
    }
}

// Convierte un valor con signo a absoluto mediante un operador ternario
double absolute(double n)
{
    return n > 0 ? n : -n;
}

// Encuentra el mayor valor absoluto en un vector y devuelve el índice en el que se encuentra
int find_max_absolute(double vec[])
{
    double max = absolute(vec[0]), current = 0;
    int position = 0;
    for (int i = 1; i < N; i++)
    {
        current = absolute(vec[i]);
        if (current > max)
        {
            position = i;
            max = current;
        }
    }
    return position;
}

// Divido cada componente de un vector por el parámetro max
void divide(double vec[], double max)
{
    for (int i = 0; i < N; i++)
    {
        vec[i] /= max;
    }
}

// Inicialización de una matriz mediante punteros
void init_matrix(double ***M_ptr)
{
    // Se declara y reserva memoria para la matriz
    double **matrix = (double **)malloc(sizeof(double *) * N);

    for (int i = 0; i < N; ++i)
    {
        matrix[i] = (double *)malloc(sizeof(double) * N);
    }

    *M_ptr = matrix;
}

void fill_matrix(double **matrix)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            // Diagonal matriz
            if (j == i)
            {
                matrix[i][j] = 1;
            }
            // Mitad triangular superior
            else if (i > j)
            {
                matrix[i][j] = (double)-50 * (i + 1) * (j + 1) / (N + N);
            }
            // Mitad triangular inferior
            else
            {
                matrix[i][j] = (double)50 * (i + 1) * (j + 1) / (N + N);
            }
        }
    }
}

// Inicialización de un vector mediante punteros, en este caso, el vector unidad
void init_unity_vector(double *x0_ptr[])
{
    double *vector_aux = (double *)malloc(sizeof(double) * N);
    for (int i = 0; i < N; i++)
    {
        vector_aux[i] = 1;
    }
    *x0_ptr = vector_aux;
}

// Función para iniciar el sistema iterativo
void init_system(double **M, double x0[], int m, FILE *result)
{

    double vector_aux[N], sum;

    for (int k = 0; k < m; k++)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                sum += M[i][j] * x0[j];
            }
            // Actualizo la componente del vector y limpio el valor de la suma, así en una nueva iteración no arrastro el valor de la componente antigua
            vector_aux[i] = sum;
            sum = 0;
        }

        if (k > 0)
        {
            // Encuentra y devuelve el índice de la componente con un mayor valor absoluto dentro de un vector
            int position = find_max_absolute(vector_aux);

            fprintf(result, "Iteracion: %d, valor max: %f, posicion: %d\n", k, vector_aux[position], position);

            // dividir los valores del vector por el máximo absoluto, respetando el signo original
            divide(vector_aux, vector_aux[position]);
        }

        // Copio los valores del vector resultante en el vector inicial para que afecte a las siguientes iteraciones
        copy_vector_values(vector_aux, x0);
    }
}

int main(int argc, char *argv[])
{
    double **M;
    double *x0;
    int m = atoi(argv[1]);
    FILE *result = fopen("informacion_sistema_iterativo.txt", "w");

    init_unity_vector(&x0);
    init_matrix(&M);

    if (argc < 3)
    {
        fill_matrix(M);
        write_matrix_data(M);
    }
    else
    {
        read_matrix_from_file(M, argv[2]);
    }

    init_system(M, x0, m, result);
    return 0;
}