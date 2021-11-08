#include <iostream>
#include <fstream>
#include <random>
#include <time.h>
#include <chrono>

#define N 15000

void print_matrix(double **M)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            std::cout << M[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void print_vector(double vector[])
{
    for (int i = 0; i < N; i++)
    {
        std::cout << vector[i] << " ";
    }
    std::cout << std::endl;
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
double abs(double n) {
    return n > 0 ? n : -n;
}

// Encuentra el mayor valor absoluto en un vector y devuelve el índice en el que se encuentra
int find_max_absolute(double vec[])
{
    double max = abs(vec[0]), current = 0;
    int position = 0;
    for (int i = 1; i < N; i++)
    {
        current = abs(vec[i]);
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
    double **matrix = new double *[N];

    for (auto i = 0; i < N; ++i)
    {
        matrix[i] = new double[N];
    }

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

    *M_ptr = matrix;
}

// Inicialización de un vector mediante punteros, en este caso, el vector unidad
void init_unity_vector(double *x0_ptr[])
{
    double *vector_aux = new double[N];
    for (int i = 0; i < N; i++)
    {
        vector_aux[i] = 1;
    }
    *x0_ptr = vector_aux;
}

// Función para iniciar el sistema iterativo
void init_system(double **M, double x0[], int m, std::ofstream &output)
{

    double vector_aux[N], sum;

    for (int k = 0; k < m; k++)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                // Multiplico el valor de la matriz por la componente del vector
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

            output << "Iteracion n: " << k << " valor maximo: " << vector_aux[position] << " indice: " << position << "\n";

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
    int m = std::stoi(argv[1]);
    std::ofstream output("informacion_sistema_iterativo.txt");

    auto global_begin = std::chrono::high_resolution_clock::now();

    auto matrix_begin = std::chrono::high_resolution_clock::now();
    init_matrix(&M);
    auto matrix_end = std::chrono::high_resolution_clock::now();

    init_unity_vector(&x0);
    
    auto system_begin = std::chrono::high_resolution_clock::now();
    init_system(M, x0, m, output);
    auto system_end = std::chrono::high_resolution_clock::now();
    
    auto global_end = std::chrono::high_resolution_clock::now();

    auto global_duration = std::chrono::duration_cast<std::chrono::microseconds>(global_end - global_begin);
    auto system_duration = std::chrono::duration_cast<std::chrono::microseconds>(system_end - system_begin);
    auto matrix_duration = std::chrono::duration_cast<std::chrono::microseconds>(matrix_end - matrix_begin);


    output << "Parametros de ejecucion: argv[0] = ejecutable, argv[1] = numero de iteraciones \n"; 
    output << "Numero de iteraciones: " << m << "\n";
    output << "Tiempo global <microsegundos>: " << global_duration.count() << "\n";
    output << "Tiempo ejecucion <microsegundos>: " << system_duration.count() << "\n";
    output << "Tiempo de generacion de la matriz <microsegundos>: " << matrix_duration.count() << "\n";
    output << "Fichero de salida: " << "informacion_sistema_iterativo.txt" << "\n";

    return 0;
}