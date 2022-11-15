#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define N 15000

void escribir_matriz_datos(double **M, char *fichero)
{
    int i, j;
    FILE *f = fopen(fichero, "wb");

    if (f == NULL)
    {
        printf("ERROR: Fallo al crear el fichero.");
        exit(1);
    }

    for (i = 0; i < N; i++)
    {
        fwrite(M[i], sizeof(*M[i]), N, f);
    }

    fclose(f);
}

void leer_matriz_datos(double **M, char *fichero)
{
    FILE *f = fopen(fichero, "rb");

    if (f == NULL)
    {
        printf("ERROR: Fallo al abrir el fichero.");
        exit(1);
    }

    for (int i = 0; i < N; i++)
    {
        fread(M[i], sizeof(*M[i]), N, f);
    }

    fclose(f);
}

int comprobar_fichero_existe(char *fichero)
{
    FILE *f = fopen(fichero, "rb");

    if (f == NULL)
    {
        return 0;
    }

    fclose(f);
    return 1;
}

// Inicialización de una matriz mediante punteros
void reserva_memoria_matriz(double ***matriz_ptr)
{
    // Se declara y reserva memoria para la matriz
    double **matriz = (double **)malloc(sizeof(double *) * N);

    if (matriz == NULL)
    {
        printf("ERROR: Fallo al reservar memoria para la matriz.");
        exit(1);
    }

    for (int i = 0; i < N; ++i)
    {
        matriz[i] = (double *)malloc(sizeof(double) * N);

        if (matriz[i] == NULL)
        {
            printf("ERROR: Fallo al reservar memoria para la matriz.");
            exit(1);
        }
    }

    *matriz_ptr = matriz;
}

void reserva_memoria_vector(double **vector_ptr)
{
    double *vector = (double *)malloc(sizeof(double) * N);

    if (vector == NULL)
    {
        printf("ERROR: Fallo al reservar memoria para el vector.");
        exit(1);
    }

    *vector_ptr = vector;
}

// Inicialización de un vector mediante punteros, en este caso, el vector unidad
void inicializar_vector_unidad(double *x0_ptr)
{
    for (int i = 0; i < N; i++)
    {
        x0_ptr[i] = 1.0;
    }
}

void inicializar_matriz(double **M)
{
    int i, j;

    double aux;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            if (i == j)
            {
                M[i][j] = 1.0;
            }
            else if (i > j)
            {
                M[i][j] = (double)50 * (i + 1) * (j + 1) / ((double)N * N * 10000); // Triangular inferior
            }
            else
            {
                M[i][j] = (double)-50 * (i + 1) * (j + 1) / ((double)N * N * 10000); // Triangular superior
            }
        }
    }
}

// Encuentra el valor máximo absoluto y devuelve el índice
int indice_maximo_absoluto(double *vector)
{
    double max = fabs(vector[0]), actual = 0;
    int indice = 0;
    for (int i = 1; i < N; i++)
    {
        actual = fabs(vector[i]);
        if (actual > max)
        {
            indice = i;
            max = actual;
        }
    }
    return indice;
}

// Se divide cada componente de un vector por un escalar
void dividir(double *vector, double maximo)
{
    for (int i = 0; i < N; i++)
    {
        vector[i] /= maximo;
    }
}

void producto_matriz_vector(double **M, double *xi, double *vector_aux)
{
    double suma;
    for (int i = 0; i < N; i++)
    {
        // Se limpia el valor de la suma
        suma = 0;
        for (int j = 0; j < N; j++)
        {
            // Multiplico el valor de la matriz por la componente del vector
            suma += (M[i][j] * xi[j]);
        }
        // Actualizo la componente del vector
        vector_aux[i] = suma;
    }
}

/*
* Si el valor absoluto es mayor a 25, se dividen las componentes del vector por ese valor máximo.
* En caso de no serlo, se guarda la información de la iteración sin aplicar ningúna transformación
*/
void convertir(double *vector, FILE *fichero_salida)
{
    int posicion = indice_maximo_absoluto(vector);
    double maximo = vector[posicion];
    if (abs(maximo) > 25.0)
    {
        fprintf(fichero_salida, "Maximo absoluto: %.1f en la posicion %d. Hago cambio\n", vector[posicion], posicion);
        dividir(vector, maximo);
    }
    else
    {
        fprintf(fichero_salida, "Maximo absoluto: %.1f en la posicion %d.\n", vector[posicion], posicion);
    }
}

void inicializar_sistema(double **M, double *xi, int m, FILE *fichero_salida)
{
    double *vector_aux;
    reserva_memoria_vector(&vector_aux);
    for (int i = 0; i < m; i++)
    {
        producto_matriz_vector(M, xi, vector_aux);
        if (i > 0) {
            convertir(vector_aux, fichero_salida);
        }
        // Se copia el vector auxiliar en el vector xi para guardar el anterior xi
        memcpy(xi, vector_aux, sizeof(double) * N);    
    }
    free(vector_aux);
}

void liberar_memoria_matriz(double **matriz)
{
    for (int i = 0; i < N; ++i)
    {
        free(matriz[i]);
    }
    free(matriz);
}

int main(int argc, char *argv[])
{
    double **M;
    double *x0;

    if (argc != 4)
    {
        printf("ERROR: Numero de argumentos incorrecto.\n1: Nombre del fichero de entrada.\n2: Nombre del fichero de salida.\n3: Numero de iteraciones.\n");
        exit(1);
    }

    reserva_memoria_matriz(&M);

    clock_t inicio_global = clock();
    if (comprobar_fichero_existe(argv[1]))
    {
        leer_matriz_datos(M, argv[1]);
    }
    else
    {
        inicializar_matriz(M);
        escribir_matriz_datos(M, argv[1]);
    }

    FILE *fichero_salida = fopen(argv[2], "w");


    reserva_memoria_vector(&x0);
    inicializar_vector_unidad(x0);

    clock_t inicio_sistema = clock();
    inicializar_sistema(M, x0, atoi(argv[3]), fichero_salida);
    clock_t fin_sistema = clock();

    fprintf(fichero_salida, "Tiempo de ejecucion sistema (segundos): %f\n", (double)(fin_sistema - inicio_sistema) / CLOCKS_PER_SEC);
    fprintf(fichero_salida, "Tiempo de ejecucion global (segundos): %f\n", (double)(fin_sistema - inicio_global) / CLOCKS_PER_SEC);

    liberar_memoria_matriz(M);
    free(x0);
    fclose(fichero_salida);

    return 0;
}