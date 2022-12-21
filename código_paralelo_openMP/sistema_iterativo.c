#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <omp.h>

#define N 15000

// Gestión de ficheros
int comprobar_fichero_existe(char *);
void leer_fichero(char *, double **, int);
void escribir_fichero(double **, int, char *);

// Gestión de memoria
void reservar_memoria_matriz(double ***, int);
void reserva_memoria_vector(double **, int);

// Inicialización de datos
void inicializar_vector_unidad(double *, int);
void inicializar_matriz(double **);

int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        printf("ERROR: Numero de argumentos incorrecto.\n1: Nombre del fichero de entrada.\n2: Nombre del fichero de salida.\n3: Numero de iteraciones.\n4. Numero de hilos.\n");
        exit(1);
    }

    // Parámetros de entrada
    char *fichero_entrada = argv[1];
    char *fichero_salida = argv[2];
    int iteraciones = atoi(argv[3]);
    int num_hilos = atoi(argv[4]);

    // Parámetros de distribución de carga
    int id_hilo;
    // int componentes_hilo = num_hilos / N;
    // int resto_hilos = num_hilos % N;
    double maximo_absoluto = 0;
    double maximo = 0;
    double minimo = 0;

    // Parámetros de tiempo
    double tiempo_inicio, tiempo_final;

    // Estructuras de datos
    double **M;
    double *x0, *x1;
    double *vector_maximos;

    // Reserva de memoria para la matriz y los vectores
    reservar_memoria_matriz(&M, N);
    reserva_memoria_vector(&x0, N);
    reserva_memoria_vector(&x1, N);
    reserva_memoria_vector(&vector_maximos, iteraciones - 1);

    // Si el fichero existe se lee, si no, se inicializa la matriz
    if (comprobar_fichero_existe(fichero_entrada))
    {
        leer_fichero(fichero_entrada, M, N);
        printf("Lectura de la matriz de datos realizada correctamente.\n");
    }
    else
    {
        inicializar_matriz(M);
        escribir_fichero(M, N, fichero_entrada);
        printf("Se ha creado el fichero de datos correctamente.\n");
    }

    // Inicialización del vector x0
    inicializar_vector_unidad(x0, N);

    tiempo_inicio = omp_get_wtime();

    #pragma omp parallel num_threads(num_hilos) private(id_hilo) shared(M, x0, x1, iteraciones, minimo, maximo, maximo_absoluto, vector_maximos)
    {
        int componentes_hilo = N / num_hilos;
        
        for (int m = 0; m < iteraciones; m++)
        {
            // Se reparten las filas de la matriz entre los hilos
            #pragma omp for schedule(static, (componentes_hilo))
            for (int i = 0; i < N; i++)
            {
                double suma = 0;
                for (int j = 0; j < N; j++)
                {
                    suma += M[i][j] * x0[j];
                }
                x1[i] = suma;
            }

            if (m > 0) {
                // Se hace un reduce con el operador max y min para obtener el máximo y el mínimo de los valores del vector
                #pragma omp for schedule(static, (componentes_hilo)) reduction(max:maximo) reduction(min:minimo)
                for(int i = 0; i < N; i++) 
                {
                    if (x1[i] > maximo)
                    {
                        maximo = x1[i];    
                    }

                    if (x1[i] < minimo)
                    {
                        minimo = x1[i];    
                    }
                }

                // El primer hilo disponible se encarga de calcular el máximo absoluto y guardarlo en el vector de máximos
                #pragma omp single
                {
                    maximo_absoluto = maximo > fabs(minimo) ? maximo : minimo;
                    vector_maximos[m - 1] = maximo_absoluto;
                    maximo = 0;
                    minimo = 0;
                }

                // Se normalizan los valores del vector si el valor es mayor que 25
                if (fabs(maximo_absoluto) > 25) {
                    #pragma omp for schedule(static, (componentes_hilo))
                    for(int i = 0; i < N; i++) 
                    {
                        x1[i] /= maximo_absoluto;
                    }
                }
            }
            
            #pragma omp for schedule(static, (componentes_hilo))
            for(int i = 0; i < N; i++) 
            {
                x0[i] = x1[i];
            }
        }
    }

    tiempo_final = omp_get_wtime();

    FILE *salida = fopen(fichero_salida, "w");

    for(int i = 0; i < iteraciones - 1; i++)
    {
        fprintf(salida, "Iteracion: %d, Maximo: %.1f\n", i + 1, vector_maximos[i]);
    }

    fprintf(salida, "Nombre fichero de datos: %s\n", fichero_entrada);
    fprintf(salida, "Nombre fichero de salida: %s\n", fichero_salida);
    fprintf(salida, "Numero de iteraciones: %d\n", iteraciones);
    fprintf(salida, "Numero de hilos: %d\n", num_hilos);
    fprintf(salida, "Tiempo: %fs\n", tiempo_final - tiempo_inicio);

    printf("Fichero de salida generado correctamente.\n");

    fclose(salida);

    return 0;
}

void reservar_memoria_matriz(double ***matriz, int dimension)
{
    *matriz = (double **)malloc(dimension * sizeof(double *));

    if (*matriz == NULL)
    {
        printf("Error al reservar memoria para la matriz\n");
        exit(1);
    }

    for (int i = 0; i < dimension; i++)
    {
        (*matriz)[i] = (double *)malloc(dimension * sizeof(double));

        if ((*matriz)[i] == NULL)
        {
            printf("Error al reservar memoria para la matriz\n");
            exit(1);
        }
    }
}

void leer_fichero(char *nombre_fichero, double **matriz, int dimension)
{
    FILE *fichero;

    fichero = fopen(nombre_fichero, "rb");

    if (fichero == NULL)
    {
        printf("Error al abrir el fichero\n");
        exit(1);
    }

    for (int i = 0; i < dimension; i++)
    {
        fread(matriz[i], sizeof(double), dimension, fichero);
    }

    fclose(fichero);
}

void escribir_fichero(double **M, int dimension, char *fichero)
{
    FILE *f = fopen(fichero, "wb");

    if (f == NULL)
    {
        printf("ERROR: Fallo al crear el fichero.");
        exit(1);
    }

    for (int i = 0; i < dimension; i++)
    {
        fwrite(M[i], sizeof(double), dimension, f);
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

void reserva_memoria_vector(double **vector_ptr, int dimension)
{
    double *vector = (double *)malloc(sizeof(double) * dimension);

    if (vector == NULL)
    {
        printf("ERROR: Fallo al reservar memoria para el vector.");
        exit(1);
    }

    *vector_ptr = vector;
}

void inicializar_vector_unidad(double *x0_ptr, int dimension)
{
    for (int i = 0; i < dimension; i++)
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