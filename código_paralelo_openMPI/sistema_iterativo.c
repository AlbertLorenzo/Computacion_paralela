#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define N 15000

// Alumno: cp_2022_50

void reservar_memoria_matriz(double ***M_ptr, int largo, int ancho)
{
    double **matriz_aux = (double **)malloc(largo * sizeof(double *));

    if (matriz_aux == NULL)
    {
        printf("Error al reservar memoria para la matriz.\n");
        exit(1);
    }

    for (int i = 0; i < largo; i++)
    {
        matriz_aux[i] = (double *)malloc(ancho * sizeof(double));

        if (matriz_aux[i] == NULL)
        {
            printf("Error al reservar memoria para la matriz.\n");
            exit(1);
        }
    }

    *M_ptr = matriz_aux;
}

void liberar_memoria_matriz(double ***M_ptr, int largo)
{
    double **matriz_aux = *M_ptr;

    for (int i = 0; i < largo; i++)
    {
        free(matriz_aux[i]);
    }

    free(matriz_aux);
}

void escribir_matriz_datos(double **M, char *fichero)
{
    FILE *f = fopen(fichero, "wb");

    if (f == NULL)
    {
        printf("ERROR: Fallo al crear el fichero.");
        exit(1);
    }

    for (int i = 0; i < N; i++)
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

void reserva_memoria_vector(double **vector_ptr, int ancho)
{
    double *vector = (double *)malloc(sizeof(double) * ancho);

    if (vector == NULL)
    {
        printf("ERROR: Fallo al reservar memoria para el vector.");
        exit(1);
    }

    *vector_ptr = vector;
}

// Inicialización de un vector mediante punteros, en este caso, el vector unidad
void inicializar_vector_unidad(double *x0_ptr, int ancho)
{
    for (int i = 0; i < ancho; i++)
    {
        x0_ptr[i] = 1.0;
    }
}

void inicializar_matriz(double **M, int largo, int ancho, int i, int j)
{
    for (i = 0; i < largo; i++)
    {
        for (j = 0; j < ancho; j++)
        {
            if (i == j)
            {
                M[i][j] = 1.0;
            }
            else if (i > j)
            {
                M[i][j] = (double)50 * (i + 1) * (j + 1) / ((double)N * N * 10000);
            }
            else
            {
                M[i][j] = (double)-50 * (i + 1) * (j + 1) / ((double)N * N * 10000);
            }
        }
    }
}

// Encuentra el valor máximo absoluto y devuelve el índice
int indice_maximo_absoluto(double *vector, int ancho)
{
    double max = fabs(vector[0]), actual = 0;
    int indice = 0;
    for (int i = 1; i < ancho; i++)
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
void dividir(double *x0, double *x0_aux, int componentes, double maximo)
{
    for (int i = 0; i < componentes; i++)
    {
        x0[i] = x0_aux[i] / maximo;
    }
}

void producto_matriz_vector(double **M, double *x0, double *x0_aux, int filas)
{
    double suma;
    for (int i = 0; i < filas; i++)
    {
        // Se limpia el valor de la suma
        x0_aux[i] = 0;
        for (int j = 0; j < N; j++)
        {
            // Multiplico el valor de la matriz por la componente del vector
            x0_aux[i] += (M[i][j] * x0[j]);
        }
    }
}

int main(int argc, char *argv[])
{
    int myrank;
    int numero_procesos;
    const int ROOT = 0;
    const int ETIQUETA_COMUNICACIONES = 50;
    clock_t tiempo_inicio_sistema, tiempo_fin_sistema, tiempo_lectura_escritura_matriz_inicio, tiempo_lectura_escritura_matriz_fin;
    double tiempo_MPI_inicio, tiempo_MPI_fin;

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numero_procesos);

    FILE *fichero_salida;

    double **matriz_datos;
    double **M;

    int m = atoi(argv[3]);
    int resto_filas = N % numero_procesos;
    int filas_procesos = (myrank == ROOT) ? N / numero_procesos + resto_filas : N / numero_procesos;

    MPI_Barrier(MPI_COMM_WORLD);
    tiempo_MPI_inicio = MPI_Wtime();

    reservar_memoria_matriz(&M, filas_procesos, N);

    if (myrank == ROOT)
    {
        if (argc != 4)
        {
            printf("La lista de argumentos no es correcta: <fichero_entrada> <fichero_resultados> <numero_iteraciones>");
            exit(1);
        }

        reservar_memoria_matriz(&matriz_datos, N, N);

        tiempo_lectura_escritura_matriz_inicio = clock();

        if (comprobar_fichero_existe(argv[1]))
        {
            leer_matriz_datos(matriz_datos, argv[1]);
            printf("Se ha completado la lectura del fichero de entrada.\n");
        }
        else
        {
            inicializar_matriz(matriz_datos, N, N, 0, 0);
            escribir_matriz_datos(matriz_datos, argv[1]);
            printf("Se ha generado el fichero.\n");
        }

        tiempo_lectura_escritura_matriz_fin = clock();

        int iteraciones = filas_procesos - resto_filas;
        for (int i = 0; i < numero_procesos; i++)
        {
            if (i > 0)
            {
                for (int j = 0; j < iteraciones; j++)
                {
                    MPI_Send(matriz_datos[(iteraciones * i) + j + resto_filas], N, MPI_DOUBLE, i, ETIQUETA_COMUNICACIONES, MPI_COMM_WORLD);
                }
            }
            else
            {
                for (int j = 0; j < filas_procesos; j++)
                {
                    memcpy(M[j], matriz_datos[j], sizeof(double) * N);
                }
            }
        }

        liberar_memoria_matriz(&matriz_datos, N);
        tiempo_inicio_sistema = clock();
    }
    else
    {
        for (int i = 0; i < filas_procesos; i++)
        {
            MPI_Recv(M[i], N, MPI_DOUBLE, ROOT, ETIQUETA_COMUNICACIONES, MPI_COMM_WORLD, &status);
        }
    }

    // Se declara, reserva e inicializan todas las variables necesarias para el sistema
    int posiciones_procesos[numero_procesos], posiciones_resultantes[m - 1], ajustes[m - 1];
    double *x0, *x0_aux, *maximos_procesos, *maximos_resultantes, maximo_local;
    reserva_memoria_vector(&x0, N);
    reserva_memoria_vector(&x0_aux, N);
    reserva_memoria_vector(&maximos_resultantes, m - 1);
    reserva_memoria_vector(&maximos_procesos, numero_procesos);

    inicializar_vector_unidad(x0, N);

    // Siguientes iteraciones
    for (int i = 0; i < m; i++)
    {
        // Se realiza el producto matriz vector
        producto_matriz_vector(M, x0, x0_aux, filas_procesos);

        // Se copian los resultados de x0_aux a x0
        memcpy(&x0[0], &x0_aux[0], sizeof(double) * filas_procesos);

        if (i > 0)
        {
            // Se calcula el índice del valor máximo absoluto para cada proceso
            int indice = indice_maximo_absoluto(x0, filas_procesos);

            // Se calcula su índice real respecto al tamaño N = 15000
            int indice_calc = myrank * filas_procesos + indice + resto_filas;

            // Se actualiza el valor máximo absoluto para cada proceso
            maximo_local = x0[indice];

            // Se recoge el máximo y el índice local de cada proceso
            MPI_Gather(&maximo_local, 1, MPI_DOUBLE, maximos_procesos, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
            MPI_Gather(&indice_calc, 1, MPI_INT, posiciones_procesos, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

            if (myrank == ROOT)
            {
                // Se calcula el índice del máximo absoluto de todos los máximos
                int indice_global_maximo = indice_maximo_absoluto(maximos_procesos, numero_procesos);
                maximo_local = maximos_procesos[indice_global_maximo];

                // Se guarda el valor y el índice del máximo absoluto global
                posiciones_resultantes[i - 1] = posiciones_procesos[indice_global_maximo];
                maximos_resultantes[i - 1] = maximo_local;
            }

            // Se distribuye el máximo global calculado por ROOT al resto de procesos
            MPI_Bcast(&maximo_local, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

            // Si el máximo > 25.0 se divide cada componente por dicho valor
            if (fabs(maximo_local) > 25.0)
            {
                dividir(x0, x0_aux, filas_procesos, maximo_local);
                ajustes[i] = 1;
            }
        }

        // Se envían, a expcepción de ROOT, los vectores resultantes del producto matriz-vector de cada proceso a ROOT para obtener el vector completo
        if (myrank >= 1)
        {
            MPI_Send(&x0[0], filas_procesos, MPI_DOUBLE, ROOT, ETIQUETA_COMUNICACIONES, MPI_COMM_WORLD);
        }
        else
        {
            for (int i = 1; i < numero_procesos; i++)
            {
                MPI_Recv(&x0[(filas_procesos * i) - (resto_filas * (i - 1))], (filas_procesos - resto_filas), MPI_DOUBLE, i, ETIQUETA_COMUNICACIONES, MPI_COMM_WORLD, &status);
            }
        }

        // Se distribuye el vector x0 completamente calculado al resto de procesos
        MPI_Bcast(x0, N, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    tiempo_MPI_fin = MPI_Wtime();

    // Fin iteraciones
    if (myrank == 0)
    {
        tiempo_fin_sistema = clock();

        fichero_salida = fopen(argv[2], "w");
        fprintf(fichero_salida, "Dimension N: %d\n", N);
        printf("Dimension N: %d\n", N);
        for (int i = 0; i < m - 1; i++)
        {
            printf("Máximo en iteracion %d: %.1f Posicion: %d Ajuste: %s\n", i + 1, maximos_resultantes[i], posiciones_resultantes[i], (ajustes[i + 1] == 1) ? "Si" : "No");
            fprintf(fichero_salida, "Máximo en iteracion %d: %.1f Posicion: %d Ajuste: %s\n", i + 1, maximos_resultantes[i], posiciones_resultantes[i], (ajustes[i + 1] == 1) ? "Si" : "No");
        }
        printf("Tiempo lectura o escritura de la matriz: %.2f seg\n", ((double)(tiempo_lectura_escritura_matriz_fin - tiempo_lectura_escritura_matriz_inicio)) / CLOCKS_PER_SEC);
        printf("Tiempo de ejecucion del sistema: %.2f seg\n", ((double)(tiempo_fin_sistema - tiempo_inicio_sistema)) / CLOCKS_PER_SEC);
        printf("Tiempo MPI: %.2f seg\n", tiempo_MPI_fin - tiempo_MPI_inicio);
        fprintf(fichero_salida, "Tiempo lectura o escritura de la matriz: %.2f seg\n", ((double)(tiempo_lectura_escritura_matriz_fin - tiempo_lectura_escritura_matriz_inicio)) / CLOCKS_PER_SEC);
        fprintf(fichero_salida, "Tiempo de ejecucion del sistema: %.2f seg\n", ((double)(tiempo_fin_sistema - tiempo_inicio_sistema)) / CLOCKS_PER_SEC);
        fprintf(fichero_salida, "Tiempo MPI: %.2f seg\n", tiempo_MPI_fin - tiempo_MPI_inicio);
        fprintf(fichero_salida, "Numero de procesos: %d\nM: %d\n", numero_procesos, m);
        
        fclose(fichero_salida);
    }

    // Se libera la memoria
    liberar_memoria_matriz(&M, filas_procesos);
    free(x0);
    free(x0_aux);
    free(maximos_procesos);
    free(maximos_resultantes);

    MPI_Finalize();

    return 0;
}