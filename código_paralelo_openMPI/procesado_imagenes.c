#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

void reservar_memoria_matriz(unsigned char ***matriz_ptr, int largo, int ancho)
{
    unsigned char **matriz_aux = (unsigned char **)malloc(sizeof(unsigned char *) * largo);

    if (matriz_aux == NULL)
    {
        printf("Error: No se pudo reservar memoria para la matriz.\n");
        exit(1);
    }

    for (int i = 0; i < largo; i++)
    {
        matriz_aux[i] = (unsigned char *)malloc(sizeof(unsigned char) * ancho);

        if (matriz_aux[i] == NULL)
        {
            printf("Error: No se pudo reservar memoria para la matriz.\n");
            exit(1);
        }
    }

    *matriz_ptr = matriz_aux;
}

void leer_fichero_matriz(unsigned char **matriz, const char *nombre_fichero, int largo, int ancho)
{
    FILE *fichero = fopen(nombre_fichero, "rb");

    if (fichero == NULL)
    {
        printf("Error: No se pudo abrir el fichero.\n");
        exit(1);
    }

    for (int i = 0; i < largo; i++)
    {
        fread(matriz[i], sizeof(unsigned char), ancho, fichero);
    }

    fclose(fichero);
}

void volcar_fichero(unsigned char **matriz, const char *nombre_fichero_salida, int largo, int ancho)
{
    FILE *fichero_salida = fopen(nombre_fichero_salida, "wb");

    if (fichero_salida == NULL)
    {
        printf("Error: No se pudo abrir el fichero.\n");
        exit(1);
    }

    for (int i = 0; i < largo; i++)
    {
        fwrite(matriz[i], sizeof(unsigned char), ancho, fichero_salida);
    }
    fclose(fichero_salida);
}

void copiar_filas(unsigned char **matriz, unsigned char **matriz_aux, int fila_inicial, int fila_final, int ancho)
{
    for (int i = fila_inicial; i < fila_final; i++)
    {
        memcpy(matriz_aux[i], matriz[i], sizeof(unsigned char) * ancho);
    }
}

void filtro_media(unsigned char **matriz, int largo, int ancho)
{
    unsigned char **matriz_aux;
    reservar_memoria_matriz(&matriz_aux, largo, ancho);
    copiar_filas(matriz, matriz_aux, 0, largo, ancho);

    for (int i = 1; i < largo - 1; i++)
    {
        for (int j = 1; j < ancho - 1; j++)
        {
            matriz_aux[i][j] = (matriz[i - 1][j - 1] + matriz[i - 1][j] + matriz[i - 1][j + 1] + matriz[i][j - 1] + matriz[i][j] + matriz[i][j + 1] + matriz[i + 1][j - 1] + matriz[i + 1][j] + matriz[i + 1][j + 1]) / 9;
        }
    }

    for (int i = 1; i < largo - 1; i++)
    {
        for (int j = 1; j < ancho - 1; j++)
        {
            matriz[i][j] = matriz_aux[i][j];
        }
    }

    for (int i = 0; i < largo; i++)
    {
        free(matriz_aux[i]);
    }
    free(matriz_aux);
}

void filtro_mediana(unsigned char **matriz, int largo, int ancho)
{
    unsigned char **matriz_aux;
    unsigned char componentes[9];
    reservar_memoria_matriz(&matriz_aux, largo, ancho);
    copiar_filas(matriz, matriz_aux, 0, largo, ancho);

    for (int i = 1; i < largo - 1; i++)
    {
        for (int j = 1; j < ancho - 1; j++)
        {
            componentes[0] = matriz[i - 1][j - 1];
            componentes[1] = matriz[i - 1][j];
            componentes[2] = matriz[i - 1][j + 1];
            componentes[3] = matriz[i][j - 1];
            componentes[4] = matriz[i][j];
            componentes[5] = matriz[i][j + 1];
            componentes[6] = matriz[i + 1][j - 1];
            componentes[7] = matriz[i + 1][j];
            componentes[8] = matriz[i + 1][j + 1];

            for (int k = 0; k < 9; k++)
            {
                for (int l = k + 1; l < 9; l++)
                {
                    if (componentes[k] > componentes[l])
                    {
                        unsigned char aux = componentes[k];
                        componentes[k] = componentes[l];
                        componentes[l] = aux;
                    }
                }
            }

            matriz_aux[i][j] = componentes[4];
        }
    }

    for (int i = 1; i < largo - 1; i++)
    {
        for (int j = 1; j < ancho - 1; j++)
        {
            matriz[i][j] = matriz_aux[i][j];
        }
    }

    for (int i = 0; i < largo; i++)
    {
        free(matriz_aux[i]);
    }
    free(matriz_aux);
}

void extension_simetrica(unsigned char **matriz_a, unsigned char **matriz_b, int largo, int ancho)
{
    // Filas
    for (int i = 1; i < largo - 1; i++)
    {
        matriz_b[i][0] = matriz_a[i - 1][1];
        matriz_b[i][ancho - 1] = matriz_a[i - 1][ancho - 4];
    }

    // Columnas
    for (int j = 1; j < ancho - 1; j++)
    {
        matriz_b[0][j] = matriz_a[1][j - 1];
        matriz_b[largo - 1][j] = matriz_a[largo - 4][j - 1];
    }

    // Copia de datos
    for (int i = 1; i < largo - 1; i++)
    {
        for (int j = 1; j < ancho - 1; j++)
        {
            matriz_b[i][j] = matriz_a[i - 1][j - 1];
        }
    }

    // Esquinas
    matriz_b[0][0] = matriz_b[2][2];
    matriz_b[0][ancho - 1] = matriz_b[2][ancho - 3];
    matriz_b[largo - 1][0] = matriz_b[largo - 3][2];
    matriz_b[largo - 1][ancho - 1] = matriz_b[largo - 3][ancho - 3];
}

void filtro_sobel(unsigned char **matriz, int largo, int ancho)
{
    unsigned char **matriz_aux;
    reservar_memoria_matriz(&matriz_aux, largo, ancho);
    extension_simetrica(matriz, matriz_aux, largo, ancho);

    int c[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
    int f[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
    unsigned char componentes[9];

    for (int i = 1; i < largo - 1; i++)
    {
        for (int j = 1; j < ancho - 1; j++)
        {
            int suma_c = 0;
            int suma_f = 0;

            componentes[0] = matriz_aux[i - 1][j - 1];
            componentes[1] = matriz_aux[i - 1][j];
            componentes[2] = matriz_aux[i - 1][j + 1];
            componentes[3] = matriz_aux[i][j - 1];
            componentes[4] = matriz_aux[i][j];
            componentes[5] = matriz_aux[i][j + 1];
            componentes[6] = matriz_aux[i + 1][j - 1];
            componentes[7] = matriz_aux[i + 1][j];
            componentes[8] = matriz_aux[i + 1][j + 1];

            // Se calcula el valor de cada suma
            for (int k = 0; k < 9; k++)
            {
                suma_c += c[k] * componentes[k];
                suma_f += f[k] * componentes[k];
            }

            // Se calcula el valor del pixel actual
            matriz[i - 1][j - 1] = sqrt(suma_c * suma_c + suma_f * suma_f);
        }
    }

    for (int i = 0; i < largo; i++)
    {
        free(matriz_aux[i]);
    }
    free(matriz_aux);
}

void filtro(const char *tipo_filtro, unsigned char **matriz, int largo, int ancho)
{
    if (strcmp(tipo_filtro, "sobel") == 0)
    {
        filtro_sobel(matriz, largo + 2, ancho + 2);
    }

    if (strcmp(tipo_filtro, "media") == 0)
    {
        filtro_media(matriz, largo, ancho);
    }

    if (strcmp(tipo_filtro, "mediana") == 0)
    {
        filtro_mediana(matriz, largo, ancho);
    }
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

int main(int argc, char *argv[])
{
    // Parámetros MPI
    int myrank;
    int numero_procesos;
    const int ROOT = 0;
    const int ETIQUETA_COMUNICACIONES = 50;
    double tiempo_inicial, tiempo_final;
    clock_t tiempo_inicial_proceso_root, tiempo_final_proceso_root;

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numero_procesos);

    // Matriz de datos
    unsigned char **matriz_datos;

    // Datos matriz
    int largo_global = atoi(argv[4]);
    int ancho_global = atoi(argv[5]);

    // Filas para los procesos
    int filas_procesos = largo_global / numero_procesos;
    int resto_filas = largo_global % numero_procesos;
    
    // Captura de tiempo global
    MPI_Barrier(MPI_COMM_WORLD);
    tiempo_inicial = MPI_Wtime();

    if (myrank == ROOT)
    {
        // Comprobaciones iniciales
        if (argc != 6)
        {
            printf("Error en el numero de argumentos. Uso: <nombre_fichero_entrada> <nombre_fichero_salida> <nombre_filtro {sobel, media, mediana} > <largo> <ancho>\n");
            return 1;
        }

        if (comprobar_fichero_existe(argv[1]) == 0)
        {
            printf("El fichero de entrada %s no existe\n", argv[1]);
            return 1;
        }

        if (strcmp(argv[3], "sobel") == 0 || strcmp(argv[3], "media") == 0 || strcmp(argv[3], "mediana") == 0)
        {
            printf("Filtro elegido: %s\n", argv[3]);
            printf("Fichero de salida: %s\n", argv[2]);
        }
        else
        {
            printf("Error en el nombre del filtro. Uso: <nombre_fichero_entrada> <nombre_fichero_salida> <nombre_filtro {sobel, media, mediana} <largo> <ancho> >\n");
            return 1;
        }

        // Reservar memoria para la matriz del proceso 0
        reservar_memoria_matriz(&matriz_datos, largo_global, ancho_global);

        // Llenar la matriz
        leer_fichero_matriz(matriz_datos, argv[1], largo_global, ancho_global);

        tiempo_inicial_proceso_root = clock();
        // Distribuir la matriz
        for (int i = 1; i < numero_procesos; i++)
        {
            if (i < (numero_procesos - 1)) // Si no se trata del último proceso
            {
                for (int j = -1; j < filas_procesos + 1; j++)
                {
                    MPI_Send(matriz_datos[(filas_procesos * i) + j], ancho_global, MPI_UNSIGNED_CHAR, i, ETIQUETA_COMUNICACIONES, MPI_COMM_WORLD);
                }
            }
            else // Si es el último proceso
            {
                for (int j = -1; j < filas_procesos + resto_filas; j++)
                {
                    MPI_Send(matriz_datos[(filas_procesos * i) + j], ancho_global, MPI_UNSIGNED_CHAR, i, ETIQUETA_COMUNICACIONES, MPI_COMM_WORLD);
                }
            }
        }
    }
    else if (myrank < (numero_procesos - 1) && myrank > 0)
    {
        // Reservar memoria para la matriz para los procesos mayores que 0 y menores que N
        reservar_memoria_matriz(&matriz_datos, filas_procesos + 2, ancho_global);

        // Recibimos los datos para la matriz
        for (int i = 0; i < filas_procesos + 2; i++)
        {
            MPI_Recv(matriz_datos[i], ancho_global, MPI_UNSIGNED_CHAR, ROOT, ETIQUETA_COMUNICACIONES, MPI_COMM_WORLD, &status);
        }

        // Se aplica el filtro y se envían los datos de vuelta al proceso ROOT
        filtro(argv[3], matriz_datos, filas_procesos + 2, ancho_global);
        for (int i = 1; i < filas_procesos + 1; i++)
        {
            MPI_Send(matriz_datos[i], ancho_global, MPI_UNSIGNED_CHAR, ROOT, ETIQUETA_COMUNICACIONES, MPI_COMM_WORLD);
        }
    }
    else
    {
        // Reservar memoria para la matriz del último proceso
        reservar_memoria_matriz(&matriz_datos, filas_procesos + resto_filas + 1, ancho_global);

        // Recibimos los datos para la matriz
        for (int i = 0; i < filas_procesos + resto_filas + 1; i++)
        {
            MPI_Recv(matriz_datos[i], ancho_global, MPI_UNSIGNED_CHAR, ROOT, ETIQUETA_COMUNICACIONES, MPI_COMM_WORLD, &status);
        }

        // Se aplica el filtro y se envían los datos de vuelta al proceso ROOT
        filtro(argv[3], matriz_datos, filas_procesos + resto_filas + 1, ancho_global);
        for (int i = 1; i < filas_procesos + resto_filas + 1; i++)
        {
            MPI_Send(matriz_datos[i], ancho_global, MPI_UNSIGNED_CHAR, ROOT, ETIQUETA_COMUNICACIONES, MPI_COMM_WORLD);
        }
    }

    if (myrank == 0)
    {
        // Si sólo hay un proceso se llama al filtro con todas las filas, sino, se enviarán las filas correspondientes + 1
        (numero_procesos == 1) ? filtro(argv[3], matriz_datos, largo_global, ancho_global) : filtro(argv[3], matriz_datos, filas_procesos + 1, ancho_global);
        
        // Se recogen los datos de los procesos restantes
        for (int i = 1; i < numero_procesos; i++)
        {
            if (i < (numero_procesos - 1)) 
            {
                for (int j = 0; j < filas_procesos; j++)
                {
                    MPI_Recv(matriz_datos[(filas_procesos * i) + j], ancho_global, MPI_UNSIGNED_CHAR, i, ETIQUETA_COMUNICACIONES, MPI_COMM_WORLD, &status);
                }
            }
            else
            {
                for (int j = 0; j < filas_procesos + resto_filas; j++)
                {
                    MPI_Recv(matriz_datos[(filas_procesos * i) + j], ancho_global, MPI_UNSIGNED_CHAR, i, ETIQUETA_COMUNICACIONES, MPI_COMM_WORLD, &status);
                }
            }
        }
        tiempo_final_proceso_root = clock();

        volcar_fichero(matriz_datos, argv[2], largo_global, ancho_global);
    }

    // Liberar memoria
    for (int i = 0; i < filas_procesos; i++)
    {
        free(matriz_datos[i]);
    }

    free(matriz_datos);

    MPI_Barrier(MPI_COMM_WORLD);
    tiempo_final = MPI_Wtime();

    MPI_Finalize();

    if (myrank == ROOT)
    {
        FILE *resultados = fopen("resultados_ejecucion.txt", "w");
        printf("Tiempo de ejecucion global sin lectura/escritura <segundos>: %f\n", (float)(tiempo_final_proceso_root - tiempo_inicial_proceso_root) / CLOCKS_PER_SEC);
        printf("Tiempo de ejecucion global MPI <segundos>: %f\n", tiempo_final - tiempo_inicial);
        fprintf(resultados, "Tiempo de ejecucion global sin lectura/escritura <segundos>: %f\n", (float)(tiempo_final_proceso_root - tiempo_inicial_proceso_root) / CLOCKS_PER_SEC); 
        fprintf(resultados, "Tiempo de ejecucion global MPI <segundos>: %f\n", tiempo_final - tiempo_inicial);
        fprintf(resultados, "Ancho matriz: %d\n", atoi(argv[5]));
        fprintf(resultados, "Largo matriz: %d\n", atoi(argv[4]));
        fprintf(resultados, "Fichero de salida (imagen): %s\n", argv[2]);
        fprintf(resultados, "Fichero de entrada (imagen): %s\n", argv[1]);
        fprintf(resultados, "Filtro elegido: %s\n", argv[3]);
        fprintf(resultados, "Numero de procesos: %d\n", numero_procesos);
        fclose(resultados);
    }

    return 0;
}