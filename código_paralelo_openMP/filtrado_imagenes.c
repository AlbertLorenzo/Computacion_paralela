#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

void reservar_memoria_matriz(unsigned char ***matriz_ptr, int largo, int ancho);
void liberar_memoria_matriz(unsigned char **matriz_ptr, int largo);

int comprobar_fichero_existe(char *fichero);
void leer_fichero_matriz(unsigned char **matriz, const char *nombre_fichero, int largo, int ancho);
void volcar_fichero(unsigned char **matriz, const char *nombre_fichero_salida, int largo, int ancho);

void procesado_imagen(const char *tipo_filtro, unsigned char **matriz, unsigned char **matriz_filtrada, int inicio, int fin, int ancho, int largo);

void filtro_media_nuevo(unsigned char **matriz, unsigned char **matriz_filtrada, int inicio, int fin, int ancho, int largo);
void filtro_mediana_nuevo(unsigned char **matriz, unsigned char **matriz_filtrada, int inicio, int fin, int ancho, int largo);
void filtro_sobel_nuevo(unsigned char **matriz, unsigned char **matriz_filtrada, int inicio, int fin, int ancho, int largo);
void extension_simetrica(unsigned char **M, unsigned char *componentes, int i, int j, int largo, int ancho);

int main(int argc, char *argv[])
{
    // Comprobaciones iniciales
    if (argc != 6)
    {
        printf("Error en el numero de argumentos. Uso: <nombre_fichero_entrada> <nombre_fichero_salida> <nombre_filtro {sobel, media, mediana} > <dimension_imagen> <numero_hilos>\n");
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
        printf("Error en el nombre del filtro. Uso: <nombre_fichero_entrada> <nombre_fichero_salida> <nombre_filtro {sobel, media, mediana} <dimension_imagen> <numero_hilos>\n");
        return 1;
    }

    int dimension_imagen = atoi(argv[4]);

    if (dimension_imagen <= 0)
    {
        printf("Error en la dimension de la imagen. Uso: <nombre_fichero_entrada> <nombre_fichero_salida> <nombre_filtro {sobel, media, mediana} <dimension_imagen> <numero_hilos>\n");
        return 1;
    }

    // Parámetros para medir el tiempo
    double inicio_openmp, fin_openmp, tiempo_openmp, inicio_lectura, fin_lectura;

    // Parámetros openMP
    int id_hilo, numero_hilos;

    // Número de hilos
    numero_hilos = atoi(argv[5]);

    // Número de largo por hilo
    int filas_hilo = dimension_imagen / numero_hilos;
    int resto_filas = dimension_imagen % numero_hilos;

    // Reserva memoria para la matriz
    unsigned char **matriz;
    unsigned char **matriz_filtrada;
    reservar_memoria_matriz(&matriz_filtrada, dimension_imagen, dimension_imagen);
    reservar_memoria_matriz(&matriz, dimension_imagen, dimension_imagen);

    inicio_lectura = omp_get_wtime();
    leer_fichero_matriz(matriz, argv[1], dimension_imagen, dimension_imagen);
    fin_lectura = omp_get_wtime();

    inicio_openmp = omp_get_wtime();

#pragma omp parallel num_threads(numero_hilos) private(id_hilo) shared(filas_hilo, resto_filas, dimension_imagen, matriz, matriz_filtrada)
    {
        id_hilo = omp_get_thread_num();

        // Se calcula el índice por el que empieza y acaba cada hilo, el hilo master se encarga de las filas restantes
        int inicio = (id_hilo == 0) ? 0 : filas_hilo * id_hilo + resto_filas;
        int fin = (id_hilo == 0) ? filas_hilo + resto_filas : filas_hilo * (id_hilo + 1) + resto_filas;

        procesado_imagen(argv[3], matriz, matriz_filtrada, inicio, fin, dimension_imagen, dimension_imagen);
    }

    fin_openmp = omp_get_wtime();

    volcar_fichero(matriz_filtrada, argv[2], dimension_imagen, dimension_imagen);

    FILE *fichero_resultados = fopen("resultados.txt", "w");   

    tiempo_openmp = fin_openmp - inicio_openmp;
    fprintf(fichero_resultados, "Tiempo lectura: %f segundos \n", fin_lectura - inicio_lectura);
    fprintf(fichero_resultados, "Tiempo openMP: %f segundos \n", tiempo_openmp);
    fprintf(fichero_resultados, "Hilos usados: %d \n", numero_hilos);
    fprintf(fichero_resultados, "Dimension imagen: %d \n", dimension_imagen);
    fprintf(fichero_resultados, "Filtro: %s\n", argv[3]);
    fprintf(fichero_resultados, "Nombre fichero entrada: %s", argv[1]);

    fclose(fichero_resultados);
    liberar_memoria_matriz(matriz, dimension_imagen);
    liberar_memoria_matriz(matriz_filtrada, dimension_imagen);

    return 0;
}

void filtro_media_nuevo(unsigned char **matriz, unsigned char **matriz_filtrada, int inicio, int fin, int ancho, int largo)
{
    unsigned char componentes[9];

    for (int i = inicio; i < fin; i++)
    {
        for (int j = 0; j < ancho; j++)
        {
            // Cuando se llega a los bordes de la imagen, se copia el valor de la matriz original
            if (i == 0 || j == 0 || i == largo - 1 || j == ancho - 1)
            {
                matriz_filtrada[i][j] = matriz[i][j];
            }
            else
            {
                // En caso contrario, se calcula la media de los 9 componentes
                matriz_filtrada[i][j] = (matriz[i - 1][j - 1] + matriz[i - 1][j] + matriz[i - 1][j + 1] + matriz[i][j - 1] + matriz[i][j] + matriz[i][j + 1] + matriz[i + 1][j - 1] + matriz[i + 1][j] + matriz[i + 1][j + 1]) / 9;
            }
        }
    }
}

void filtro_mediana_nuevo(unsigned char **matriz, unsigned char **matriz_filtrada, int inicio, int fin, int ancho, int largo)
{
    unsigned char componentes[9];

    for (int i = inicio; i < fin; i++)
    {
        for (int j = 0; j < ancho; j++)
        {
            if (i == 0 || j == 0 || i == largo - 1 || j == ancho - 1)
            {
                matriz_filtrada[i][j] = matriz[i][j];
            }
            else
            {
                // Guarda las componentes de la matriz en un vector
                componentes[0] = matriz[i - 1][j - 1];
                componentes[1] = matriz[i - 1][j];
                componentes[2] = matriz[i - 1][j + 1];
                componentes[3] = matriz[i][j - 1];
                componentes[4] = matriz[i][j];
                componentes[5] = matriz[i][j + 1];
                componentes[6] = matriz[i + 1][j - 1];
                componentes[7] = matriz[i + 1][j];
                componentes[8] = matriz[i + 1][j + 1];

                // Ordena las componentes del vector
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

                matriz_filtrada[i][j] = componentes[4];
            }
        }
    }
}

void filtro_sobel_nuevo(unsigned char **matriz, unsigned char **matriz_filtrada, int inicio, int fin, int ancho, int largo)
{
    int c[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
    int f[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
    unsigned char componentes[9];

    for (int i = inicio; i < fin; i++)
    {
        for (int j = 0; j < ancho; j++)
        {
            if (i == 0 || j == 0 || i == largo - 1 || j == ancho - 1)
            {
                // Se calculan los valores de los componentes de la matriz teniendo en cuenta la extensión simétrica cuando iteramos sobre los bordes
                extension_simetrica(matriz, componentes, i, j, largo, ancho);
            }
            else
            {
                // En caso contrario, se guardan los valores como de costumbre
                componentes[0] = matriz[i - 1][j - 1];
                componentes[1] = matriz[i - 1][j];
                componentes[2] = matriz[i - 1][j + 1];
                componentes[3] = matriz[i][j - 1];
                componentes[4] = matriz[i][j];
                componentes[5] = matriz[i][j + 1];
                componentes[6] = matriz[i + 1][j - 1];
                componentes[7] = matriz[i + 1][j];
                componentes[8] = matriz[i + 1][j + 1];
            }

            int suma_c = 0;
            int suma_f = 0;

            // Se calcula el valor de cada suma
            for (int k = 0; k < 9; k++)
            {
                suma_c += c[k] * componentes[k];
                suma_f += f[k] * componentes[k];
            }

            // Se calcula el valor del pixel actual
            matriz_filtrada[i][j] = sqrt(suma_c * suma_c + suma_f * suma_f);
        }
    }
}

void extension_simetrica(unsigned char **matriz, unsigned char *componentes, int i, int j, int largo, int ancho)
{
    if (i == 0 && j == 0)
    {
        // Esquina superior izquierda
        componentes[0] = matriz[i + 1][j + 1];
        componentes[1] = matriz[i + 1][j];
        componentes[2] = matriz[i + 1][j + 1];

        componentes[3] = matriz[i][j + 1];
        componentes[4] = matriz[i][j];
        componentes[5] = matriz[i][j + 1];

        componentes[6] = matriz[i + 1][j + 1];
        componentes[7] = matriz[i + 1][j];
        componentes[8] = matriz[i + 1][j + 1];
    }
    else if (i == 0 && j == ancho - 1)
    {
        // Esquina superior derecha
        componentes[0] = matriz[i + 1][j - 1];
        componentes[1] = matriz[i + 1][j];
        componentes[2] = matriz[i + 1][j - 1];

        componentes[3] = matriz[i][j - 1];
        componentes[4] = matriz[i][j];
        componentes[5] = matriz[i][j - 1];

        componentes[6] = matriz[i + 1][j - 1];
        componentes[7] = matriz[i + 1][j];
        componentes[8] = matriz[i + 1][j - 1];
    }
    else if (i == largo - 1 && j == 0)
    {
        // Esquina inferior izquierda
        componentes[0] = matriz[i - 1][j + 1];
        componentes[1] = matriz[i - 1][j];
        componentes[2] = matriz[i - 1][j + 1];

        componentes[3] = matriz[i][j + 1];
        componentes[4] = matriz[i][j];
        componentes[5] = matriz[i][j + 1];

        componentes[6] = matriz[i - 1][j + 1];
        componentes[7] = matriz[i - 1][j];
        componentes[8] = matriz[i - 1][j + 1];
    }
    else if (i == largo - 1 && j == ancho - 1)
    {
        // Esquina inferior derecha
        componentes[0] = matriz[i - 1][j - 1];
        componentes[1] = matriz[i - 1][j];
        componentes[2] = matriz[i - 1][j - 1];

        componentes[3] = matriz[i][j - 1];
        componentes[4] = matriz[i][j];
        componentes[5] = matriz[i][j - 1];

        componentes[6] = matriz[i - 1][j - 1];
        componentes[7] = matriz[i - 1][j];
        componentes[8] = matriz[i - 1][j - 1];
    }
    else if (i == 0)
    {
        // Borde superior
        componentes[0] = matriz[i + 1][j - 1];
        componentes[1] = matriz[i + 1][j];
        componentes[2] = matriz[i + 1][j + 1];

        componentes[3] = matriz[i][j - 1];
        componentes[4] = matriz[i][j];
        componentes[5] = matriz[i][j + 1];

        componentes[6] = matriz[i + 1][j - 1];
        componentes[7] = matriz[i + 1][j];
        componentes[8] = matriz[i + 1][j + 1];
    }
    else if (j == 0)
    {
        // Borde izquierdo
        componentes[0] = matriz[i - 1][j + 1];
        componentes[1] = matriz[i - 1][j];
        componentes[2] = matriz[i - 1][j + 1];

        componentes[3] = matriz[i][j + 1];
        componentes[4] = matriz[i][j];
        componentes[5] = matriz[i][j + 1];

        componentes[6] = matriz[i + 1][j + 1];
        componentes[7] = matriz[i + 1][j];
        componentes[8] = matriz[i + 1][j + 1];
    }
    else if (i == largo - 1 && j != 0 && j != ancho - 1)
    {
        // Borde inferior
        componentes[0] = matriz[i - 1][j - 1];
        componentes[1] = matriz[i - 1][j];
        componentes[2] = matriz[i - 1][j + 1];

        componentes[3] = matriz[i][j - 1];
        componentes[4] = matriz[i][j];
        componentes[5] = matriz[i][j + 1];

        componentes[6] = matriz[i - 1][j - 1];
        componentes[7] = matriz[i - 1][j];
        componentes[8] = matriz[i - 1][j + 1];
    }
    else if (j == ancho - 1)
    {
        // Borde derecho
        componentes[0] = matriz[i - 1][j - 1];
        componentes[1] = matriz[i - 1][j];
        componentes[2] = matriz[i - 1][j - 1];

        componentes[3] = matriz[i][j - 1];
        componentes[4] = matriz[i][j];
        componentes[5] = matriz[i][j - 1];

        componentes[6] = matriz[i + 1][j - 1];
        componentes[7] = matriz[i + 1][j];
        componentes[8] = matriz[i + 1][j - 1];
    }
}

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

void liberar_memoria_matriz(unsigned char **matriz_ptr, int largo)
{
    for (int i = 0; i < largo; i++)
    {
        free(matriz_ptr[i]);
    }
    free(matriz_ptr);
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

void procesado_imagen(const char *tipo_filtro, unsigned char **matriz, unsigned char **matriz_filtrada, int inicio, int fin, int ancho, int largo)
{
    if (strcmp(tipo_filtro, "sobel") == 0)
    {
        filtro_sobel_nuevo(matriz, matriz_filtrada, inicio, fin, ancho, largo);
    }

    if (strcmp(tipo_filtro, "media") == 0)
    {
        filtro_media_nuevo(matriz, matriz_filtrada, inicio, fin, ancho, largo);
    }

    if (strcmp(tipo_filtro, "mediana") == 0)
    {
        filtro_mediana_nuevo(matriz, matriz_filtrada, inicio, fin, ancho, largo);
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