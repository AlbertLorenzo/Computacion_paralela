#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

void liberar_memoria_matriz(unsigned char **matriz, int dimension)
{
    for (int i = 0; i < dimension; ++i)
    {
        free(matriz[i]);
    }

    free(matriz);
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

// Si pasamos la dimensión por parámetro no se usará esta función
int obtener_dimension_fichero(FILE *fichero)
{
    fseek(fichero, 0, SEEK_END);
    long tam_fichero = ftell(fichero);
    fseek(fichero, 0, SEEK_SET);
    int dimension = sqrt(tam_fichero / sizeof(unsigned char));

    return dimension;
}

void intercambio(unsigned char *xp, unsigned char *yp)
{
    unsigned char temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void ordenar_burbuja(unsigned char arr[], int n)
{
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < n - i - 1; j++)
        {
            if (arr[j] > arr[j + 1])
            {
                intercambio(&arr[j], &arr[j + 1]);
            }
        }
    }
}

void reservar_memoria_matriz(unsigned char ***matriz_ptr, int dimension)
{
    // Se declara y reserva memoria para la matriz
    unsigned char **matriz = (unsigned char **)malloc(sizeof(unsigned char *) * dimension);

    if (matriz == NULL)
    {
        printf("ERROR: Fallo al reservar memoria para la matriz.");
        exit(1);
    }

    for (int i = 0; i < dimension; ++i)
    {
        matriz[i] = (unsigned char *)malloc(sizeof(unsigned char) * dimension);

        if (matriz[i] == NULL)
        {
            printf("ERROR: Fallo al reservar memoria para la matriz.");
            exit(1);
        }
    }

    *matriz_ptr = matriz;
}

void rellenar_matriz(unsigned char **matriz_datos, int dimension, FILE *fichero)
{
    for (int i = 0; i < dimension; i++)
    {
        fread(&matriz_datos[i][0], sizeof(unsigned char), dimension, fichero);
    }
}

void volcar_fichero(unsigned char **matriz, const char *nombre_fichero_salida, int dimension)
{
    FILE *fichero_salida = fopen(nombre_fichero_salida, "wb");
    for (int i = 0; i < dimension; i++)
    {
        fwrite(matriz[i], sizeof(unsigned char), dimension, fichero_salida);
    }
    fclose(fichero_salida);
}

void extension_simetrica(unsigned char **matriz_a, unsigned char **matriz_b, int dimension)
{
    // Extensión de filas
    for (int i = 1; i < dimension - 1; i++)
    {
        matriz_b[i][0] = matriz_a[i - 1][1];
        matriz_b[i][dimension - 1] = matriz_a[i - 1][dimension - 4];
    }

    // Extensión de columnas
    for (int i = 1; i < dimension - 1; i++)
    {
        matriz_b[0][i] = matriz_a[1][i - 1];
        matriz_b[dimension - 1][i] = matriz_a[dimension - 4][i - 1];
    }

    // Extensión de esquinas
    matriz_b[0][0] = matriz_b[2][2];
    matriz_b[0][dimension - 1] = matriz_b[2][dimension - 3];
    matriz_b[dimension - 1][0] = matriz_b[dimension - 3][2];
    matriz_b[dimension - 1][dimension - 1] = matriz_b[dimension - 3][dimension - 3];
}

void copiar_valores(unsigned char **matriz_a, unsigned char **matriz_b, int dimension)
{
    for (int i = 1; i < dimension - 1; i++)
    {
        /*
        * Se copian los valores de la matriz_a en la matriz_b. Se acceden a los 
        * valores de la matriz_b desde la posición [1][1] hasta la posición [n - 1][1]
        * copiando sólo n-2 valores.
        * 
        * Para la matriz b, empezamos desde la [0][0] hasta la [n-2][n-2]
        */
        memcpy(&matriz_b[i][1], &matriz_a[i - 1][0], sizeof(unsigned char) * (dimension - 2));
    }
}

void operador_sobel(unsigned char **matriz_a, unsigned char **matriz_b, int dimension)
{
    int c_sum, f_sum;
    int c[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
    int f[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
    unsigned char componentes[9];

    for (int i = 1; i < dimension - 1; i++)
    {
        for (int j = 1; j < dimension - 1; j++)
        {
            int c_sum = 0, f_sum = 0;

            // Se recogen los valores de los 9 pixeles que rodean al pixel actual
            componentes[0] = matriz_b[i - 1][j - 1];
            componentes[1] = matriz_b[i - 1][j];
            componentes[2] = matriz_b[i - 1][j + 1];
            componentes[3] = matriz_b[i][j - 1];
            componentes[4] = matriz_b[i][j];
            componentes[5] = matriz_b[i][j + 1];
            componentes[6] = matriz_b[i + 1][j - 1];
            componentes[7] = matriz_b[i + 1][j];
            componentes[8] = matriz_b[i + 1][j + 1];

            // Se calcula el valor de cada suma
            for (int k = 0; k < 9; k++)
            {
                c_sum += componentes[k] * c[k];
                f_sum += componentes[k] * f[k];
            }

            /*
            * Se cambia el valor por el nuevo procesado. 
            * Además, se accede a la matriz original con los índices - 1 ya que empezaríamos
            * desde el 0 0 en lugar del 1 1 del bucle
            */
            matriz_a[i - 1][j - 1] = sqrt(c_sum * c_sum + f_sum * f_sum);
        }
    }
}

void filtro_mediana(unsigned char **matriz_a, unsigned char **matriz_b, int dimension)
{
    unsigned char vector_aux[9];

    for (int i = 1; i < dimension - 1; i++)
    {
        for (int j = 1; j < dimension - 1; j++)
        {
            vector_aux[0] = matriz_a[i - 1][j - 1];
            vector_aux[1] = matriz_a[i - 1][j];
            vector_aux[2] = matriz_a[i - 1][j + 1];
            vector_aux[3] = matriz_a[i][j - 1];
            vector_aux[4] = matriz_a[i][j];
            vector_aux[5] = matriz_a[i][j + 1];
            vector_aux[6] = matriz_a[i + 1][j - 1];
            vector_aux[7] = matriz_a[i + 1][j];
            vector_aux[8] = matriz_a[i + 1][j + 1];

            ordenar_burbuja(vector_aux, 9);

            matriz_b[i][j] = vector_aux[4];
        }
    }
}

void filtro_media(unsigned char **matriz_a, unsigned char **matriz_b, int dimension)
{
    for (int i = 1; i < dimension - 1; i++)
    {
        for (int j = 1; j < dimension - 1; j++)
        {
            matriz_b[i][j] = (matriz_a[i - 1][j - 1] +
                              matriz_a[i - 1][j] +
                              matriz_a[i - 1][j + 1] +
                              matriz_a[i][j - 1] +
                              matriz_a[i][j] +
                              matriz_a[i][j + 1] +
                              matriz_a[i + 1][j - 1] +
                              matriz_a[i + 1][j] +
                              matriz_a[i + 1][j + 1]) /
                             9;
        }
    }
}

void copiar_bordes(unsigned char **matriz_a, unsigned char **matriz_b, int dimension)
{
    // Filas
    memcpy(matriz_b[0], matriz_a[0], sizeof(unsigned char) * dimension);
    memcpy(matriz_b[dimension - 1], matriz_a[dimension - 1], sizeof(unsigned char) * dimension);

    // Columnas
    for (int i = 0; i < dimension; i++)
    {
        matriz_b[i][0] = matriz_a[i][0];
        matriz_b[i][dimension - 1] = matriz_a[i][dimension - 1];
    }
}

int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        printf("Error en el numero de argumentos. Uso: <nombre_fichero_entrada> <nombre_fichero_salida> <nombre_filtro {sobel, media, mediana} > <dimension_imagen>\n");
        return 1;
    }

    if (comprobar_fichero_existe(argv[1]) == 0)
    {
        printf("El fichero %s no existe\n", argv[1]);
        return 1;
    }

    if (strcmp(argv[3], "sobel") == 0 || strcmp(argv[3], "media") == 0 || strcmp(argv[3], "mediana") == 0)
    {
        printf("Filtro elegido: %s\n", argv[3]);
        printf("Fichero de salida: %s\n", argv[2]);
    }
    else
    {
        printf("Error en el nombre del filtro. Uso: <nombre_fichero_entrada> <nombre_fichero_salida> <nombre_filtro {sobel, media, mediana} >\n");
        return 1;
    }

    unsigned char **matriz_a, **matriz_b;

    FILE *fichero = fopen(argv[1], "rb");

    const int DIM = atoi(argv[4]);

    if (DIM == 0)
    {
        printf("Error al obtener la dimensión del fichero\n");
        return 1;
    }

    clock_t inicio, fin;
    inicio = clock();
    if (strcmp(argv[3], "media") == 0 || strcmp(argv[3], "mediana") == 0)
    {
        /*
         * Reservamos memoria para las dos matrices
         */
        reservar_memoria_matriz(&matriz_a, DIM);
        reservar_memoria_matriz(&matriz_b, DIM);

        /*
         * Leemos la matriz del fichero
         */
        rellenar_matriz(matriz_a, DIM, fichero);

        /*
         * Copiamos los bordes de la matriz_a en la matriz_b
         */
        copiar_bordes(matriz_a, matriz_b, DIM);

        /*
         * Aplicamos el filtro
         */
        if (strcmp(argv[3], "media") == 0)
        {
            filtro_media(matriz_a, matriz_b, DIM);
        }
        else
        {
            filtro_mediana(matriz_a, matriz_b, DIM);
        }

        /*
         * Se vuelca la matriz_b en el fichero de salida
         */
        volcar_fichero(matriz_b, argv[2], DIM);

        liberar_memoria_matriz(matriz_a, DIM);
        liberar_memoria_matriz(matriz_b, DIM);
    }


    /*
    * Filtro sobel
    *
    * Para el proceso sobel será un poco diferente, ya que primero rellenaremos la matriz_a con los 
    * datos del fichero, para posteriormente rellenar la matriz_b con los datos de la matriz_a y hacer
    * la extensión simétrica.
    * 
    * Finalmente se aplicará el filtro y se grabarán los datos sobre el fichero de salida.
    */
    if (strcmp(argv[3], "sobel") == 0)
    {
        reservar_memoria_matriz(&matriz_a, DIM);
        reservar_memoria_matriz(&matriz_b, DIM + 2);

        rellenar_matriz(matriz_a, DIM, fichero);
        copiar_valores(matriz_a, matriz_b, DIM + 2);
        extension_simetrica(matriz_a, matriz_b, DIM + 2);
        operador_sobel(matriz_a, matriz_b, DIM + 2);
        volcar_fichero(matriz_a, argv[2], DIM);

        liberar_memoria_matriz(matriz_a, DIM);
        liberar_memoria_matriz(matriz_b, DIM + 2);

        reservar_memoria_matriz(&matriz_a, DIM);
        reservar_memoria_matriz(&matriz_b, DIM + 2);
    }

    fin = clock();
    FILE *fichero_resultados = fopen("resultados.txt", "w");
    fprintf(fichero_resultados, "Filtro seleccionado: %s\n", argv[3]);
    fprintf(fichero_resultados, "Dimension de la imagen: %d\n", DIM);
    fprintf(fichero_resultados, "Fichero de salida: %s\n", argv[2]);
    fprintf(fichero_resultados, "Tiempo de ejecución: %f\n", (double)(fin - inicio) / CLOCKS_PER_SEC);

    fclose(fichero_resultados);
    fclose(fichero);

    return 0;
}