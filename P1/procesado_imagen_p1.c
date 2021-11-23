#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void bubble_sort(int array[], int size)
{
    for (int i = 0; i < size - 1; i++)
    {
        for (int j = 0; j < size - i - 1; j++)
        {
            if (array[j] > array[i + 1])
            {
                int aux = array[j];
                array[j] = array[i];
                array[i] = aux;
            }
        }
    }
}

// Se declara y reserva memoria cualquier matriz
void init_matrix(unsigned char ***M_ptr, int width, int height)
{
    unsigned char **matrix;
    matrix = (unsigned char **)malloc(sizeof(unsigned char *) * height);

    for (int i = 0; i < width; i++)
    {
        matrix[i] = (unsigned char *)malloc(sizeof(unsigned char) * height);
    }

    *M_ptr = matrix;
}

// Función para inicializar y copiar los datos del fichero .raw en una matriz
void read_binary_data(unsigned char **matrix, const char *input_file, int width, int height)
{
    FILE *raw_data;

    raw_data = fopen(input_file, "rb");

    for (int i = 1; i < height + 1; ++i)
    {
        fread(matrix[i] + 1, sizeof(unsigned char), width, raw_data);
    }

    fclose(raw_data);
}

// Vuelca los datos en un fichero .raw con formato binario
void write_binary_data(unsigned char *array, const char *output_file, int width, int height)
{
    FILE *output;

    output = fopen(output_file, "wb");

    fwrite(array, sizeof(unsigned char), width * height, output);

    fclose(output);
}

void padding(unsigned char **matrix, int width, int height)
{
    // Esquina izquierda superior
    matrix[0][0] = matrix[2][2];

    // Esquina derecha superior
    matrix[0][width + 1] = matrix[2][width - 1];

    // Esquina izquierda inferior
    matrix[height + 1][0] = matrix[height - 1][2];

    // Esquina derecha inferior
    matrix[height + 1][width + 1] = matrix[height - 1][width - 1];

    // Fila superior
    for (int i = 1; i <= width; i++)
    {
        matrix[0][i] = matrix[2][i];
    }

    // Fila inferior
    for (int i = 1; i <= width; i++)
    {
        matrix[height + 1][i] = matrix[height - 1][i];
    }

    // Columna izquierda
    for (int i = 1; i <= height; i++)
    {
        matrix[i][0] = matrix[i][2];
    }

    // Columna derecha
    for (int i = 0; i <= height; i++)
    {
        matrix[i][width + 1] = matrix[i][width - 1];
    }
}

// Filtro de media
void average_filter(unsigned char **matrix, unsigned char **result, int width, int height)
{
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            result[i][j] = (matrix[i][j] +
                            matrix[i][j + 1] +
                            matrix[i][j + 2] +
                            matrix[i + 1][j] +
                            matrix[i + 1][j + 1] +
                            matrix[i + 1][j + 2] +
                            matrix[i + 2][j] +
                            matrix[i + 2][j + 1] +
                            matrix[i + 2][j + 2]) /
                           9;
        }
    }
}

// Filtro para la mediana
void median_filter(unsigned char **matrix, unsigned char **result, int width, int height)
{
    int data[9];

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            data[0] = matrix[i][j];
            data[1] = matrix[i][j + 1];
            data[2] = matrix[i][j + 2];
            data[3] = matrix[i + 1][j];
            data[4] = matrix[i + 1][j + 1];
            data[5] = matrix[i + 1][j + 2];
            data[6] = matrix[i + 2][j];
            data[7] = matrix[i + 2][j + 1];
            data[8] = matrix[i + 2][j + 2];

            // Se ordena el vector con los 9 datos y se escoge el del medio para el bit procesado
            bubble_sort(data, 9);

            result[i][j] = data[4];
        }
    }
}

void sobel_filter(unsigned char **matrix, unsigned char *result, int width, int height)
{
    int Gx[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1}, Gy[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
    unsigned char data[9];

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int gx_sum = 0, gy_sum = 0;

            data[0] = matrix[i][j];
            data[1] = matrix[i][j + 1];
            data[2] = matrix[i][j + 2];
            data[3] = matrix[i + 1][j];
            data[4] = matrix[i + 1][j + 1];
            data[5] = matrix[i + 1][j + 2];
            data[6] = matrix[i + 2][j];
            data[7] = matrix[i + 2][j + 1];
            data[8] = matrix[i + 2][j + 2];

            for (int k = 0; k < 9; k++)
            {
                gx_sum += data[k] * Gx[k];
                gy_sum += data[k] * Gy[k];
            }

            // Esta parte del código es opcional ya que sirve para reducir la intensidad de los bordes
            gx_sum /= 4;
            gy_sum /= 4;

            int processed_bit = sqrt(gx_sum * gx_sum + gy_sum * gy_sum);

            result[i * width + j] = processed_bit;
        }
    }
}

// Los parámetros recogidos son: anchura, altura, fichero de entrada, fichero de salida y filtro. No hay validación de errores.
int main(int argc, char *argv[])
{
    int width = atoi(argv[1]), height = atoi(argv[2]), filter = atoi(argv[5]);
    char *input_file = argv[3], *output_file = argv[4];

    unsigned char **matrix, result[width * height];

    init_matrix(&matrix, width + 2, height + 2);

    read_binary_data(matrix, input_file, width, height);

    padding(matrix, width, height);

    sobel_filter(matrix, result, width, height);

    write_binary_data(result, output_file, width, height);

    return 0;
}