#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <math.h>
#include <chrono>

// Algoritmo de ordenación básico, al ser siempre vectores de dimensión 9, la complejidad alta del algoritmo es despreciable
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

// Vuelca los datos de la matriz en un fichero .raw con formato binario
void dump_data(unsigned char **matrix, int width, int height, std::string output_file)
{
    std::ofstream output(output_file, std::ios::binary);

    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            output << matrix[i][j];
        }
    }

    output.close();
}

// Se declara y reserva memoria para la matriz
void init_matrix(unsigned char ***M_ptr, int width, int height)
{
    unsigned char **matrix = new unsigned char *[width];

    for (auto i = 0; i < height; ++i)
    {
        matrix[i] = new unsigned char[height];
    }

    *M_ptr = matrix;
}

// Función para inicializar y copiar los datos del fichero .raw en una matriz
void read_file(unsigned char ***matrix_ptr, std::string file_name, int width, int height)
{
    std::ifstream raw_data;
    raw_data.open(file_name, std::ios::binary);
    unsigned char **matrix = new unsigned char *[width];

    // Reserva espacio para la matriz
    for (auto i = 0; i < height; ++i)
    {
        matrix[i] = new unsigned char[height];
    }

    // Se añade la información de la imagen a la matriz
    while (!raw_data.eof())
    {
        for (int i = 0; i < width; i++)
        {
            for (int j = 0; j < height; j++)
            {
                raw_data >> matrix[i][j];
            }
        }
    }

    raw_data.close();

    *matrix_ptr = matrix;
}

// Copia los bordes de la matriz original en la resultante ya que estos no se alteran durante la operación de convolución
void copy_borders(unsigned char **matrix, unsigned char **result, int width, int height)
{
    for (int i = 0; i < width; i++)
    {
        result[0][i] = matrix[0][i];
    }

    for (int i = 0; i < height; i++)
    {
        result[i][0] = matrix[i][0];
    }

    for (int i = 0; i < width; i++)
    {
        result[width - 1][i] = matrix[width - 1][i];
    }

    for (int i = 0; i < height; i++)
    {
        result[i][height - 1] = matrix[i][height - 1];
    }
}

/* 
Filtros: Para aplicar filtros se utiliza una ventana o kernel que se desplaza por dentro de la matriz original con otra 3x3
por dentro de ésta última. Los bordes se copiarán en la matriz resultante/imagen previamente ya que a la hora de iterar se desprecian pero no para el cálculo 
de cada nuevo bit procesado.
*/

// Filtro de media
void average_filter(unsigned char **matrix, unsigned char **result, int width, int height)
{
    for (int i = 1; i < width - 1; i++)
    {
        for (int j = 1; j < height - 1; j++)
        {
            result[i][j] = (matrix[i - 1][j - 1] +
                            matrix[i - 1][j] +
                            matrix[i - 1][j + 1] +
                            matrix[i][j - 1] +
                            matrix[i][j] +
                            matrix[i][j + 1] +
                            matrix[i + 1][j - 1] +
                            matrix[i + 1][j] +
                            matrix[i + 1][j + 1]) /
                           9;
        }
    }
}

// Filtro para la mediana
void median_filter(unsigned char **matrix, unsigned char **result, int width, int height)
{
    int data[9];

    for (int i = 1; i < width - 1; i++)
    {
        for (int j = 1; j < height - 1; j++)
        {
            data[0] = matrix[i - 1][j - 1];
            data[1] = matrix[i - 1][j];
            data[2] = matrix[i - 1][j + 1];
            data[3] = matrix[i][j - 1];
            data[4] = matrix[i][j];
            data[5] = matrix[i][j + 1];
            data[6] = matrix[i + 1][j - 1];
            data[7] = matrix[i + 1][j];
            data[8] = matrix[i + 1][j + 1];

            // Se ordena el vector con los 9 datos y se escoge el del medio para el bit procesado
            bubble_sort(data, 9);

            result[i][j] = data[4];
        }
    }
}

// Operador sobel
void sobel_filter(unsigned char **matrix, unsigned char **result, int width, int height)
{

    int Gx[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
    int Gy[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
    unsigned char data[9];

    for (int i = 1; i < width - 1; i++)
    {
        for (int j = 1; j < height - 1; j++)
        {

            int gx_sum = 0, gy_sum = 0;

            data[0] = matrix[i - 1][j - 1];
            data[1] = matrix[i - 1][j];
            data[2] = matrix[i - 1][j + 1];
            data[3] = matrix[i][j - 1];
            data[4] = matrix[i][j];
            data[5] = matrix[i][j + 1];
            data[6] = matrix[i + 1][j - 1];
            data[7] = matrix[i + 1][j];
            data[8] = matrix[i + 1][j + 1];

            for (int k = 0; k < 9; k++)
            {
                gx_sum += data[k] * Gx[k];
                gy_sum += data[k] * Gy[k];
            }

            // Esta parte del código es opcional ya que sirve para reducir la intensidad de los bordes
            gx_sum /= 4;
            gy_sum /= 4;

            int processed_bit = sqrt(gx_sum * gx_sum + gy_sum * gy_sum);

            result[i][j] = processed_bit;
        }
    }
}

// Los parámetros recogidos son: anchura, altura, fichero de entrada, fichero de salida y filtro. No hay validación de errores.
int main(int argc, char *argv[])
{
    int width = std::stoi(argv[1]), height = std::stoi(argv[2]), filter = std::stoi(argv[5]);
    std::string in(argv[3]), out(argv[4]), used_filter;
    unsigned char **matrix, **result;
    
    // Inicio del programa
    read_file(&matrix, in, width, height);

    init_matrix(&result, width, height);

    copy_borders(matrix, result, width, height);

    auto begin = std::chrono::high_resolution_clock::now();
    switch (filter)
    {
    case 1:
        average_filter(matrix, result, width, height);
        used_filter = "Media";
        break;

    case 2:
        median_filter(matrix, result, width, height);
        used_filter = "Mediana";
        break;

    case 3:
        sobel_filter(matrix, result, width, height);
        used_filter = "Sobel";
        break;

    default:
        break;
    }
    auto end = std::chrono::high_resolution_clock::now();

    dump_data(result, width, height, out);
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);

    // Volcado información final con los resultados
    std::ofstream output("informacion_procesado_imagen.txt");
    output << "Parametros de ejecucion: " << "argv[0] = ejecutable, argv[1] = anchura, argv[2] = altura, argv[3] = fichero de entrada, argv[4] = fichero de salida, argv[5] = filtro seleccionado {1, 2, 3} \n"; 
    output << "Fichero de entrada: " << in << "\n";
    output << "Fichero de salida: " << out << "\n";
    output << "Filtro utilizado: " << used_filter << "\n";
    output << "Tiempo de ejecucion <microsegundos>: " << duration.count() << "\n";
    output.close();

    return 0;
}