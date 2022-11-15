#include <stdio.h>
#include <stdlib.h>

#define DIM 512

void read_file(const char *file_name, unsigned char **m)
{
  FILE *file = fopen(file_name, "rb");
  for (int i = 0; i < DIM; i++)
  {
    fread(&m[i][0], sizeof(unsigned char), DIM, file);
  }
  fclose(file);
}

void allocate_matrix_memory(unsigned char ***m_ptr)
{
  unsigned char **m = (unsigned char **)malloc(sizeof(unsigned char *) * DIM);
  for (int i = 0; i < DIM; i++)
  {
    m[i] = (unsigned char *)malloc(sizeof(unsigned char) * DIM);
  }
  *m_ptr = m;
}

void compare_matrix(unsigned char **mA, unsigned char **mB)
{
  int fallos = 0;
  for (int i = 0; i < DIM; i++)
  {
    for (int j = 0; j < DIM; j++)
    {
      if (mA[i][j] != mB[i][j])
      {
        fallos++;
        printf("p: [%d][%d] vmA: %d, vmH %d\n", i, j, mA[i][j], mB[i][j]);
      }
    }
  }
  printf("Fallos: %d\n", fallos);
}

int main(void)
{
  unsigned char **mA, **mH;
  allocate_matrix_memory(&mA);
  allocate_matrix_memory(&mH);

  read_file("matriz_loquesea_a.raw", mA);
  read_file("matriz_loquesea_b.raw", mH);

  compare_matrix(mA, mH);

  return 0;
}
