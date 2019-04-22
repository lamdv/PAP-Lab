#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

int MatrixMultiply(int n, int *a, int *b, int *c);
void print_matrix(int n, int *a);

int main(int argc, char *argv[])
{
    int n = 4, *a, *b, *c, myrank;
    MPI_Init(&argc, &argv);
    a = (int *)malloc(n * n * sizeof(int));
    b = (int *)malloc(n * n * sizeof(int));
    c = (int *)malloc(n * n * sizeof(int));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            a[i * n + j] = i * n + j;
            if (i == j)
                b[i * n + j] = 1;
            else
                b[i * n + j] = 0;
        }
    }
    printf("A=\n");
    print_matrix(n, a);
    printf("B=\n");
    print_matrix(n, b);

    MatrixMultiply(n, a, b, c);
    MPI_Finalize();

    printf("====\n");
    printf("C=\n");
    print_matrix(n, c);
}

int MatrixMultiply(int n, int *a, int *b, int *c)
{
    int i, j, k;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < n; k++)
            {
                c[i * n + j] += a[i * n + k] * b[k * n + j];
            }
        }
    }
}

void print_matrix(int n, int *a)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%3d, ", a[i * n + j]);
        }
        printf("\n");
    }
}