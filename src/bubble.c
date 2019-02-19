#include <stdio.h>
#include <omp.h>

#include <x86intrin.h>

#define NBEXPERIMENTS 7
static long long unsigned int experiments[NBEXPERIMENTS];

/* 
  bubble sort -- sequential, parallel -- 
*/

static unsigned int N;

typedef int *array_int;

/* the X array will be sorted  */

static array_int X;

void swap(int *i, int *j)
{
  int t = *i;
  *i = *j;
  *j = t;
}

long long unsigned int average(long long unsigned int *exps)
{
  unsigned int i;
  long long unsigned int s = 0;

  for (i = 2; i < (NBEXPERIMENTS - 2); i++)
  {
    s = s + exps[i];
  }

  return s / (NBEXPERIMENTS - 2);
}

void init_array(array_int T)
{
  register int i;

  for (i = 0; i < N; i++)
  {
    T[i] = N - i;
    // i = 0;
  }
}

void print_array(array_int T)
{
  register int i;

  for (i = 0; i < N; i++)
  {
    printf("%d ", T[i]);
  }
  printf("\n");
}

/*
  test if T is properly sorted
 */
int is_sorted(array_int T)
{
  register int i;

  for (i = 1; i < N; i++)
  {
    /* test designed specifically for our usecase */
    if (T[i - 1] + 1 != T[i])
      return 0;
  }
  return 1;
}

void sequential_bubble_sort(int *T, const int size)
{
  /* TODO: sequential implementation of bubble sort */
  int i, is_swap = 0;
  while (is_swap == 0)
  {
    is_swap = 1;
    for (i = 0; i < size - 1; i++)
    {
      if (T[i] > T[i + 1])
      {
        swap(T + i, T + i + 1);
        is_swap = 0;
      }
    }
  }

  return;
}

void parallel_bubble_sort(int *T, const int size)
{
  /* TODO: parallel implementation of bubble sort */
  int i, is_swap = 0;
  while (is_swap == 0)
  {
    #pragma omp parallel default(none) shared(T, i, is_swap)
    is_swap = 1;
    for (i = 0; i < size - 1; i++)
    {
      if (T[i] > T[i + 1])
      {
        swap(T + i, T + i + 1);
        is_swap = 0;
      }
    }
  }
  return;
}

int main(int argc, char **argv)
{
  unsigned long long int start, end, residu;
  unsigned long long int av;
  unsigned int exp;

  /* the program takes one parameter N which is the size of the array to
     be sorted. The array will have size 2^N */
  if (argc != 2)
  {
    fprintf(stderr, "bubble N \n");
    exit(-1);
  }

  N = 1 << (atoi(argv[1]));
  X = (int *)malloc(N * sizeof(int));

  printf("--> Sorting an array of size %u\n", N);

  start = _rdtsc();
  end = _rdtsc();
  residu = end - start;

  printf("=======================\n");
  printf("\nSequential Bubble Sort: ");

  for (exp = 0; exp < NBEXPERIMENTS; exp++)
  {
    init_array(X);

    start = _rdtsc();

    sequential_bubble_sort(X, N);

    end = _rdtsc();
    experiments[exp] = end - start;

    /* verifying that X is properly sorted */
    if (!is_sorted(X))
    {
      fprintf(stderr, "ERROR: the array is not properly sorted\n");
      exit(-1);
    }
    // print_array(X);
  }

  av = average(experiments);

  printf("%Ld cycles\n\n", av - residu);

  // for (exp = 0; exp < NBEXPERIMENTS; exp++)
  // {
  //   init_array(X);
  //   start = _rdtsc();

  //   parallel_bubble_sort(X, N);

  //   end = _rdtsc();
  //   experiments[exp] = end - start;

  //   /* verifying that X is properly sorted */
  //   if (!is_sorted(X))
  //   {
  //     fprintf(stderr, "ERROR: the array is not properly sorted\n");
  //     exit(-1);
  //   }
  // }

  // av = average(experiments);
  // printf("\n bubble parallel \t %Ld cycles\n\n", av - residu);

  // print_array (X) ;
}
