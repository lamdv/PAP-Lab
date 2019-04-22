#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

/* Support Functions Declaration */
void Init(int n, int *a, int *b, int *c);
void MPIMatrixMultiply(int n, int *a, int *b, int *c, MPI_Comm comm);
int MatrixMultiply(int n, int *a, int *b, int *c);
void print_matrix(int n, int *a);

int main(int argc, char *argv[])
{
	int n = 8, *a, *b, *c, myrank;
	MPI_Init(&argc, &argv);
	a = (int *)malloc(n * n * sizeof(int));
	b = (int *)malloc(n * n * sizeof(int));
	c = (int *)malloc(n * n * sizeof(int));
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			a[i * n + j] = i * n + j;
			b[i * n + j] = 1;
		}
	}
	print_matrix(n, a);
	print_matrix(n, b);
	MPI_Barrier(MPI_COMM_WORLD);
	MPIMatrixMultiply(n, a, b, c, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();
	print_matrix(n, c);
}

void MPIMatrixMultiply(int n, int *a, int *b, int *c, MPI_Comm comm)
{
	int i;
	int nlocal;
	int size, dims[2], periods[2];
	int myrank, my2drank, mycoords[2];
	int up, down, left, right;
	// , coords[2];
	int shiftsource, shiftdest;
	MPI_Status status;
	MPI_Comm comm_2d;
	/* Get the communicator related information */
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &myrank);

	/* Set up the Cartesian topology */
	dims[0] = dims[1] = sqrt(size);

	/* Set the periods for wraparound connections */
	periods[0] = periods[1] = 1;

	/* Create the Cartesian topology */
	MPI_Cart_create(comm, 2, dims, periods, 1, &comm_2d);

	/* Get the rank and coordinates with respect to the Cartesian topology */
	MPI_Comm_rank(comm_2d, &my2drank);
	MPI_Cart_coords(comm_2d, my2drank, 2, mycoords);

	/* Compute ranks of the up and left shifts */
	MPI_Cart_shift(comm_2d, 0, -1, &right, &left);
	MPI_Cart_shift(comm_2d, 1, -1, &down, &up);
	nlocal = n / dims[0];

	/* Initial matrix skewing */
	MPI_Cart_shift(comm_2d, 0, -mycoords[0], &shiftsource, &shiftdest);

	// /* Shift matrix a left by one d */
	// MPI_Sendrecv_replace(a, nlocal * nlocal, MPI_INT, left, 1, right, 1, comm_2d, &status);

	// /* Shift matrix b up by one d */
	// MPI_Sendrecv_replace(b, nlocal * nlocal, MPI_INT, up, 1, down, 1, comm_2d, &status);
	MPI_Sendrecv_replace(a, nlocal * nlocal, MPI_INT, shiftdest, 1, shiftsource, 1, comm_2d, &status);

	MPI_Cart_shift(comm_2d, 1, -mycoords[1], &shiftsource, &shiftdest);
	MPI_Sendrecv_replace(b, nlocal * nlocal, MPI_INT, shiftdest, 1, shiftsource, 1, comm_2d, &status);

	/* Main computation loop */
	for (i = 0; i < dims[0]; i++)
	{
		/* Perform Serials multiplycation at location
            c=c+a*b*/
		MatrixMultiply(nlocal, a, b, c);

		/* Shift matrix a left by one d */
		MPI_Sendrecv_replace(a, nlocal * nlocal, MPI_INT, left, 1, right, 1, comm_2d, &status);

		/* Shift matrix b up by one d */
		MPI_Sendrecv_replace(b, nlocal * nlocal, MPI_INT, up, 1, down, 1, comm_2d, &status);
	}

	/* Restore a and b */
	MPI_Cart_shift(comm_2d, 0, +mycoords[0], &shiftsource, &shiftdest);
	MPI_Sendrecv_replace(a, nlocal * nlocal, MPI_INT, shiftdest, 1, shiftsource, 1, comm_2d, &status);

	MPI_Cart_shift(comm_2d, 1, +mycoords[1], &shiftsource, &shiftdest);
	MPI_Sendrecv_replace(b, nlocal * nlocal, MPI_INT, shiftdest, 1, shiftsource, 1, comm_2d, &status);

	MPI_Comm_free(&comm_2d);
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