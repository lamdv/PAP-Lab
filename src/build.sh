

if [[ $* == *"1"* ]] || [ $# == 0 ];
then
echo "=================================================================="
echo "Building Task 1: OpenMP Loop Scheduling with Matrix Multiplication"
echo "gcc -fopenmp triangular_matrix.c -o build/triangular_matrix"
gcc -fopenmp triangular_matrix.c -o build/triangular_matrix 
./build/triangular_matrix
fi
if [[ $* == *"2"* ]] || [ $# == 0 ];
then
echo "=================================================================="
echo "Building Task 2: Bubble Sort - Sequential and Parallel Implementation"
echo "gcc -fopenmp bubble.c -o build/bubble"
gcc -fopenmp bubble.c -o build/bubble 
build/bubble 10
fi
if [[ $* == *"3"* ]] || [ $# == 0 ];
then
echo "=================================================================="
echo "Building Task 3: Quick Sort"
echo "gcc -fopenmp qsort.c -o build/qsort"
gcc -fopenmp qsort.c -o build/qsort
build/qsort
fi
if [[ $* == *"4"* ]] || [ $# == 0 ];
then
echo "=================================================================="
echo "Building Task 2: Merge Sort"
echo "gcc -fopenmp mergesort.c -o build/mergesort"
gcc -fopenmp mergesort.c -o build/mergesort
build/mergesort
fi
echo "=================================================================="
