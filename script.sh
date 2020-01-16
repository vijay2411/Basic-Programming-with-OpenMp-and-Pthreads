export OMP_NUM_THREADS=$2
g++ $1 -o out.o -fopenmp -g
./out.o $3
