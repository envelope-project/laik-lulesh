module load openmpi
make clean
make TARGET=REPARTITIONING -j
echo "running REPARTITIONING"
mpirun --oversubscribe -host i10se1 -n 13 ./lulesh2.0 -s 30 -repart 0 : -host i10se4 -n 14 ./lulesh2.0 -s 30 -repart 0 > repart.log

make clean
make TARGET=PERFORMANCE -j
echo "running PERFORMANCE"
mpirun --oversubscribe -host i10se1 -n 13 ./lulesh2.0 -s 30 -repart 0 : -host i10se4 -n 14 ./lulesh2.0 -s 30 -repart 0 > perf.log
