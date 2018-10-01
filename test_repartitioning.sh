echo "testing 8, no repartitioning"
mpirun -np 8 --oversubscribe ./lulesh2.0 -q -s 2 -repart 0
echo "testing 8 -> 1, @ 10"
mpirun -np 8 --oversubscribe ./lulesh2.0 -q -s 2 -repart 1 -repart_cycle 10
echo "testing 8 -> 1, @ 20"
mpirun -np 8 --oversubscribe ./lulesh2.0 -q -s 2 -repart 1 -repart_cycle 20
echo "testing 27 -> 8, @ 20"
mpirun -np 27 --oversubscribe ./lulesh2.0 -q -s 2 -repart 8 -repart_cycle 20
echo "testing 27 -> 1, @ 20"
mpirun -np 27 --oversubscribe ./lulesh2.0 -q -s 2 -repart 1 -repart_cycle 20
echo "testing 64 -> 8, @ 20"
mpirun -np 64 --oversubscribe ./lulesh2.0 -q -s 2 -repart 8 -repart_cycle 20

