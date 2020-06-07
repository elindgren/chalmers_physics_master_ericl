export PROJ222=step1
mpirun -np 1 /home/bin/lammps_mpi  < input.${PROJ222} > output.${PROJ222}
