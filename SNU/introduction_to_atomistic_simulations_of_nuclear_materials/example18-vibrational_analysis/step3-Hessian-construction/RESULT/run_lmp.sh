export PROJ222=step3_minus
mpirun -np 1 /home/bin/lammps_mpi  < input.${PROJ222} > output.${PROJ222}
export PROJ222=step3_plus
mpirun -np 1 /home/bin/lammps_mpi  < input.${PROJ222} > output.${PROJ222}

