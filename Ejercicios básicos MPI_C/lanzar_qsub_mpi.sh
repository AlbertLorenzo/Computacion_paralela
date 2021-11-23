#!/bin/bash
#$ -N NOMBRE_CORTO
#$ -o Salida.out
#$ -e Error.err
#$ -q cp6.q
#$ -pe mpi 5
#$ -cwd
. /etc/profile.d/modules.sh
mpirun ./ej 
