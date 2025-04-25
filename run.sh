#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --job-name=PIvan
#SBATCH --comment="Run student mpi from config"
#SBATCH --output=complete.txt
mpirun ./a.out "$@"