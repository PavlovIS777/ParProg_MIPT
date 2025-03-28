#!/bin/bash
#
#SBATCH --nodes=1 # Кол-во узлов
#SBATCH --ntasks-per-node=6 # Кол-во процессов на узел
#SBATCH --cpus-per-task=1 # Кол-во CPU на процесс (для многопоточности)
#SBATCH --partition=RT #группа машин кластера
#SBATCH --job-name=PIvan # Имя задачи (обязательно напишите что-нибудь своё)
#SBATCH --comment="Run student mpi from config" # Также обязательно к заполнению
#SBATCH --output=out.txt # Файл для печати вывода
#SBATCH --error=error.txt # Файл для печати ошибок
mpirun ./a.out "$@"