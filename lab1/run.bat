@echo off
setlocal enabledelayedexpansion

set EXE=transport_parallel.exe
set PARAM=50000 12000 1 1 1 0
set MPI_PROCESSES=4

echo ===== Running with parameters: !PARAM! ===
mpiexec -n !MPI_PROCESSES! "%EXE%" !PARAM!
echo ===== End of run =====================================

endlocal