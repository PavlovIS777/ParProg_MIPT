@echo off
setlocal enabledelayedexpansion

set EXE=transport.exe
set PARAM=5000 1200 1 1 1 1
set MPI_PROCESSES=4

echo ===== Running with parameters: !PARAM! ===
mpiexec -n !MPI_PROCESSES! "%EXE%" !PARAM!
echo ===== End of run =====================================

endlocal