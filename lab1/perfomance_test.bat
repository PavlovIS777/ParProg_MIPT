@echo off
setlocal enabledelayedexpansion

:: Имя исполняемого файла
set EXE=transport_parallel.exe
:: Лог-файл
set LOG_FILE=output.log

:: Очищаем лог, если он уже есть
if exist "%LOG_FILE%" del "%LOG_FILE%"

:: Параметры для запуска (просто список, без скобок)
set PARAM_1=50000 2000 1 1 1
set PARAM_2=50000 3000 1 1 1
set PARAM_3=50000 4000 1 1 1
set PARAM_4=50000 5000 1 1 1
set PARAM_5=50000 6000 1 1 1
set PARAM_6=50000 7000 1 1 1
set PARAM_7=50000 8000 1 1 1
set PARAM_8=50000 9000 1 1 1
set PARAM_9=50000 10000 1 1 1
set PARAM_10=50000 12000 1 1 1
set PARAM_11=50000 14000 1 1 1
set PARAM_12=50000 16000 1 1 1
set PARAM_13=50000 18000 1 1 1
set PARAM_14=50000 20000 1 1 1
set PARAM_15=50000 24000 1 1 1
set PARAM_16=50000 28000 1 1 1
set PARAM_17=50000 32000 1 1 1
set PARAM_18=50000 35000 1 1 1
set PARAM_19=50000 40000 1 1 1

:: Количество процессов MPI
set MPI_PROCESSES=4

:: Запускаем программу с каждым набором параметров
for /l %%i in (1,1,19) do (
    echo Running test with parameters: !PARAM_%%i!
    echo ===== Test with parameters: !PARAM_%%i! ===== >> "%LOG_FILE%"
    
    mpiexec -n !MPI_PROCESSES! "%EXE%" !PARAM_%%i! >> "%LOG_FILE%" 2>&1
    
    echo. >> "%LOG_FILE%"
    echo ===== End of test ===== >> "%LOG_FILE%"
    echo. >> "%LOG_FILE%"
)

echo All tests completed. Results saved to "%LOG_FILE%"
endlocal