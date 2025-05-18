@echo off
setlocal

:: Настройка среды Visual Studio — путь может отличаться!
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
if errorlevel 1 (
    echo Не удалось запустить vcvars64.bat. Проверь путь.
    pause
    exit /b 1
)
::_parallel
set filename=transport
:: Имя исходного файла (замените на ваше)
set SOURCE_FILE=%filename%.cpp
:: Имя выходного исполняемого файла
set OUTPUT_FILE=%filename%.exe
set OBJ_FILE=%filename%.obj

:: Компиляция с MSVC и поддержкой MPI
echo Compiling %SOURCE_FILE%...
cl /EHsc /std:c++20 /I"%MSMPI_INC%" "%SOURCE_FILE%" /Fe:"%OUTPUT_FILE%" /link /LIBPATH:"%MSMPI_LIB%" msmpi.lib
:: Удаление объектного файла
del /Q "%OBJ_FILE%"
echo Compilation successful!

endlocal