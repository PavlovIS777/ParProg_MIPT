@echo off
setlocal

:: ��������� ����� Visual Studio � ���� ����� ����������!
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
if errorlevel 1 (
    echo �� ������� ��������� vcvars64.bat. ������� ����.
    pause
    exit /b 1
)
::_parallel
set filename=transport
:: ��� ��������� ����� (�������� �� ����)
set SOURCE_FILE=%filename%.cpp
:: ��� ��������� ������������ �����
set OUTPUT_FILE=%filename%.exe
set OBJ_FILE=%filename%.obj

:: ���������� � MSVC � ���������� MPI
echo Compiling %SOURCE_FILE%...
cl /EHsc /std:c++20 /I"%MSMPI_INC%" "%SOURCE_FILE%" /Fe:"%OUTPUT_FILE%" /link /LIBPATH:"%MSMPI_LIB%" msmpi.lib
:: �������� ���������� �����
del /Q "%OBJ_FILE%"
echo Compilation successful!

endlocal