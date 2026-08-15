@echo off
setlocal

set "SCRIPT_DIR=%~dp0"
for %%I in ("%SCRIPT_DIR%..") do set "ROOT=%%~fI"
set BUILD_DIR=%ROOT%\build
set VCVARS=D:\Program Files\Microsoft Visual Studio\18\Insiders\VC\Auxiliary\Build\vcvars64.bat
set CMAKE_EXE=D:\Program Files\CMake\bin\cmake.exe

if not exist "%BUILD_DIR%" mkdir "%BUILD_DIR%"
cd /d "%BUILD_DIR%"

if exist CMakeCache.txt del /f /q CMakeCache.txt
if exist CMakeFiles rmdir /s /q CMakeFiles
if exist decom_h.pyd del /f /q decom_h.pyd

call "%VCVARS%"
if errorlevel 1 goto :err

"%CMAKE_EXE%" -S "%ROOT%\csrc" -B "%BUILD_DIR%" -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=D:/vcpkg-master/installed/x64-windows -DPython3_EXECUTABLE=D:/Anaconda3/python.exe -DPYTHON_EXECUTABLE=D:/Anaconda3/python.exe
if errorlevel 1 goto :err

"%CMAKE_EXE%" --build .
if errorlevel 1 goto :err

echo.
echo Build succeeded.
exit /b 0

:err
echo.
echo Build failed with code %errorlevel%.
exit /b %errorlevel%
