@echo off
rem examples/windows/planetgen-orbits.cmd
rem
rem Windows counterpart of examples/maintenance/planetgen-orbits@.service:
rem runs planetgen.cli.orbits once against one database. Schedule it
rem monthly with Task Scheduler (docs/deployment/windows.md), one task
rem per database:
rem     planetgen-orbits.cmd planetgen
rem Credentials come from config.json like everything else.

if "%~1"=="" (
    echo usage: %~nx0 DATABASE
    exit /b 2
)
set "PLANETGEN_MYSQL_DATABASE=%~1"
set "NLTK_DATA=C:\ProgramData\planetgen\nltk_data"
cd /d C:\srv\planetGen\src
"C:\srv\planetgen-venv\Scripts\python.exe" -m planetgen.cli.orbits
