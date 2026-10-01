<#
.SYNOPSIS
    Installs planetGen's web interface on native Windows.

.DESCRIPTION
    install.ps1 is install.sh's Windows counterpart and follows the same
    steps (install.sh on Linux and macOS; a change to one belongs in all
    three), with the layout of docs/deployment/windows.md:

      1. A virtual environment (-VenvDir, default C:\srv\planetgen-venv)
         made with the Python launcher's newest Python 3, holding the
         libraries and waitress from requirements-server.lock, installed
         with --require-hashes.
      2. The NLTK 'words' corpus in <DataDir>\nltk_data, with NLTK_DATA
         set machine-wide (before the database step, which imports
         stellarObjects and would otherwise fetch it into the admin's
         own profile).
      3. config.json from examples\windows\config.json.example when there
         is none (this install's folders, a new secret_key; you fill in
         the database settings), then src\migrateDb.py with the same
         migrate-or-delete question as install.sh (y/N, 30 seconds,
         default N) and its progress bar. Skipped, with a note, while
         mysql.password is still the example's CHANGE-ME. Then it offers
         (y/N, 30 seconds, default N; skipped with no console) to run the
         population pass, generate.py population.
      4. The tile cache, jobs and log folders (<DataDir>\tiles, jobs, logs
         by default, or wherever config.json points).
      5. Permissions with icacls, as set-permissions.sh does on Linux: the
         app's account (-ServiceAccount) reads the code, venv and corpus,
         changes only the runtime folders, and config.json is readable by
         Administrators, SYSTEM and that account only. Then it checks
         that account can write both logs' folders (the debug log's and
         the activity log's, default or configured); one it can't only
         warns, with the New-Item and icacls commands that fix it.
      6. Imports the web app with the venv's Python, so a problem shows
         here instead of as a 500.

    The web server in front (IIS, Caddy or Apache Lounge, each with
    waitress) is never set up automatically; the closing message points
    at the guide. Run from an elevated PowerShell in the checkout:

        powershell -ExecutionPolicy Bypass -File .\install.ps1

.PARAMETER VenvDir
    Where the virtual environment goes. Default C:\srv\planetgen-venv.

.PARAMETER DataDir
    Where the runtime folders and the NLTK corpus go. Default
    C:\ProgramData\planetgen.

.PARAMETER ServiceAccount
    The account the app runs as: NT SERVICE\planetgen for the waitress
    service (the default), or IIS AppPool\planetgen for IIS.

.PARAMETER SkipDatabase
    Leaves out the database step, for a host whose database isn't set up
    yet (and for CI). Run update.ps1 once it is.

.PARAMETER Population
    Runs the population pass (generate.py population: species,
    civilizations, territories) after the database step without asking.
    Off by default.
#>
[CmdletBinding()]
param(
    [string]$VenvDir = "C:\srv\planetgen-venv",
    [string]$DataDir = (Join-Path $env:ProgramData "planetgen"),
    [string]$ServiceAccount = "NT SERVICE\planetgen",
    [switch]$SkipDatabase,
    [switch]$Population
)

$ErrorActionPreference = "Stop"
$Root = $PSScriptRoot
. (Join-Path $Root "scripts\deploy-common.ps1")

Assert-Administrator

Write-Step "1/6: Installing the Python libraries and waitress into $VenvDir"
Install-PythonDeps

Write-Step "2/6: Fetching the NLTK 'words' corpus"
Install-NltkWords

Write-Step "3/6: Migrating the configured MySQL database to the current schema"
Initialize-Config
if ($SkipDatabase) {
    Write-Host "Skipped (-SkipDatabase). Run update.ps1 once the database is set up."
} elseif (Test-DatabaseUnconfigured) {
    Write-Host "Skipped: config.json still has the example's database password. Set the mysql settings in"
    Write-Host "  $(Join-Path $Root 'config.json'), then run update.ps1."
} else {
    Invoke-MigrateOrReset
    Invoke-OptionalPopulation -Run:$Population
}

Write-Step "4/6: Creating the tile cache, jobs and log folders"
New-RuntimeDirs

Write-Step "5/6: Setting permissions for $ServiceAccount"
Set-PlanetGenPermissions

Write-Step "6/6: Checking that the web app imports"
Test-AppImports

$waitress = Join-Path $VenvDir "Scripts\waitress-serve.exe"
Write-Host @"

------------------------------------------------------------------------
Install steps complete. What remains is running waitress and the web
server in front of it, never set up automatically. Try waitress by hand:

    cd $(Join-Path $Root 'src\html')
    $waitress --listen=127.0.0.1:8000 --threads=5 --no-clear-untrusted-proxy-headers wsgi:application

then pick an option in docs\deployment\windows.md (IIS, Caddy or Apache
Lounge), and log in at https://<your site>/login with the admin username
and password printed once in step 3/6 above, and change both.
------------------------------------------------------------------------
"@
