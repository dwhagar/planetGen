# scripts/deploy-common.ps1
#
# Steps shared by install.ps1 and update.ps1 on Windows, dot-sourced by
# both so the two can't drift apart. The Windows counterpart of
# scripts/deploy-common.sh (and of install-python-deps.sh and the
# examples/apache/*.sh helpers): a change to a step there belongs here
# too. Each step looks first and changes only what is missing, so running
# it again on a working server does nothing.
#
# The layout is docs/deployment/windows.md's: a venv (default
# C:\srv\planetgen-venv) with the libraries and waitress from
# requirements-server.lock, and the runtime folders (tile cache, jobs,
# logs, NLTK corpus) under C:\ProgramData\planetgen. Runs in Windows
# PowerShell 5.1 and PowerShell 7.
#
# Expects $Root (the checkout), $VenvDir, $DataDir and $ServiceAccount to
# be set by the calling script.

Set-StrictMode -Version 2.0

# Every runtime requirement, as in scripts/install-python-deps.sh's
# REQUIREMENTS (setup.py's install_requires plus the 'api' extra), plus
# waitress. src/tests/test_install_python_deps.py checks these match.
$script:Requirements = @(
    "nltk>=3.9.1",
    "pymysql>=1.1.1",
    "dbutils>=3.1.0",
    "werkzeug>=3.0.0",
    "rich>=13.7.0",
    "flask>=3.0.3",
    "flask-limiter>=3.7.0",
    "redis>=5.0.0",
    "rq>=1.16.0",
    "pyotp>=2.9.0",
    "segno>=1.6.0",
    "markdown>=3.6",
    "cachetools>=5.3.0",
    "sqlalchemy>=2.0.30",
    "alembic>=1.13.0",
    "pydantic>=2.7.0",
    "numpy>=1.26.0",
    "scipy>=1.13.0",
    "astropy>=6.0.0",
    "scikit-image>=0.22.0"
)
$script:ServerRequirement = "waitress>=3.0.1"

function Write-Step([string]$Text) {
    Write-Host ""
    Write-Host "== $Text =="
}

function Assert-Administrator {
    $identity = [Security.Principal.WindowsIdentity]::GetCurrent()
    $principal = New-Object Security.Principal.WindowsPrincipal($identity)
    if (-not $principal.IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)) {
        throw "Run this from an elevated PowerShell (Run as administrator)."
    }
}

# Runs a program (the first argument, with the rest as its arguments) and
# throws if it exits non-zero. A plain function, not an advanced one, so
# arguments such as -m or -c reach the program instead of being taken as
# this function's parameters. Output goes straight to the console
# (stderr is never redirected: Windows PowerShell would turn git's or
# pip's progress lines into errors).
function Invoke-Checked {
    $exe = $args[0]
    $rest = @($args | Select-Object -Skip 1)
    & $exe @rest
    if ($LASTEXITCODE -ne 0) {
        throw "$exe $($rest -join ' ') exited with status $LASTEXITCODE."
    }
}

function Get-VenvPython {
    Join-Path $VenvDir "Scripts\python.exe"
}

# The Python to make the venv with: $env:PYTHON, then the Python launcher's
# newest Python 3, then python.exe on the PATH (not the Microsoft Store
# stub under WindowsApps, which opens the Store instead of running).
function Find-BasePython {
    if ($env:PYTHON) { return $env:PYTHON }
    $launcher = Get-Command py.exe -ErrorAction SilentlyContinue
    if ($launcher) {
        $path = & $launcher.Source -3 -c "import sys; print(sys.executable)"
        if ($LASTEXITCODE -eq 0 -and $path) { return $path.Trim() }
    }
    foreach ($candidate in @(Get-Command python.exe -All -ErrorAction SilentlyContinue)) {
        if ($candidate.Source -notlike "*\WindowsApps\*") { return $candidate.Source }
    }
    throw "No Python found. Install Python 3.12 or later from python.org 'for all users', then run this again."
}

# Probe lines ("<state> <spec> <version> <directory>") for every
# requirement plus waitress, from the venv's python.
function Get-RequirementStates([string]$Python) {
    $specs = $script:Requirements + $script:ServerRequirement
    $lines = & $Python (Join-Path $Root "scripts\probe_requirements.py") @specs
    if ($LASTEXITCODE -ne 0) { return @() }
    @($lines | Where-Object { $_ })
}

function Write-RequirementReport([string]$Python) {
    $failed = $false
    foreach ($line in Get-RequirementStates $Python) {
        $state, $spec, $detail, $where = $line -split " ", 4
        if ($state -eq "ok") {
            Write-Host ("  {0,-10} {1,-24} {2,-8} (pip, in {3})" -f "ok", $spec, $detail, $VenvDir)
        } else {
            Write-Host ("  {0,-10} {1} ({2}: {3})" -f "failed", $spec, $state, $detail)
            $failed = $true
        }
    }
    if ($failed) { throw "Some Python libraries are still unusable (see above)." }
}

# The libraries and waitress in the venv, from requirements-server.lock
# with --require-hashes, so pip installs only the exact files the lock
# names. With -Check (update.ps1) nothing is installed when every
# requirement already imports at or above its floor.
function Install-PythonDeps([switch]$Check) {
    $python = Get-VenvPython
    if ($Check -and (Test-Path $python)) {
        $notOk = @(Get-RequirementStates $python | Where-Object { -not $_.StartsWith("ok ") })
        if ((Get-RequirementStates $python).Count -gt 0 -and $notOk.Count -eq 0) {
            Write-RequirementReport $python
            return
        }
        Write-Host "Missing, too old or not importable in ${VenvDir}: installing from requirements-server.lock."
    }
    if (-not (Test-Path $python)) {
        $base = Find-BasePython
        Invoke-Checked $base -c "import sys; sys.exit(0 if sys.version_info >= (3, 9) else 'Python 3.9 or later is needed')"
        Write-Host "Creating the virtual environment $VenvDir with $base."
        New-Item -ItemType Directory -Force -Path (Split-Path $VenvDir) | Out-Null
        Invoke-Checked $base -m venv $VenvDir
        Invoke-Checked $python -m pip install --quiet --upgrade pip
    }
    Write-Host "Installing the locked libraries and waitress into $VenvDir."
    Invoke-Checked $python -m pip install --quiet --require-hashes -r (Join-Path $Root "requirements-server.lock")
    Write-RequirementReport $python
}

# The NLTK 'words' corpus in a folder every account can read, with
# NLTK_DATA set machine-wide (and in this session) so the services and the
# command-line tools find it. Downloaded only when it's missing.
function Install-NltkWords {
    $dir = Join-Path $DataDir "nltk_data"
    New-Item -ItemType Directory -Force -Path $dir | Out-Null
    $python = Get-VenvPython
    & $python -c @"
import sys
import nltk
try:
    nltk.data.find('corpora/words', paths=[sys.argv[1]])
except LookupError:
    sys.exit(1)
"@ $dir | Out-Null
    if ($LASTEXITCODE -eq 0) {
        Write-Host "NLTK 'words' corpus: present in $dir"
    } else {
        Invoke-Checked $python -c "import nltk, sys; sys.exit(0 if nltk.download('words', download_dir=sys.argv[1]) else 1)" $dir
        Write-Host "NLTK 'words' corpus: installed in $dir"
    }
    if ([Environment]::GetEnvironmentVariable("NLTK_DATA", "Machine") -ne $dir) {
        [Environment]::SetEnvironmentVariable("NLTK_DATA", $dir, "Machine")
        Write-Host "Set NLTK_DATA=$dir for the whole machine."
    }
    $env:NLTK_DATA = $dir
}

# install.ps1 only: config.json from examples/windows/config.json.example
# when there is none, with this install's folders and venv and a new
# secret_key. The database settings are left for the admin. An existing
# config.json is never touched.
function Initialize-Config {
    $config = Join-Path $Root "config.json"
    if (Test-Path $config) {
        Write-Host "config.json: present"
        return
    }
    $bytes = New-Object byte[] 32
    [Security.Cryptography.RandomNumberGenerator]::Create().GetBytes($bytes)
    $secret = -join ($bytes | ForEach-Object { $_.ToString("x2") })
    $text = Get-Content -Raw (Join-Path $Root "examples\windows\config.json.example")
    $text = $text.Replace("C:/srv/planetgen-venv/Scripts/python.exe", ((Get-VenvPython) -replace "\\", "/"))
    $text = $text.Replace("C:/ProgramData/planetgen", ($DataDir -replace "\\", "/"))
    $text = $text.Replace("CHANGE-ME-to-64-random-hex-characters", $secret)
    [IO.File]::WriteAllText($config, $text)
    Write-Host "config.json: written from examples\windows\config.json.example with a new secret_key."
    Write-Host "  Set mysql.password (and the other mysql settings) in $config."
}

# Whether the database settings are still the example's placeholder.
function Test-DatabaseUnconfigured {
    if ($env:PLANETGEN_MYSQL_PASSWORD) { return $false }
    $config = Join-Path $Root "config.json"
    if (-not (Test-Path $config)) { return $false }
    try {
        $mysql = (Get-Content -Raw $config | ConvertFrom-Json).mysql
        return ($null -ne $mysql -and $mysql.password -eq "CHANGE-ME")
    } catch {
        return $false
    }
}

# y/N with a timeout. Returns "" when nobody answers in time, or when
# there is no console to ask on (a scheduled task).
function Read-AnswerWithTimeout([string]$Prompt, [int]$Seconds) {
    try {
        if ([Console]::IsInputRedirected) { return $null }
    } catch {
        return $null
    }
    Write-Host -NoNewline $Prompt
    $answer = ""
    $deadline = (Get-Date).AddSeconds($Seconds)
    while ((Get-Date) -lt $deadline) {
        if (-not [Console]::KeyAvailable) {
            Start-Sleep -Milliseconds 100
            continue
        }
        $key = [Console]::ReadKey($true)
        if ($key.Key -eq "Enter") {
            Write-Host ""
            return $answer
        } elseif ($key.Key -eq "Backspace") {
            if ($answer.Length -gt 0) {
                $answer = $answer.Substring(0, $answer.Length - 1)
                Write-Host -NoNewline "`b `b"
            }
        } elseif ($key.KeyChar -ne [char]0) {
            $answer += $key.KeyChar
            Write-Host -NoNewline $key.KeyChar
        }
    }
    Write-Host ""
    return ""
}

# Brings the configured database up to the current schema. When a
# migration is pending, first asks whether to delete the galaxy data
# instead: y wipes every generated sector and system (src/resetDb.py;
# admin logins are kept) and the empty database is brought to the
# current schema. Anything else, no answer within 30 seconds, or no
# console to ask on keeps the data and migrates it (with migrateDb.py's
# progress bar). Nothing is asked when the database is already current.
function Invoke-MigrateOrReset {
    $python = Get-VenvPython
    Push-Location $Root
    try {
        $status = @(& $python "src\migrateDb.py" --status)
        if ($LASTEXITCODE -ne 0) { throw "src\migrateDb.py --status failed (see above)." }
        $current, $target, $pending, $database = ("$($status[-1])".Trim() -split "\s+")
        if ([int]$pending -gt 0) {
            Write-Host "Database '$database' is at schema v$current; this version needs v$target ($pending migration step(s))."
            $answer = Read-AnswerWithTimeout "Delete all galaxy data in '$database' instead of migrating it? [y/N] (default N in 30s): " 30
            if ($null -eq $answer) {
                Write-Host "(No console to ask on: keeping the data and migrating it.)"
            } elseif ($answer -match "^(y|yes)$") {
                Write-Host "Deleting the galaxy data in '$database' (admin logins are kept)."
                Invoke-Checked $python "src\resetDb.py" --yes
            } else {
                Write-Host "Keeping the data and migrating it."
            }
        }
        Invoke-Checked $python "src\migrateDb.py"
    } finally {
        Pop-Location
    }
}

# Optionally runs the population pass (generate.py population: species,
# civilizations and territories, docs\design\population-and-politics.md).
# Off by default (Boss, 2026-10-01): it runs with -Population, or when
# someone answers y to the prompt (y/N, 30 seconds, default N). With no
# console it is skipped. generate.py population can always be run later.
# The math check (generate.py check-math, TEST.68): known answers from
# real astronomy, identities and sampler distributions. A failure only
# warns -- the update carries on and the site keeps serving -- and
# returns $false so update.ps1 skips the population pass (which would
# refuse anyway) and repeats the warning at the end.
function Test-MathCheck {
    $python = Get-VenvPython
    Push-Location $Root
    try {
        & $python "generate.py" check-math | Out-Host
        $ok = ($LASTEXITCODE -eq 0)
    } finally {
        Pop-Location
    }
    if (-not $ok) {
        Write-Warning "The math check failed (above). Bulk generation refuses to start until it passes; run generate.py check-math -v for every check."
    }
    return $ok
}

function Invoke-OptionalPopulation([switch]$Run) {
    if ($Run) {
        $answer = "y"
    } else {
        $answer = Read-AnswerWithTimeout "Run the population pass now (species, civilizations, territories)? [y/N] (default N in 30s): " 30
    }
    if ($answer -match "^(y|yes)$") {
        Write-Host "Running the population pass."
        $python = Get-VenvPython
        Push-Location $Root
        try {
            Invoke-Checked $python "generate.py" population
        } finally {
            Pop-Location
        }
    } else {
        Write-Host "Skipping the population pass (run generate.py population any time, or pass -Population)."
    }
}

# The two logs (OPS.5), worked out the way the program does
# (stellarObjects\appconfig.py): the debug log is PLANETGEN_LOG_FILE, else
# "log_file" in config.json, else its default; the always-on activity log
# is planetgen.log in PLANETGEN_LOG_DIR, else "log_dir", else logs under
# the checkout. One object per log: Name, File, Dir and Setting (what to
# change to move it). Empty when config.json can't be read (a warning
# says so).
$script:LogLocations = $null

function Get-LogLocations {
    if ($null -ne $script:LogLocations) { return $script:LogLocations }
    $python = Get-VenvPython
    $lines = @(& $python -I -c @"
import importlib.util, os, sys
spec = importlib.util.spec_from_file_location('appconfig', os.path.join(sys.argv[1], 'src', 'planetgen', 'util', 'appconfig.py'))
appconfig = importlib.util.module_from_spec(spec)
spec.loader.exec_module(appconfig)
config = appconfig.load_config()
print(os.path.abspath(appconfig.log_file_path(config)))
print(appconfig.activity_log_path(config))
"@ $Root)
    if ($LASTEXITCODE -ne 0 -or $lines.Count -lt 2) {
        $script:LogLocations = @()
        Write-Warning "Couldn't work out where the logs go (see above), so their folders were not checked. Check config.json's `"log_file`" and `"log_dir`" (each a path, or leave them out for the defaults), then run update.ps1 again."
        return @()
    }
    $debugFile = "$($lines[0])".Trim()
    $activityFile = "$($lines[1])".Trim()
    $script:LogLocations = @(
        [PSCustomObject]@{ Name = "debug log"; File = $debugFile; Dir = (Split-Path $debugFile); Setting = "`"log_file`" in config.json (or PLANETGEN_LOG_FILE)" },
        [PSCustomObject]@{ Name = "activity log"; File = $activityFile; Dir = (Split-Path $activityFile); Setting = "`"log_dir`" in config.json (or PLANETGEN_LOG_DIR)" }
    )
    $script:LogLocations
}

# The tile cache and jobs folders config.json names (read the way
# create-cache-dir.sh reads them, with examples\apache\deploy-paths.py),
# the debug log's folder, the always-on activity log's folder (logs under
# the checkout unless "log_dir" says otherwise; the app rotates it
# itself), and the logs folder the services write to.
function Get-RuntimeDirs {
    $python = Get-VenvPython
    $paths = @(& $python -I (Join-Path $Root "examples\apache\deploy-paths.py") $Root)
    if ($LASTEXITCODE -ne 0) { throw "Could not read the folders from config.json (see above)." }
    $dirs = @()
    $logDirs = @(Get-LogLocations | ForEach-Object { $_.Dir })
    foreach ($path in @($paths[0], $paths[1]) + $logDirs + @(Join-Path $DataDir "logs")) {
        if ($path -and ($dirs -notcontains $path)) { $dirs += $path }
    }
    $dirs
}

# Whether a folder is one of the logs' (a problem there only warns).
function Test-IsLogDir([string]$Dir) {
    $logDirs = @(Get-LogLocations | ForEach-Object { $_.Dir })
    return (($logDirs -contains $Dir) -or ($Dir -eq (Join-Path $DataDir "logs")))
}

# The warning for a log that can't be set up, with the exact commands that
# fix it (OPS.5). Never throws: a log never stops an install or update.
function Write-LogWarning($Log, [string]$Problem) {
    Write-Warning "The $($Log.Name) $($Log.File) isn't fully set up: $Problem"
    Write-Host "  The install/update carries on; the site still runs, but this log may not be written. Fix it"
    Write-Host "  from an elevated PowerShell with:"
    Write-Host "    New-Item -ItemType Directory -Force -Path `"$($Log.Dir)`""
    Write-Host "    icacls `"$($Log.Dir)`" /grant `"${ServiceAccount}:(OI)(CI)M`""
    Write-Host "  or point $($Log.Setting) at a folder this server can write, then run update.ps1 again."
}

function New-RuntimeDirs {
    foreach ($dir in Get-RuntimeDirs) {
        if (Test-IsLogDir $dir) {
            $reason = "it isn't there afterwards"
            try {
                New-Item -ItemType Directory -Force -Path $dir -ErrorAction Stop | Out-Null
            } catch {
                $reason = $_.Exception.Message
            }
            if (Test-Path -LiteralPath $dir -PathType Container) {
                Write-Host "Runtime folder: $dir"
            } else {
                Write-Warning "Couldn't create the log folder ${dir}: $reason (the next step says how to fix it)."
            }
            continue
        }
        New-Item -ItemType Directory -Force -Path $dir | Out-Null
        Write-Host "Runtime folder: $dir"
    }
}

# Whether $Account (or a group every account is in) may change files in
# $Dir, from the folder's access rules: an Allow rule with write rights
# and no Deny rule for them.
function Test-AccountCanWrite([string]$Dir, [string]$Account) {
    try {
        $sid = (New-Object Security.Principal.NTAccount($Account)).Translate([Security.Principal.SecurityIdentifier]).Value
        $rules = (Get-Acl -Path $Dir -ErrorAction Stop).Access
    } catch {
        return $false
    }
    # The account itself, Everyone, Authenticated Users, BUILTIN\Users.
    $sids = @($sid, "S-1-1-0", "S-1-5-11", "S-1-5-32-545")
    $write = [Security.AccessControl.FileSystemRights]::WriteData -bor [Security.AccessControl.FileSystemRights]::AppendData
    $allowed = $false
    foreach ($rule in $rules) {
        try {
            $ruleSid = $rule.IdentityReference.Translate([Security.Principal.SecurityIdentifier]).Value
        } catch {
            continue
        }
        if ($sids -notcontains $ruleSid) { continue }
        if (($rule.FileSystemRights -band $write) -eq 0) { continue }
        if ($rule.AccessControlType -eq "Deny") { return $false }
        $allowed = $true
    }
    return $allowed
}

# OPS.5: both logs' folders exist and the app's account can write in them.
# Only ever warns, with the commands that fix it for this machine.
function Test-LogLocations {
    $accountExists = Test-AccountExists $ServiceAccount
    foreach ($log in Get-LogLocations) {
        if (-not (Test-Path -LiteralPath $log.Dir -PathType Container)) {
            Write-LogWarning $log "its folder $($log.Dir) doesn't exist and couldn't be created (a missing drive, or no rights)."
        } elseif (-not $accountExists) {
            Write-LogWarning $log "the account '$ServiceAccount' doesn't exist yet, so it couldn't be given the folder."
        } elseif (-not (Test-AccountCanWrite $log.Dir $ServiceAccount)) {
            Write-LogWarning $log "'$ServiceAccount' can't write in $($log.Dir)."
        } else {
            $label = $log.Name.Substring(0, 1).ToUpper() + $log.Name.Substring(1)
            Write-Host "${label}: $($log.File) ($ServiceAccount can write in $($log.Dir))"
        }
    }
}

function Test-AccountExists([string]$Account) {
    try {
        (New-Object Security.Principal.NTAccount($Account)).Translate([Security.Principal.SecurityIdentifier]) | Out-Null
        return $true
    } catch {
        return $false
    }
}

# The same rules set-permissions.sh applies on Linux, with icacls: the
# account the app runs as can read the checkout, the venv and the corpus
# but not change them, can change only the runtime folders, and
# config.json (the database password and secret_key) is readable only by
# Administrators, SYSTEM and that account. Skipped with a warning when
# the account doesn't exist yet (NT SERVICE\planetgen appears once the
# service does); run update.ps1 again afterwards.
function Set-PlanetGenPermissions {
    if (-not (Test-AccountExists $ServiceAccount)) {
        Write-Warning "The account '$ServiceAccount' doesn't exist yet, so permissions were not set. Create the service (docs\deployment\windows.md), then run update.ps1 (or pass -ServiceAccount 'IIS AppPool\planetgen' for IIS)."
        Test-LogLocations
        return
    }
    $read = @($Root, $VenvDir, (Join-Path $DataDir "nltk_data"))
    foreach ($dir in $read) {
        Invoke-Checked icacls.exe $dir /grant "${ServiceAccount}:(OI)(CI)RX" /Q
    }
    foreach ($dir in Get-RuntimeDirs) {
        if (Test-IsLogDir $dir) {
            # A log folder that can't be granted only warns (OPS.5);
            # Test-LogLocations below says how to fix it.
            if (Test-Path -LiteralPath $dir -PathType Container) {
                & icacls.exe $dir /grant "${ServiceAccount}:(OI)(CI)M" /Q
            }
            continue
        }
        Invoke-Checked icacls.exe $dir /grant "${ServiceAccount}:(OI)(CI)M" /Q
    }
    $config = Join-Path $Root "config.json"
    if (Test-Path $config) {
        Invoke-Checked icacls.exe $config /inheritance:r /grant:r "*S-1-5-32-544:F" "*S-1-5-18:F" "${ServiceAccount}:R" /Q
    }
    Write-Host "Permissions: $ServiceAccount reads the code, venv and corpus, changes only the runtime folders; config.json is Administrators, SYSTEM and $ServiceAccount only."
    Test-LogLocations
}

# Imports the web app and the generator package with the venv's Python,
# so a library that doesn't import (or the checkout's own code failing)
# shows up here rather than as a 500. Runs as the current administrator:
# Windows can't start a process as a virtual service account from here.
function Test-AppImports {
    $python = Get-VenvPython
    & $python -c @"
import os, sys
root = sys.argv[1]
sys.path[:0] = [os.path.join(root, 'src', 'html'), os.path.join(root, 'src')]
import stellarObjects
from api.app import create_app
import web
"@ $Root
    if ($LASTEXITCODE -ne 0) { throw "The web app does not import with $python (see above)." }
    Write-Host "The web app and stellarObjects import cleanly with $python."
}

# The Redis server at config.json's redis.url (OPS.21), which the work
# queue and the rate limits will use. Redis has no supported native
# Windows build, so nothing is installed here: it only checks that one
# answers, and otherwise says to run Memurai (a Redis-compatible Windows
# service) or Redis in WSL (docs/deployment/windows.md). Nothing uses
# Redis yet, so a server that doesn't answer only warns.
function Test-Redis {
    $python = Get-VenvPython
    $url = & $python -c "import sys; sys.path.insert(0, sys.argv[1]); from planetgen.util.appconfig import load_config; print(load_config()['redis']['url'])" (Join-Path $Root "src")
    if ($LASTEXITCODE -ne 0 -or -not $url) {
        Write-Warning "Couldn't read redis.url from config.json; skipping the Redis check."
        return
    }
    # The failure is caught in Python, not redirected here: Windows
    # PowerShell turns a native program's stderr into errors.
    & $python -c "import sys, redis`ntry: redis.Redis.from_url(sys.argv[1], socket_connect_timeout=3).ping()`nexcept Exception: sys.exit(1)" $url
    if ($LASTEXITCODE -eq 0) {
        Write-Host "Redis: answering at $url"
        return
    }
    Write-Warning ("No Redis server answers at $url (config.json's redis.url). Redis has no native " +
        "Windows build: install Memurai (https://www.memurai.com/) or run Redis in WSL " +
        "(docs/deployment/windows.md). Nothing needs it yet; the work queue will.")
    # A warning, not a failure: don't let the probe's exit code become the
    # script's (powershell -command exits with the last $LASTEXITCODE).
    $global:LASTEXITCODE = 0
}

# What restarts the app after an update, for the closing message.
function Get-RestartHint {
    if (Get-Service -Name planetgen -ErrorAction SilentlyContinue) { return "Restart-Service planetgen" }
    if (Get-Command Restart-WebAppPool -ErrorAction SilentlyContinue) { return "Restart-WebAppPool planetgen   (IIS)" }
    return "restart the planetgen service or IIS app pool"
}
