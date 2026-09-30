<#
.SYNOPSIS
    Schedules planetGen's monthly maintenance with Task Scheduler on
    native Windows.

.DESCRIPTION
    The Windows counterpart of install-maintenance-timer.sh (systemd
    timers on Linux, launchd on macOS), on the same schedule:

      - "planetGen update" runs update.ps1 at 03:00 on the 1st of each
        month (unless -SkipUpdateTask). Whatever is on the tracked branch
        gets deployed with no human review in between; leave it out if
        only a person should update this site. With no console to ask
        on, a pending migration keeps the data and migrates it.
      - "planetGen orbits (<database>)" runs src\updateOrbits.py against
        each database at 03:30 on the 1st.

    Both run as SYSTEM, and run as soon as possible after a start missed
    while the machine was off (like systemd's Persistent=true). Each task
    runs a small .cmd written to <DataDir>\tasks with this install's
    paths, which keeps the task's command line short. Credentials come
    from config.json. Re-running replaces the tasks.

    Run from an elevated PowerShell:

        powershell -ExecutionPolicy Bypass -File .\examples\maintenance\install-maintenance-task.ps1 [-Database planetgen,other] [-SkipUpdateTask]

.PARAMETER Database
    The databases to schedule an orbit update for. Default:
    $env:PLANETGEN_MYSQL_DATABASE, else planetgen.

.PARAMETER SkipUpdateTask
    Only schedule the orbit updates, not update.ps1.
#>
[CmdletBinding()]
param(
    [string[]]$Database = @($(if ($env:PLANETGEN_MYSQL_DATABASE) { $env:PLANETGEN_MYSQL_DATABASE } else { "planetgen" })),
    [switch]$SkipUpdateTask,
    [string]$VenvDir = "C:\srv\planetgen-venv",
    [string]$DataDir = (Join-Path $env:ProgramData "planetgen"),
    [string]$ServiceAccount = "NT SERVICE\planetgen"
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
. (Join-Path $Root "scripts\deploy-common.ps1")

Assert-Administrator

$tasksDir = Join-Path $DataDir "tasks"
New-Item -ItemType Directory -Force -Path $tasksDir | Out-Null
New-Item -ItemType Directory -Force -Path (Join-Path $DataDir "logs") | Out-Null
$python = Get-VenvPython
$nltk = Join-Path $DataDir "nltk_data"
$logs = Join-Path $DataDir "logs"

# A monthly task on the 1st at $Time running $Cmd as SYSTEM, replacing
# any task of the same name. schtasks.exe because New-ScheduledTaskTrigger
# has no monthly trigger.
# The .cmd path goes to schtasks unquoted (Windows PowerShell mangles
# quotes inside a native command's argument), so it must have no spaces.
function Register-MonthlyTask([string]$Name, [string]$Time, [string]$Cmd) {
    if ($Cmd -match "\s") { throw "-DataDir must not contain spaces (the task runs $Cmd)." }
    Invoke-Checked schtasks.exe /Create /F /TN $Name /SC MONTHLY /D 1 /ST $Time /RU SYSTEM /TR $Cmd
    $settings = New-ScheduledTaskSettingsSet -StartWhenAvailable -AllowStartIfOnBatteries `
        -DontStopIfGoingOnBatteries -ExecutionTimeLimit (New-TimeSpan -Hours 12)
    Set-ScheduledTask -TaskName $Name -Settings $settings | Out-Null
    Write-Host "  $Name at $Time on the 1st of each month ($Cmd)"
}

Write-Step "Scheduling planetGen's maintenance tasks"

if (-not $SkipUpdateTask) {
    $cmd = Join-Path $tasksDir "update.cmd"
    Set-Content -Encoding ASCII -Path $cmd -Value @(
        "@echo off",
        "rem Written by examples\maintenance\install-maintenance-task.ps1.",
        "powershell.exe -NoProfile -ExecutionPolicy Bypass -File `"$Root\update.ps1`" -VenvDir `"$VenvDir`" -DataDir `"$DataDir`" -ServiceAccount `"$ServiceAccount`" >> `"$logs\update.log`" 2>&1"
    )
    Register-MonthlyTask "planetGen update" "03:00" $cmd
} else {
    Write-Host "  -SkipUpdateTask given: not scheduling update.ps1"
}

foreach ($db in $Database) {
    if ($db -notmatch '^[A-Za-z0-9_]+$') { throw "Not a database name: '$db'" }
    $cmd = Join-Path $tasksDir "orbits-$db.cmd"
    Set-Content -Encoding ASCII -Path $cmd -Value @(
        "@echo off",
        "rem Written by examples\maintenance\install-maintenance-task.ps1.",
        "set `"PLANETGEN_MYSQL_DATABASE=$db`"",
        "set `"NLTK_DATA=$nltk`"",
        "cd /d `"$Root`"",
        "`"$python`" src\updateOrbits.py >> `"$logs\orbits-$db.log`" 2>&1"
    )
    Register-MonthlyTask "planetGen orbits ($db)" "03:30" $cmd
}

Write-Host ""
Write-Host "Done. Check them in Task Scheduler, or with: schtasks /Query /TN `"planetGen orbits ($($Database[0]))`""
Write-Host "Run one now with: schtasks /Run /TN `"planetGen orbits ($($Database[0]))`""
Write-Host "Logs: $logs"
