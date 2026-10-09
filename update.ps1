<#
.SYNOPSIS
    Pulls the latest planetGen and checks everything the site needs, on
    native Windows, without reinstalling what's already there.

.DESCRIPTION
    update.ps1 is update.sh's Windows counterpart and follows the same
    steps (a change to one belongs in both):

      1. Pulls: git fetch, then git reset --hard to origin's branch tip,
         so the checkout matches the branch exactly. Uncommitted changes
         to tracked files are listed and overwritten; untracked files,
         config.json among them, are never touched.
      2. Checks every Python library in the venv (and waitress) by
         importing it; only when one is missing, too old or broken does
         pip install from requirements-server.lock (--require-hashes).
      3. The NLTK 'words' corpus: fetched only if it's missing.
      4. planetgen.cli.generate check-math, the math check (TEST.68). A failure
         only warns (repeated at the end).
      5. planetgen.cli.migrate, a no-op when the database is current. When a
         migration is pending it first asks (y/N, 30 seconds, default N)
         whether to delete the galaxy data instead; a scheduled run with
         no console keeps the data and migrates it. The update never runs
         the population pass (OPS.7); run planetgen.cli.generate population by hand
         when wanted.
      6. The tile cache, jobs and log folders. A log folder that can't be
         made only warns.
      7. Permissions for the app's account (icacls), in case new folders
         came in, then a check that it can write both logs' folders (a
         warning with the fixing commands if not; never a stop).
      8. Imports the web app, so anything unusable fails here instead of
         as a 500.

    Run from an elevated PowerShell in the checkout:

        powershell -ExecutionPolicy Bypass -File .\update.ps1

    then restart the app (the closing message says how). Takes the same
    -VenvDir, -DataDir and -ServiceAccount as install.ps1.
#>
[CmdletBinding()]
param(
    [string]$VenvDir = "C:\srv\planetgen-venv",
    [string]$DataDir = (Join-Path $env:ProgramData "planetgen"),
    [string]$ServiceAccount = "NT SERVICE\planetgen"
)

$ErrorActionPreference = "Stop"
$Root = $PSScriptRoot
. (Join-Path $Root "scripts\deploy-common.ps1")

Assert-Administrator
Set-Location $Root

$inside = & git rev-parse --is-inside-work-tree
if ($LASTEXITCODE -ne 0 -or "$inside".Trim() -ne "true") {
    throw "$Root is not a git checkout; can't pull an update here."
}

Write-Step "1/8: Pulling the latest changes"
$dirty = @(& git status --porcelain)
if ($dirty.Count -gt 0) {
    Write-Warning "Uncommitted local changes in $Root will be overwritten:"
    $dirty | ForEach-Object { Write-Host "  $_" }
}
$branch = "$(& git rev-parse --abbrev-ref HEAD)".Trim()
if ($branch -eq "HEAD") {
    throw "The repository is in a detached HEAD state; check out a branch first."
}
$before = "$(& git rev-parse HEAD)".Trim()
Invoke-Checked git fetch origin $branch
# reset --hard rather than pull: tracked files always match origin's tip,
# whatever was committed or edited here. Never git clean: config.json is
# untracked and must survive.
Invoke-Checked git reset --hard "origin/$branch"
$after = "$(& git rev-parse HEAD)".Trim()
if ($before -eq $after) {
    Write-Host "Already up to date ($before)."
} else {
    Write-Host "Updated $before..$after`:"
    & git log --oneline "$before..$after"
}
# Again, now that the pull may have changed it.
. (Join-Path $Root "scripts\deploy-common.ps1")

Write-Step "2/8: Checking the Python libraries in $VenvDir"
Install-PythonDeps -Check

Write-Step "3/8: Checking the NLTK 'words' corpus"
Install-NltkWords

Write-Step "4/8: Checking the generator's math"
$mathOk = Test-MathCheck

Write-Step "5/8: Migrating the configured MySQL database to the current schema"
if (Test-DatabaseUnconfigured) {
    Write-Host "Skipped: config.json still has the example's database password. Set the mysql settings in"
    Write-Host "  $(Join-Path $Root 'config.json'), then run this again."
} else {
    Invoke-MigrateOrReset
    Invoke-RecordVersionKey
}

Write-Step "6/8: Checking the tile cache, jobs and log folders"
New-RuntimeDirs

Write-Step "7/8: Setting permissions for $ServiceAccount"
Set-PlanetGenPermissions

Write-Step "8/8: Checking that the web app imports"
Test-AppImports
Test-Redis
Start-GalaxyMapWarmup

Write-Host ""
if (-not $mathOk) {
    Write-Warning "The math check failed (step 4): bulk generation refuses to start until it passes."
}
Write-Host "The update doesn't run the population pass (species, civilizations, territories); when wanted:"
Write-Host "  $(Get-VenvPython) -m planetgen.cli.generate population"
if ($before -ne $after) {
    Write-Host "Done. Restart the app so the site runs the new code:"
    Write-Host "  $(Get-RestartHint)"
} else {
    Write-Host "Done. Nothing new was pulled."
}
