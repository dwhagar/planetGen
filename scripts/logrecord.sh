# scripts/logrecord.sh
#
# The shell's side of the logging standard (UX.123, docs/design/operations/logging-standard.md): the display line
# that `planetgen/util/logrecord.py` draws for a record, so update.sh's own output reads like the Generate page's and the
# command line's debug log:
#
#     2026-10-24 14:40:00 EDT INFO     update: Already up to date (abc123).
#     2026-10-24 14:40:00 EDT INFO     update: === Stage 4 of 9: Checking the generator's math ===
#
# local time with its zone, the severity upper case and padded to eight, the component and the message. Severities are
# the RFC 5424 names the standard uses (debug, info, notice, warning, error, critical); warning and above go to
# stderr. A message of several lines is several records, one to a line. The component is $LOG_COMPONENT (default
# `update`). Written for macOS's bash 3.2 as well: no associative arrays.
#
#     log_record SEVERITY MESSAGE     one record
#     log_info, log_notice, log_warning, log_error MESSAGE
#     log_stage FIRST TOTAL LABEL     the divider of a stage
#     log_lines [SEVERITY]            every non-empty line read from stdin, as a record (a heredoc of instructions)
#     log_tail PREFIX                 the lines a command wrote to stdout, as records (OPS.98)
#     log_summary_note SEVERITY MSG   a line for the closing summary, and what log_summary repeats (OPS.100)
#     log_summary HEADLINE            the closing block of a script
#
# The severities, the stage word and the divider are generated from the module (the block below, scripts/build_log_vocabulary.py), and
# tests/test_log_shell.py runs these against `logrecord.display_line`, so the two cannot drift apart.

# BEGIN generated vocabulary (scripts/build_log_vocabulary.py from planetgen/util/logrecord.py; do not edit)
LOG_SEVERITIES="debug info notice warning error critical"
LOG_STDERR_SEVERITIES="warning error critical"
LOG_STAGE_WORD="Stage"
LOG_DIVIDER="==="
LOG_COUNTING_WORDS="Stage Part Phase"
LOG_COMPONENTS="update install warm tiles orbits generate queue web job activity command lib logs permissions cache-dir swap python-deps"
# END generated vocabulary

log_record() {
    local severity="$1" message="$2" stamp upper line record
    [[ " $LOG_SEVERITIES " == *" $severity "* ]] || severity=info
    upper="$(printf '%s' "$severity" | tr '[:lower:]' '[:upper:]')"
    stamp="$(date '+%Y-%m-%d %H:%M:%S %Z')"
    while IFS= read -r line || [[ -n "$line" ]]; do
        # A warning or an error is kept for the closing summary (OPS.100): its first line only, the lines under it are its fix.
        if [[ -n "${LOG_WARNING_FILE:-}" && -w "${LOG_WARNING_FILE:-}" && " $LOG_STDERR_SEVERITIES " == *" $severity "* && "$line" != [[:space:]]* ]]; then
            printf '%s\n' "$line" >> "$LOG_WARNING_FILE"
        fi
        record="$(printf '%s %-8s %s: %s' "$stamp" "$upper" "${LOG_COMPONENT:-update}" "$line")"
        if [[ " $LOG_STDERR_SEVERITIES " == *" $severity "* ]]; then
            printf '%s\n' "$record" >&2
        else
            printf '%s\n' "$record"
        fi
        # The same record goes to the script's log file when it has one (OPS.101).
        if [[ -n "${LOG_TEE:-}" ]]; then
            printf '%s\n' "$record" >> "$LOG_TEE" 2>/dev/null || true
        fi
    done <<< "$message"
}

log_info() { log_record info "$1"; }
log_notice() { log_record notice "$1"; }
log_warning() { log_record warning "$1"; }
log_error() { log_record error "$1"; }

# The divider of one stage, as `logrecord.stage_start` writes it: "=== Stage 4 of 9: label ===", an info record.
log_stage() {
    log_record info "$LOG_DIVIDER $LOG_STAGE_WORD $1 of $2: $3 $LOG_DIVIDER"
    UPDATE_STEP="$(log_lower "$LOG_STAGE_WORD") $1 of $2, $(log_lower "${3:0:1}")${3:1}"
}

# Where a script is, in one of the three counting words (OPS.103): `log_count Part 1 2` is "Part 1 of 2". The only place a script words a count;
# src/tests/test_no_own_log_formatting.py fails on "stage 1 of 2" typed out by hand anywhere else.
log_count() {
    local word="$1"
    [[ " $LOG_COUNTING_WORDS " == *" $word "* ]] || word="$LOG_STAGE_WORD"
    printf '%s %s of %s' "$word" "$2" "$3"
}

# The text in lower case, for the middle of a sentence ("part 1 of 2, ...").
log_lower() { printf '%s' "$1" | tr '[:upper:]' '[:lower:]'; }

# Every non-empty line read from stdin, as a record of the given severity (default info): the instructions a script prints in a heredoc.
log_lines() {
    local severity="${1:-info}" line
    while IFS= read -r line || [[ -n "$line" ]]; do
        [[ -z "$line" ]] && continue
        log_record "$severity" "$line"
    done
}

# The closing summary of a long script (OPS.100). `log_summary_begin` starts the clock and a file that keeps every warning and error the
# script prints (a file, so a warning printed inside `$(...)` or a pipe is kept too); `log_summary TITLE DETAILS` ends with one line,
# "=== Update done after 20:53: DETAILS ===", then every warning repeated, so the last thing on the screen is what went wrong.
log_summary_begin() {
    LOG_STARTED=$SECONDS
    LOG_WARNING_TMP="$(mktemp "${TMPDIR:-/tmp}/planetgen-warnings.XXXXXX")"
    LOG_WARNING_FILE="$LOG_WARNING_TMP"
    export LOG_WARNING_FILE
}

# Forgets the warnings file; the script's EXIT trap calls it so an abort leaves nothing in the temporary folder.
log_summary_forget() {
    [[ -n "${LOG_WARNING_TMP:-}" ]] && rm -f "$LOG_WARNING_TMP"
    [[ -n "${LOG_TEE_BUFFER:-}" ]] && rm -f "$LOG_TEE_BUFFER"
    LOG_WARNING_FILE=""
    LOG_WARNING_TMP=""
}

# A length of time the way the summary shows it: 4:07, or 1:04:07 from an hour on.
log_duration() {
    local total="$1"
    if (( total >= 3600 )); then
        printf '%d:%02d:%02d' $(( total / 3600 )) $(( total % 3600 / 60 )) $(( total % 60 ))
    else
        printf '%d:%02d' $(( total / 60 )) $(( total % 60 ))
    fi
}

log_summary() {
    local title="$1" details="${2:-}" warnings="" count=0 line
    if [[ -n "${LOG_WARNING_FILE:-}" && -r "$LOG_WARNING_FILE" ]]; then
        warnings="$(cat "$LOG_WARNING_FILE")"
    fi
    LOG_WARNING_FILE=""      # the repeated warnings are not new warnings
    log_summary_forget
    log_info "$LOG_DIVIDER $title done after $(log_duration $(( SECONDS - ${LOG_STARTED:-0} )))${details:+: $details} $LOG_DIVIDER"
    if [[ -z "$warnings" ]]; then
        log_info "No warnings during this $(printf '%s' "$title" | tr '[:upper:]' '[:lower:]')."
        return 0
    fi
    while IFS= read -r line; do
        [[ -n "$line" ]] && count=$(( count + 1 ))
    done <<< "$warnings"
    log_warning "$count warning$([[ $count == 1 ]] || printf s) during this $(printf '%s' "$title" | tr '[:upper:]' '[:lower:]'), repeated here:"
    while IFS= read -r line; do
        [[ -n "$line" ]] && log_warning "  $line"
    done <<< "$warnings"
}

# The script's own log file (OPS.101): `log_file_begin` starts keeping every record in a temporary file, because the log folder is not known
# until the Python environment is (the script may have just pulled the code that says where it is); `log_file_to PATH` appends what was kept
# to PATH and from then on records are appended there directly. A path that cannot be written leaves the screen as it was.
log_file_begin() {
    LOG_TEE_BUFFER="$(mktemp "${TMPDIR:-/tmp}/planetgen-log.XXXXXX")"
    LOG_TEE="$LOG_TEE_BUFFER"
}

log_file_to() {
    local path="$1"
    [[ -n "${LOG_TEE_BUFFER:-}" && -n "$path" ]] || return 0
    mkdir -p "$(dirname "$path")" 2>/dev/null || return 0
    if cat "$LOG_TEE_BUFFER" >> "$path" 2>/dev/null; then
        LOG_TEE="$path"
        rm -f "$LOG_TEE_BUFFER"
        LOG_TEE_BUFFER=""
    fi
}
