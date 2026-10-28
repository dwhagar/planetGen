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
#
# tests/test_log_shell.py runs these against `logrecord.display_line`, so the two cannot drift apart.

log_record() {
    local severity="$1" message="$2" stamp upper line
    case "$severity" in
        debug|info|notice|warning|error|critical) ;;
        *) severity=info ;;
    esac
    upper="$(printf '%s' "$severity" | tr '[:lower:]' '[:upper:]')"
    stamp="$(date '+%Y-%m-%d %H:%M:%S %Z')"
    while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ "$severity" == warning || "$severity" == error || "$severity" == critical ]]; then
            printf '%s %-8s %s: %s\n' "$stamp" "$upper" "${LOG_COMPONENT:-update}" "$line" >&2
        else
            printf '%s %-8s %s: %s\n' "$stamp" "$upper" "${LOG_COMPONENT:-update}" "$line"
        fi
    done <<< "$message"
}

log_info() { log_record info "$1"; }
log_notice() { log_record notice "$1"; }
log_warning() { log_record warning "$1"; }
log_error() { log_record error "$1"; }

# The divider of one stage, as `logrecord.stage_start` writes it: "=== Stage 4 of 9: label ===", an info record.
log_stage() {
    log_record info "=== Stage $1 of $2: $3 ==="
}
