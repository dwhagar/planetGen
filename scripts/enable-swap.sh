#!/usr/bin/env bash
# scripts/enable-swap.sh
#
# Turns on a swap file and keeps it across reboots (PERF.181). A server with no swap has nothing between a memory
# spike and the kernel's out-of-memory killer, which ends MySQL, Redis or the web server instead of slowing down.
# Safe to run again: it does nothing when swap is already on, and says what it did.
#
# Usage (as root):
#   sudo scripts/enable-swap.sh            # a 8G file at /swapfile
#   sudo scripts/enable-swap.sh 16G        # another size (a number and K, M or G)
#   sudo scripts/enable-swap.sh 8G /mnt/disk/swapfile
#
# It does four things: makes the file (fallocate, else dd), sets mode 600, runs mkswap and swapon, adds one line to
# /etc/fstab so it comes back after a reboot, and sets vm.swappiness to 10 (use swap only under pressure) in
# /etc/sysctl.d/99-planetgen-swap.conf.

set -euo pipefail

SIZE="${1:-8G}"
FILE="${2:-/swapfile}"
SYSCTL_FILE="/etc/sysctl.d/99-planetgen-swap.conf"

if [[ $EUID -ne 0 ]]; then
    echo "Run this as root: sudo $0 $*" >&2
    exit 1
fi

if ! [[ "$SIZE" =~ ^[0-9]+[KMGkmg]$ ]]; then
    echo "The size is a number and K, M or G (for example 8G), not '$SIZE'." >&2
    exit 1
fi

if swapon --show=NAME --noheadings | grep -q .; then
    echo "Swap is already on:"
    swapon --show
    if ! grep -qE '^[^#]*[[:space:]]swap[[:space:]]' /etc/fstab; then
        echo "Warning: no swap line in /etc/fstab, so it will not come back after a reboot." >&2
        echo "Add one for the file or partition shown above, or turn it off (swapoff -a) and run this again." >&2
    fi
    exit 0
fi

# The file's folder needs the room: the size, and a little over for the filesystem.
DIR="$(dirname "$FILE")"
mkdir -p "$DIR"
need_kb=$(numfmt --from=iec --to-unit=1024 "${SIZE^^}")
free_kb=$(df --output=avail -k "$DIR" | tail -n 1 | tr -d ' ')
if (( free_kb < need_kb + 1048576 )); then
    echo "Not enough room in $DIR: $((free_kb / 1024)) MB free, ${SIZE} wanted plus 1 GB to spare." >&2
    exit 1
fi

if [[ -e "$FILE" ]]; then
    echo "$FILE already exists but is not in use as swap; leaving it alone. Name another file or remove it first." >&2
    exit 1
fi

echo "Making $FILE ($SIZE)..."
if ! fallocate -l "$SIZE" "$FILE" 2>/dev/null; then
    # Some filesystems (and copy-on-write ones) cannot fallocate a swap file: write it out instead.
    dd if=/dev/zero of="$FILE" bs=1M count="$(( need_kb / 1024 ))" status=none
fi
chmod 600 "$FILE"
if ! mkswap "$FILE" >/dev/null; then
    rm -f "$FILE"
    echo "mkswap failed on $FILE (the filesystem may not allow swap files); nothing was changed." >&2
    exit 1
fi
swapon "$FILE"

if ! grep -qF "$FILE" /etc/fstab; then
    printf '%s none swap sw 0 0\n' "$FILE" >> /etc/fstab
    echo "Added $FILE to /etc/fstab, so it stays on after a reboot."
fi

printf 'vm.swappiness=10\n' > "$SYSCTL_FILE"
sysctl -q -p "$SYSCTL_FILE" || true

echo "Swap is on:"
swapon --show
free -m
