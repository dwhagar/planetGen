**Headline:** When update.sh cannot reload Apache or start Redis, it now prints why.

### Changed
- **update.sh keeps the reason for a failed Apache, systemctl or Redis step (OPS.166).** The Apache configuration test's output, the output of a failed `systemctl restart/reload` and of the Redis install and start now appear, indented, under the warning that names the failure. Before, it was thrown away. generation-output: unchanged
