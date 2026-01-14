#!/bin/bash
# =============================================================================
# entrypoint.sh - Docker entrypoint for processVCF container
# =============================================================================
# Handles permission issues when mounting host directories.
# Runs as root initially to fix permissions, then drops to 'user'.
# =============================================================================

set -e

# Target user/group (created in Dockerfile)
TARGET_USER="user"
TARGET_GROUP="user"

# If running as root, fix permissions and switch to user
if [ "$(id -u)" = "0" ]; then
    # Ensure /data directory exists and is writable by user
    if [ -d "/data" ]; then
        chown "$TARGET_USER:$TARGET_GROUP" /data 2>/dev/null || true

        # Also fix ownership of any existing subdirectories that need write access
        # (output, annotationTMSP, annotationCEBNX, SnapShots, html_reports, IgvBed)
        for subdir in output vcf cebpa; do
            if [ -d "/data/$subdir" ]; then
                chown -R "$TARGET_USER:$TARGET_GROUP" "/data/$subdir" 2>/dev/null || true
            fi
        done
    fi

    # Execute command as user using gosu
    exec gosu "$TARGET_USER" "$@"
fi

# If already running as non-root, just execute the command
exec "$@"
