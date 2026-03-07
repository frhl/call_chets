#!/bin/bash

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Root directory is one level up from scripts/
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Get the version number from git tags
VERSION=$(git describe --tags --abbrev=0 2>/dev/null || echo "1.0.8")

# Use environment variables if available, otherwise try git commands
if [ -n "$GIT_COMMIT" ]; then
    COMMIT=$GIT_COMMIT
else
    COMMIT=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")
fi

if [ -n "$GIT_DATE" ]; then
    DATE=$GIT_DATE
else
    DATE=$(git log -1 --format=%cd --date=short 2>/dev/null || date +"%Y-%m-%d")
fi

# Write version to file in the root directory
echo "$VERSION" > "$ROOT_DIR/.version"

echo "Version set to $VERSION / commit = $COMMIT / release = $DATE"
