#!/bin/bash

# Get the version number from major.minor.patch
VERSION="0.3.3"

# Get the git commit hash
COMMIT=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")

# Get the date
DATE=$(git log -1 --format=%cd --date=short 2>/dev/null || date +"%Y-%m-%d")

# Write version to file
echo "$VERSION" > .version

echo "Version set to $VERSION / commit = $COMMIT / release = $DATE"
