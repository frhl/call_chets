#!/bin/bash

# Get the version number from major.minor.patch
VERSION="1.0.4"

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

# Write version to file
echo "$VERSION" > .version

echo "Version set to $VERSION / commit = $COMMIT / release = $DATE"
