#!/bin/bash

# Get version from make_version.sh
VERSION=$(grep 'VERSION=' scripts/make_version.sh | cut -d'"' -f2)

# Get git information
GIT_COMMIT=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")
GIT_DATE=$(git log -1 --format=%cd --date=short 2>/dev/null || date +"%Y-%m-%d")

echo "Building Docker image fhlassen/call_chets:$VERSION with commit $GIT_COMMIT from $GIT_DATE"

# Build the Docker image with git information
docker build \
  --platform linux/amd64 \
  --build-arg GIT_COMMIT=$GIT_COMMIT \
  --build-arg GIT_DATE=$GIT_DATE \
  -t fhlassen/call_chets:$VERSION .

echo "Docker image fhlassen/call_chets:$VERSION built successfully"

# Push to DockerHub
echo "Pushing to DockerHub..."
docker push fhlassen/call_chets:$VERSION
echo "Pushed fhlassen/call_chets:$VERSION"

