#!/usr/bin/env bash
set -euo pipefail

RELEASE_URL="https://github.com/dial481/allelix/releases/download/v1.1.1/test_data.tar.gz"
DEST="test_data"

if [ -d "$DEST/real" ] && [ -d "$DEST/transcoded" ]; then
    echo "Test data already present."
    exit 0
fi

echo "Downloading test data..."
curl -L -o test_data.tar.gz "$RELEASE_URL"
tar -xzf test_data.tar.gz
rm test_data.tar.gz
echo "Done. $(find test_data/real -type f | wc -l) files extracted."
