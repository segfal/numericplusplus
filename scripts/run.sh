#!/bin/bash

# Run script - cleans, builds, runs, and then cleans again if successful

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="$PROJECT_ROOT/build"

# Clean and build first
echo "Cleaning and building..."
"$SCRIPT_DIR/build.sh"

if [ $? -ne 0 ]; then
    echo "Build failed. Exiting."
    exit 1
fi

# Run the executable
echo "Running executable..."
cd "$BUILD_DIR"
./NumericPlusPlus

RUN_EXIT_CODE=$?

if [ $RUN_EXIT_CODE -eq 0 ]; then
    echo "Execution successful!"
    # Clean again if successful
    echo "Cleaning up..."
    "$SCRIPT_DIR/clean.sh"
    echo "Cleanup complete."
else
    echo "Execution failed with exit code: $RUN_EXIT_CODE"
    exit $RUN_EXIT_CODE
fi

