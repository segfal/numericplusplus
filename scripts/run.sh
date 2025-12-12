#!/bin/bash

# Run script - cleans, builds, runs, and then cleans again if successful
# USAGE: ./run.sh <EXECUTABLE_NAME>

if [ -z "$1" ]; then
    echo "Error: Executable name not provided."
    echo "Usage: ./run.sh <EXECUTABLE_NAME>"
    exit 1
fi

EXECUTABLE_NAME=$1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="$PROJECT_ROOT/build/$EXECUTABLE_NAME"

# Clean and build first
echo "Cleaning and building..."
"$SCRIPT_DIR/build.sh"

if [ $? -ne 0 ]; then
    echo "Build failed. Exiting."
    exit 1
fi

# Run the executable
echo "Running executable: $EXECUTABLE_NAME"
cd "$BUILD_DIR"

if [ ! -f "./$EXECUTABLE_NAME" ]; then
    echo "Error: Executable '$EXECUTABLE_NAME' not found in $BUILD_DIR. Did the build fail?"
    exit 1
fi

"./$EXECUTABLE_NAME"

RUN_EXIT_CODE=$?

if [ $RUN_EXIT_CODE -eq 0 ]; then
    echo "Execution successful!"
    # Keep the output file for plotting, so skip cleanup.
else
    echo "Execution failed with exit code: $RUN_EXIT_CODE"
    exit $RUN_EXIT_CODE
fi

