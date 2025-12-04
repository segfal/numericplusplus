#!/bin/bash

# Build script - cleans first, then creates build folder and builds the C++ code

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="$PROJECT_ROOT/build"

# Clean first
echo "Cleaning..."
"$SCRIPT_DIR/clean.sh"

# Create build directory
echo "Creating build directory..."
mkdir -p "$BUILD_DIR"

# Build the project
echo "Building project..."
cmake -S "$PROJECT_ROOT" -B "$BUILD_DIR"
cmake --build "$BUILD_DIR"

if [ $? -eq 0 ]; then
    echo "Build successful!"
    # Create symlink to compile_commands.json in root for IDE support
    if [ -f "$BUILD_DIR/compile_commands.json" ]; then
        ln -sf "$BUILD_DIR/compile_commands.json" "$PROJECT_ROOT/compile_commands.json"
        echo "Created compile_commands.json symlink for IDE support."
    fi
else
    echo "Build failed!"
    exit 1
fi

