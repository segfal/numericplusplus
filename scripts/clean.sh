#!/bin/bash

# Clean script - removes the build directory if it exists
# If it doesn't exist, it ignores it

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="$PROJECT_ROOT/build"

if [ -d "$BUILD_DIR" ]; then
    echo "Removing build directory..."
    rm -rf "$BUILD_DIR"
    echo "Build directory removed."
else
    echo "Build directory does not exist. Nothing to clean."
fi

# Remove compile_commands.json symlink if it exists
if [ -L "$PROJECT_ROOT/compile_commands.json" ]; then
    echo "Removing compile_commands.json symlink..."
    rm -f "$PROJECT_ROOT/compile_commands.json"
    echo "Symlink removed."
fi

