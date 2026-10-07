#!/bin/sh
# Builds vendor/libmeshoptimizer.dylib from the pinned meshoptimizer checkout (needs clang++).
set -e
cd "$(dirname "$0")/vendor"
[ -d meshoptimizer ] || git clone --depth 1 https://github.com/zeux/meshoptimizer.git
cd meshoptimizer
echo "meshoptimizer $(git rev-parse --short HEAD)"
clang++ -O2 -std=c++11 -fPIC -shared -o ../libmeshoptimizer.dylib src/*.cpp
echo "built vendor/libmeshoptimizer.dylib"
