#!/bin/bash
echo "Current working directory:"
pwd
echo "Make bin directory"
mkdir -p bin
mkdir -p logs
echo "Check formatting/layout of source files"
python make.py 1
#mpif90 -O3 -o 
echo "Compile source"
mpif90 -O3 -g -fcheck=all -Wall -mtune=native \
-o bin/activitylite.x src/activitylite.f90 &> logs/make.log
sleep 1 
echo "Check compile result"
python make.py 2