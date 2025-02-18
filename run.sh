#! /bin/bash
g++ main.cpp ./include/SOPMatrix.hpp ./include/SOPModel.hpp ./include/SOPlpsolver.hpp -o main
./main
