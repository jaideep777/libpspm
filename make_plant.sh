#!/bin/bash 

cd tests/plant
g++ -O3 -std=c++11 -c plant.cpp
g++ -O3 -std=c++11 -c environment.cpp
cd ../..
