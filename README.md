# Simplex method
A well-known method for solving linear programming, simplex method.

For implementing issue, I used C++ with Eigen Library to do complicate linear algebra more efficient.

The ultimate goal is to build a geneal MIP Solver.

In current version, only include "primal simplex method" & "dual simplex method"
During the developing, I am still studying some classical textbook to strengthen my professional knowledge about discrete & combinatorial opitmization.

## Dependency
1. Eigen3+
```cmd
git clone https://gitlab.com/libeigen/eigen.git
cd eigen
mkdir build
cd build
cmake ..
sudo make install
sudo cp -r /usr/local/include/eigen3/Eigen usr/local/include
```	

## How to use?
```cmd
git clone https://github.com/ytsao/simplex-method.git
cd simplex-method 
mkdir build 
cd build 
cmake .. 
make 
./simplex_method
```
