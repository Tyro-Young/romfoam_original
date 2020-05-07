#!/bin/bash

echo "Regression test: Linear ROM..."
./AllrunLinear.sh > reg_file_linearROM
python check_values.py 0.402314020616319 1.0e-6 

# clean up
./Allclean.sh
sleep 1

echo "Regression test: Nonlinear ROM..."
./AllrunNonlinear.sh > reg_file_nonlinearROM
python check_values.py 0.402261736832066 1.0e-6
