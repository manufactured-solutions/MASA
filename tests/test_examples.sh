#!/bin/bash
#
# This script runs each example in the examples directory
#
# It is designed to check that the masa api 
# has not changed
#
# Return 0 for success
#        1 upon failure

err=0

echo " "
echo "-------------------------------------------------------"
echo "Initializing MASA Examples Tests"
echo "-------------------------------------------------------"

# move to example directory
cd ../examples

# cpp
 ./heat_example
 err=$(($err+$?))

 ./euler_example
 err=$(($err+$?))

 ./euler_transient
 err=$(($err+$?))

 ./euler_chem
 err=$(($err+$?))

 ./navierstokes_example
 err=$(($err+$?))

 ./rans_sa_example
 err=$(($err+$?))

# # c
 ./c_euler_example
 err=$(($err+$?))

# # f90
 ./f_cns
 err=$(($err+$?))

# return error code 
exit $err

#
# written by: Nicholas Malaya
#             2/21/2011