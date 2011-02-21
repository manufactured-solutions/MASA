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

# start tests
./heat_example
err=$(($err+$?))

exit $err

#
# written by: Nicholas Malaya
#             2/21/2011