#!/bin/bash
#
# hack to invoke python script
#

python test_swig.py
if [ $? != 0 ] 
then
    echo 'Test SWIG interfaces:: FAILED'
    exit 1
else
    exit 0
fi

#
# nick 
# 7/15/13
#