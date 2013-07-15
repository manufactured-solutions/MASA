#!/bin/bash
#
# hack to invoke python script
#

python test_swig.py
if($? != 0)
    echo 'SWIG interfaces:: FAILED'
    return 1
else
    return 0 
fi

#
# nick 
# 7/15/13
#