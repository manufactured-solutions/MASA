#!/bin/bash
# this script will revert your files after running a masa_import

echo "Reverting back before masa_import.pl execution"
cp backup/Makefile.am ../
cp backup/masa_core.cpp ../
cp backup/masa_internal.h ../

#error handling:                                                                
if [ $? != 0 ]; then
    echo "FAILED: revert not completed"
    exit 1
else
    echo "Done"
fi

# nick
# 10/24/11