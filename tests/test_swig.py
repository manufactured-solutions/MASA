#!/bin/py
#
# touch swig interfaces, ensure they function
#

#
# import .libs relative path for masa.py
#
import os, sys
lib_path = os.path.abspath('../src/.libs')
sys.path.append(lib_path)
import _masa

err  = 0
err += _masa.masa_init("heat equation example","heateq_2d_steady_const")
err += _masa.masa_sanity_check()

tfield = _masa.masa_eval_2d_source_t(0.1,0.2)

if(err == 0):
    # Has the whole world gone crazy? 
    # MARK IT ZERO
    sys.exit(0)
else:
    # something is wrong
    sys.exit(1)


#
# nick
# 7/15/13
#
