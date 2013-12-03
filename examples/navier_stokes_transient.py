import sys

# change this to your masa module path, 
# or ensure masa exists on your python path
# 
sys.path.append("/scratch/nick/masa_install/lib/python2.7/site-packages/masa")

import masa

# no errors
err  = 0

#
# Source terms available from: 
# https://github.com/manufactured-solutions/analytical/tree/master/navier_stokes/navier_stokes_Sutherland_transient/latex
#

err += masa.masa_init("3d Navier Stokes transient sutherland","navierstokes_3d_transient_sutherland")
err += masa.masa_sanity_check()

# evaluate source terms at some point in spacetime
x = 1.1
y = 1.3
z = 0.8
t = 1.2

# source terms
rho_field  = masa.masa_eval_4d_source_rho(x,y,z,t) # rho
rhou_field = masa.masa_eval_4d_source_u  (x,y,z,t) # rho*u
rhov_field = masa.masa_eval_4d_source_v  (x,y,z,t) # rho*v
rhow_field = masa.masa_eval_4d_source_w  (x,y,z,t) # rho*w
e_field    = masa.masa_eval_4d_source_e  (x,y,z,t) # e

# analytical terms
mms_rho_field  = masa.masa_eval_4d_exact_rho(x,y,z,t) # rho
mms_rhou_field = masa.masa_eval_4d_exact_u  (x,y,z,t) # rho*u
mms_rhov_field = masa.masa_eval_4d_exact_v  (x,y,z,t) # rho*v
mms_rhow_field = masa.masa_eval_4d_exact_w  (x,y,z,t) # rho*w
mms_p_field    = masa.masa_eval_4d_exact_p  (x,y,z,t) # pressure

if(err != 0):
    'Error encountered'
    sys.exit(1)
else:
    # steady as she goes
    sys.exit(0)

#                                                                                   
# nick                                                                              
# 12/3/13                                                                           
#      
