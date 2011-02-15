!! -*-f90-*-
!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! MASA - Manufactured Analytical Solutions Abstraction Library
!!
!! Copyright (C) 2010 The PECOS Development Team
!!
!! This library is free software; you can redistribute it and/or
!! modify it under the terms of the Version 2.1 GNU Lesser General
!! Public License as published by the Free Software Foundation.
!!
!! This library is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!! Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public
!! License along with this library; if not, write to the Free Software
!! Foundation, Inc. 51 Franklin Street, Fifth Floor,
!! Boston, MA  02110-1301  USA
!!
!!-----------------------------------------------------------------------el-
!!
!! f_euler1d.f90: 
!!
!! $Id: 
!! -------------------------------------------------------------------------
!! -------------------------------------------------------------------------
function K(T)
  implicit none
  
  real(8),intent(in) :: T
  real(8)            :: K

  !! hackish functional here
  !! This is an eyeballed fit (focusing on the 5000K-6000K range) 
  !! for the equilibrium constant for N2->N+N dissociation

  K = exp(4+(T-6000)/500)

end subroutine temp_function


program main
  use euler_source_interface
  use masa
  implicit none

  real(8) :: MASA_DEFAULT = -12345.67d0;

  real(8) :: u_0
  real(8) :: u_x
  real(8) :: a_ux
  real(8) :: L

  real(8) :: R
  real(8) :: Cf1_N
  real(8) :: Cf1_N2
  real(8) :: etaf1_N
  real(8) :: etaf1_N2
  real(8) :: Ea_N
  real(8) :: Ea_N2
  real(8) :: R_N
  real(8) :: R_N2
  real(8) :: theta_v_N2
  real(8) :: M_N
  real(8) :: h0_N
  real(8) :: h0_N2
  real(8) :: K

  real(8) :: rho_N_0
  real(8) :: rho_N_x
  real(8) :: a_rho_N_x

  real(8) :: rho_N2_0
  real(8) :: rho_N2_x
  real(8) :: a_rho_N2_x

  real(8) :: T_0
  real(8) :: T_x
  real(8) :: a_Tx

  !problem size
  integer :: nx
  integer :: lx  
  real(8) ::  x
  real(8) :: dx

  ! solutions
  real(8) :: ufield,ufield2,ufield3;
  real(8) :: efield,efield2,efield3;
  real(8) :: N,N2,N3;
  real(8) :: Ntwo,Ntwo2,Ntwo3;
  real(8) :: gradx,grady,gradz,gradp,gradrho;

  real(8) :: exact_t,exact_t2,exact_t3;
  real(8) :: exact_u,exact_u2,exact_u3;
  real(8) :: exact_rho,exact_rho2,exact_rho3;
  real(8) :: exact_N,exact_N2,exact_N3;
  real(8) :: exact_Ntwo,exact_Ntwo2,exact_Ntwo3;

  ! iterators
  integer :: i

  nx = 200
  lx = 10

  ! initialize the problem
  dx = real(lx)/real(nx)

  call masa_init("mytest","euler_chem_1d")

  ! initialize the default parameters
  call masa_init_param()

  ! get defaults for comparison to source terms
  u_0  = masa_get_param("u_0");
  u_x  = masa_get_param("u_x");
  a_ux = masa_get_param("a_ux");
  L    = masa_get_param("L");
  R    = masa_get_param("R");

  Cf1_N   = masa_get_param("Cf1_N");
  Cf1_N2  = masa_get_param("Cf1_N2");
  etaf1_N  = masa_get_param("etaf1_N");
  etaf1_N2 = masa_get_param("etaf1_N2");

  Ea_N  = masa_get_param("Ea_N");
  Ea_N2 = masa_get_param("Ea_N2");

  R_N   = masa_get_param("R_N");
  R_N2  = masa_get_param("R_N2");

  theta_v_N2 = masa_get_param("theta_v_N2");
  M_N   = masa_get_param("M_N");
  h0_N  = masa_get_param("h0_N");
  h0_N2 = masa_get_param("h0_N2");
  K     = masa_get_param("K");

  rho_N_0   = masa_get_param("rho_N_0");
  rho_N_x   = masa_get_param("rho_N_x");
  a_rho_N_x = masa_get_param("a_rho_N_x");

  rho_N2_0   = masa_get_param("rho_N2_0");  
  rho_N2_x   = masa_get_param("rho_N2_x");
  a_rho_N2_x = masa_get_param("a_rho_N2_x");
  
  T_0  = masa_get_param("T_0");
  T_x  = masa_get_param("T_x");
  a_Tx = masa_get_param("a_Tx");

  ! call the sanity check routine 
  ! (tests that all variables have been initialized)
  call masa_sanity_check()

  ! evaluate source terms (2D)
  do i=0, nx

     x = i * dx
     
     ! evalulate source terms
     ufield = masa_eval_1d_source_rho_u (x);
     efield = masa_eval_1d_source_rho_e (x);
     N      = masa_eval_1d_source_rho_N (x);
     Ntwo   = masa_eval_1d_source_rho_N2(x);
     
     ! evaluate analytical solution terms
     exact_t    = masa_eval_1d_exact_t     (x);
     exact_u    = masa_eval_1d_exact_u     (x);
     exact_rho  = masa_eval_1d_exact_rho   (x);
     exact_N    = masa_eval_1d_exact_rho_N (x);
     exact_Ntwo = masa_eval_1d_exact_rho_N2(x);
     
  enddo

  ! steady as she goes...
  call exit(0)

end program main
