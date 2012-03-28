//--------------------------------------------------------------------------
//
// Copyright (C) 2011,2012 Rhys Ulerich
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at http://mozilla.org/MPL/2.0/.
//
// This is one of two files implementing the manufactured solution from
// @inproceedings{Ulerich2012,
//   title={A transient manufactured solution for the compressible
//          {N}avier–-{S}tokes equations with a power law viscosity},
//   author={Rhys Ulerich and Kemelli C. Estacio-Hiroms
//           and Nicholas Malaya and Robert D. Moser},
//   booktitle={{10th World Congress on Computational Mechanics}},
//   address={{S}\~{a}o {P}aulo, {B}razil},
//   month={{J}uly},
//   year={2012}
// }
// Please refer to the paper for more details on the manufactured solution.
// Please consider citing the paper if you find the solution useful.
//
//--------------------------------------------------------------------------

#include "nsctpl_fwd.hpp"

#ifndef NSCTPL_HPP
#define NSCTPL_HPP

namespace nsctpl {

// ---------------------------------------------------------------------------
// primitive<Scalar> member implementations
// ---------------------------------------------------------------------------

template <typename Scalar>
const Scalar primitive<Scalar>::twopi = 8 * std::atan(Scalar(1));

template <typename Scalar>
template <typename T1, typename T2, typename T3, typename T4>
Scalar primitive<Scalar>::operator()(T1 x, T2 y, T3 z, T4 t) const
{
    using ::std::cos;
    using ::std::sin;
    const Scalar twopi_invLx = twopi / *Lx;
    const Scalar twopi_invLy = twopi / *Ly;
    const Scalar twopi_invLz = twopi / *Lz;

    return a_0*cos(g_0 + f_0*t) + a_x*cos(c_x + b_x*twopi_invLx*x)*cos(g_x + f_x*t) + a_y*cos(g_y + f_y*t)*cos(c_y + b_y*twopi_invLy*y) + a_z*cos(g_z + f_z*t)*cos(c_z + b_z*twopi_invLz*z) + a_xy*cos(g_xy + f_xy*t)*cos(c_xy + b_xy*twopi_invLx*x)*cos(e_xy + d_xy*twopi_invLy*y) + a_xz*cos(g_xz + f_xz*t)*cos(c_xz + b_xz*twopi_invLx*x)*cos(e_xz + d_xz*twopi_invLz*z) + a_yz*cos(c_yz + b_yz*twopi_invLy*y)*cos(e_yz + d_yz*twopi_invLz*z)*cos(g_yz + f_yz*t);
}

template <typename Scalar>
template <typename T1, typename T2, typename T3, typename T4>
Scalar primitive<Scalar>::_t(T1 x, T2 y, T3 z, T4 t) const
{
    using ::std::cos;
    using ::std::sin;
    const Scalar twopi_invLx = twopi / *Lx;
    const Scalar twopi_invLy = twopi / *Ly;
    const Scalar twopi_invLz = twopi / *Lz;

    return -a_0*f_0*sin(g_0 + f_0*t) - a_x*f_x*cos(c_x + b_x*twopi_invLx*x)*sin(g_x + f_x*t) - a_y*f_y*cos(c_y + b_y*twopi_invLy*y)*sin(g_y + f_y*t) - a_z*f_z*cos(c_z + b_z*twopi_invLz*z)*sin(g_z + f_z*t) - a_xy*f_xy*cos(c_xy + b_xy*twopi_invLx*x)*cos(e_xy + d_xy*twopi_invLy*y)*sin(g_xy + f_xy*t) - a_xz*f_xz*cos(c_xz + b_xz*twopi_invLx*x)*cos(e_xz + d_xz*twopi_invLz*z)*sin(g_xz + f_xz*t) - a_yz*f_yz*cos(c_yz + b_yz*twopi_invLy*y)*cos(e_yz + d_yz*twopi_invLz*z)*sin(g_yz + f_yz*t);
}

template <typename Scalar>
template <typename T1, typename T2, typename T3, typename T4>
Scalar primitive<Scalar>::_x(T1 x, T2 y, T3 z, T4 t) const
{
    using ::std::cos;
    using ::std::sin;
    const Scalar twopi_invLx = twopi / *Lx;
    const Scalar twopi_invLy = twopi / *Ly;
    const Scalar twopi_invLz = twopi / *Lz;

    return -a_x*b_x*twopi_invLx*cos(g_x + f_x*t)*sin(c_x + b_x*twopi_invLx*x) - a_xy*b_xy*twopi_invLx*cos(g_xy + f_xy*t)*cos(e_xy + d_xy*twopi_invLy*y)*sin(c_xy + b_xy*twopi_invLx*x) - a_xz*b_xz*twopi_invLx*cos(g_xz + f_xz*t)*cos(e_xz + d_xz*twopi_invLz*z)*sin(c_xz + b_xz*twopi_invLx*x);
}

template <typename Scalar>
template <typename T1, typename T2, typename T3, typename T4>
Scalar primitive<Scalar>::_xx(T1 x, T2 y, T3 z, T4 t) const
{
    using ::std::cos;
    using ::std::sin;
    const Scalar twopi_invLx = twopi / *Lx;
    const Scalar twopi_invLy = twopi / *Ly;
    const Scalar twopi_invLz = twopi / *Lz;

    return -a_x*b_x*b_x*twopi_invLx*twopi_invLx*cos(c_x + b_x*twopi_invLx*x)*cos(g_x + f_x*t) - a_xy*b_xy*b_xy*twopi_invLx*twopi_invLx*cos(g_xy + f_xy*t)*cos(c_xy + b_xy*twopi_invLx*x)*cos(e_xy + d_xy*twopi_invLy*y) - a_xz*b_xz*b_xz*twopi_invLx*twopi_invLx*cos(g_xz + f_xz*t)*cos(c_xz + b_xz*twopi_invLx*x)*cos(e_xz + d_xz*twopi_invLz*z);
}

template <typename Scalar>
template <typename T1, typename T2, typename T3, typename T4>
Scalar primitive<Scalar>::_xy(T1 x, T2 y, T3 z, T4 t) const
{
    using ::std::cos;
    using ::std::sin;
    const Scalar twopi_invLx = twopi / *Lx;
    const Scalar twopi_invLy = twopi / *Ly;

    return a_xy*b_xy*d_xy*twopi_invLx*twopi_invLy*cos(g_xy + f_xy*t)*sin(c_xy + b_xy*twopi_invLx*x)*sin(e_xy + d_xy*twopi_invLy*y);
}

template <typename Scalar>
template <typename T1, typename T2, typename T3, typename T4>
Scalar primitive<Scalar>::_xz(T1 x, T2 y, T3 z, T4 t) const
{
    using ::std::cos;
    using ::std::sin;
    const Scalar twopi_invLx = twopi / *Lx;
    const Scalar twopi_invLz = twopi / *Lz;

    return a_xz*b_xz*d_xz*twopi_invLx*twopi_invLz*cos(g_xz + f_xz*t)*sin(c_xz + b_xz*twopi_invLx*x)*sin(e_xz + d_xz*twopi_invLz*z);
}

template <typename Scalar>
template <typename T1, typename T2, typename T3, typename T4>
Scalar primitive<Scalar>::_y(T1 x, T2 y, T3 z, T4 t) const
{
    using ::std::cos;
    using ::std::sin;
    const Scalar twopi_invLx = twopi / *Lx;
    const Scalar twopi_invLy = twopi / *Ly;
    const Scalar twopi_invLz = twopi / *Lz;

    return -a_y*b_y*twopi_invLy*cos(g_y + f_y*t)*sin(c_y + b_y*twopi_invLy*y) - a_xy*d_xy*twopi_invLy*cos(g_xy + f_xy*t)*cos(c_xy + b_xy*twopi_invLx*x)*sin(e_xy + d_xy*twopi_invLy*y) - a_yz*b_yz*twopi_invLy*cos(e_yz + d_yz*twopi_invLz*z)*cos(g_yz + f_yz*t)*sin(c_yz + b_yz*twopi_invLy*y);
}

template <typename Scalar>
template <typename T1, typename T2, typename T3, typename T4>
Scalar primitive<Scalar>::_yy(T1 x, T2 y, T3 z, T4 t) const
{
    using ::std::cos;
    using ::std::sin;
    const Scalar twopi_invLx = twopi / *Lx;
    const Scalar twopi_invLy = twopi / *Ly;
    const Scalar twopi_invLz = twopi / *Lz;

    return -a_y*b_y*b_y*twopi_invLy*twopi_invLy*cos(g_y + f_y*t)*cos(c_y + b_y*twopi_invLy*y) - a_xy*d_xy*d_xy*twopi_invLy*twopi_invLy*cos(g_xy + f_xy*t)*cos(c_xy + b_xy*twopi_invLx*x)*cos(e_xy + d_xy*twopi_invLy*y) - a_yz*b_yz*b_yz*twopi_invLy*twopi_invLy*cos(c_yz + b_yz*twopi_invLy*y)*cos(e_yz + d_yz*twopi_invLz*z)*cos(g_yz + f_yz*t);
}

template <typename Scalar>
template <typename T1, typename T2, typename T3, typename T4>
Scalar primitive<Scalar>::_yz(T1 x, T2 y, T3 z, T4 t) const
{
    using ::std::cos;
    using ::std::sin;
    const Scalar twopi_invLy = twopi / *Ly;
    const Scalar twopi_invLz = twopi / *Lz;

    return a_yz*b_yz*d_yz*twopi_invLy*twopi_invLz*cos(g_yz + f_yz*t)*sin(c_yz + b_yz*twopi_invLy*y)*sin(e_yz + d_yz*twopi_invLz*z);
}

template <typename Scalar>
template <typename T1, typename T2, typename T3, typename T4>
Scalar primitive<Scalar>::_z(T1 x, T2 y, T3 z, T4 t) const
{
    using ::std::cos;
    using ::std::sin;
    const Scalar twopi_invLx = twopi / *Lx;
    const Scalar twopi_invLy = twopi / *Ly;
    const Scalar twopi_invLz = twopi / *Lz;

    return -a_z*b_z*twopi_invLz*cos(g_z + f_z*t)*sin(c_z + b_z*twopi_invLz*z) - a_xz*d_xz*twopi_invLz*cos(g_xz + f_xz*t)*cos(c_xz + b_xz*twopi_invLx*x)*sin(e_xz + d_xz*twopi_invLz*z) - a_yz*d_yz*twopi_invLz*cos(c_yz + b_yz*twopi_invLy*y)*cos(g_yz + f_yz*t)*sin(e_yz + d_yz*twopi_invLz*z);
}

template <typename Scalar>
template <typename T1, typename T2, typename T3, typename T4>
Scalar primitive<Scalar>::_zz(T1 x, T2 y, T3 z, T4 t) const
{
    using ::std::cos;
    using ::std::sin;
    const Scalar twopi_invLx = twopi / *Lx;
    const Scalar twopi_invLy = twopi / *Ly;
    const Scalar twopi_invLz = twopi / *Lz;

    return -a_z*b_z*b_z*twopi_invLz*twopi_invLz*cos(g_z + f_z*t)*cos(c_z + b_z*twopi_invLz*z) - a_xz*d_xz*d_xz*twopi_invLz*twopi_invLz*cos(g_xz + f_xz*t)*cos(c_xz + b_xz*twopi_invLx*x)*cos(e_xz + d_xz*twopi_invLz*z) - a_yz*d_yz*d_yz*twopi_invLz*twopi_invLz*cos(c_yz + b_yz*twopi_invLy*y)*cos(e_yz + d_yz*twopi_invLz*z)*cos(g_yz + f_yz*t);
}

// ---------------------------------------------------------------------------
// manufactured_solution<Scalar,IndexBase,Primitive> analytical member
// implementations See comments in nsctpl_fwd.hpp regarding SWIG-related
// #ifdefs
// ---------------------------------------------------------------------------

#ifndef SWIG
#define NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION     \
        template <typename Scalar,                     \
                  int IndexBase,                       \
                  template <typename> class Primitive>
#define NSCTPL_MANUFACTURED_SOLUTION_TPNAMES           \
        Scalar,IndexBase,Primitive
#else
#define NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION     \
        template <typename Scalar,                     \
                  int IndexBase>
#define NSCTPL_MANUFACTURED_SOLUTION_TPNAMES           \
        Scalar,IndexBase
#endif

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
const int manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::index_base = IndexBase;

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::grad_rho(
        T1 x, T2 y, T3 z, T4 t, int index) const
{
    switch (index) {
        case (IndexBase + 0):  return this->rho._x(x, y, z, t);
        case (IndexBase + 1):  return this->rho._y(x, y, z, t);
        case (IndexBase + 2):  return this->rho._z(x, y, z, t);
        default: return ::std::numeric_limits<Scalar>::signaling_NaN();
    }
}

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::grad_u(
        T1 x, T2 y, T3 z, T4 t, int index) const
{
    switch (index) {
        case (IndexBase + 0):  return this->u._x(x, y, z, t);
        case (IndexBase + 1):  return this->u._y(x, y, z, t);
        case (IndexBase + 2):  return this->u._z(x, y, z, t);
        default: return ::std::numeric_limits<Scalar>::signaling_NaN();
    }
}

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::grad_v(
        T1 x, T2 y, T3 z, T4 t, int index) const
{
    switch (index) {
        case (IndexBase + 0):  return this->v._x(x, y, z, t);
        case (IndexBase + 1):  return this->v._y(x, y, z, t);
        case (IndexBase + 2):  return this->v._z(x, y, z, t);
        default: return ::std::numeric_limits<Scalar>::signaling_NaN();
    }
}

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::grad_w(
        T1 x, T2 y, T3 z, T4 t, int index) const
{
    switch (index) {
        case (IndexBase + 0):  return this->w._x(x, y, z, t);
        case (IndexBase + 1):  return this->w._y(x, y, z, t);
        case (IndexBase + 2):  return this->w._z(x, y, z, t);
        default: return ::std::numeric_limits<Scalar>::signaling_NaN();
    }
}

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::grad_T(
        T1 x, T2 y, T3 z, T4 t, int index) const
{
    switch (index) {
        case (IndexBase + 0):  return this->T._x(x, y, z, t);
        case (IndexBase + 1):  return this->T._y(x, y, z, t);
        case (IndexBase + 2):  return this->T._z(x, y, z, t);
        default: return ::std::numeric_limits<Scalar>::signaling_NaN();
    }
}

// ---------------------------------------------------------------------------
// manufactured_solution<Scalar,IndexBase,Primitive> member
// implementations derived from analytical results
// ---------------------------------------------------------------------------

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::e(
        T1 x, T2 y, T3 z, T4 t) const
{
    const Scalar T = this->T(x, y, z, t);
    const Scalar u = this->u(x, y, z, t);
    const Scalar v = this->v(x, y, z, t);
    const Scalar w = this->w(x, y, z, t);
    const Scalar e = R * T   / (gamma - 1) + (u*u   + v*v   + w*w  ) / 2;
    return e;
}

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::p(
        T1 x, T2 y, T3 z, T4 t) const
{
    const Scalar T   = this->T(x, y, z, t);
    const Scalar rho = this->rho(x, y, z, t); // shadow
    const Scalar p   = rho * R * T;
    return p;
}

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::mu(
        T1 x, T2 y, T3 z, T4 t) const
{
    const Scalar T  = this->T(x, y, z, t);
    const Scalar mu = mu_r * ::std::pow(T / T_r, beta);
    return mu;
}

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::rhou(
        T1 x, T2 y, T3 z, T4 t) const
{
    return rho(x, y, z, t) * u(x, y, z, t);
}

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::rhov(
        T1 x, T2 y, T3 z, T4 t) const
{
    return rho(x, y, z, t) * v(x, y, z, t);
}

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::rhow(
        T1 x, T2 y, T3 z, T4 t) const
{
    return rho(x, y, z, t) * w(x, y, z, t);
}

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::rhoe(
        T1 x, T2 y, T3 z, T4 t) const
{
    return rho(x, y, z, t) * e(x, y, z, t);
}

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::grad_e(
        T1 x, T2 y, T3 z, T4 t, int index) const
{
    const Scalar u = this->u(x, y, z, t);
    const Scalar v = this->v(x, y, z, t);
    const Scalar w = this->w(x, y, z, t);
    Scalar dT, du, dv, dw;
    switch (index) {
        case (IndexBase + 0):
            dT = this->T._x(x, y, z, t);
            du = this->u._x(x, y, z, t);
            dv = this->v._x(x, y, z, t);
            dw = this->w._x(x, y, z, t);
            break;
        case (IndexBase + 1):
            dT = this->T._y(x, y, z, t);
            du = this->u._y(x, y, z, t);
            dv = this->v._y(x, y, z, t);
            dw = this->w._y(x, y, z, t);
            break;
        case (IndexBase + 2):
            dT = this->T._z(x, y, z, t);
            du = this->u._z(x, y, z, t);
            dv = this->v._z(x, y, z, t);
            dw = this->w._z(x, y, z, t);
            break;
        default: return ::std::numeric_limits<Scalar>::signaling_NaN();
    }
    const Scalar de = R * dT / (gamma - 1) + (u*du + v*dv + w*dw);
    return de;
}

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::grad_p(
        T1 x, T2 y, T3 z, T4 t, int index) const
{
    const Scalar rho = this->rho(x, y, z, t); // shadow
    const Scalar T   = this->T(x, y, z, t);
    Scalar drho, dT;
    switch (index) {
        case (IndexBase + 0):
            dT   = this->T._x(x, y, z, t);
            drho = this->rho._x(x, y, z, t);
            break;
        case (IndexBase + 1):
            dT   = this->T._y(x, y, z, t);
            drho = this->rho._y(x, y, z, t);
            break;
        case (IndexBase + 2):
            dT   = this->T._z(x, y, z, t);
            drho = this->rho._z(x, y, z, t);
            break;
        default: return ::std::numeric_limits<Scalar>::signaling_NaN();
    }
    const Scalar dp = drho * R * T + rho * R * dT;
    return dp;
}

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::grad_mu(
        T1 x, T2 y, T3 z, T4 t, int index) const
{
    const Scalar T = this->T(x, y, z, t);
    Scalar dT;
    switch (index) {
        case (IndexBase + 0):
            dT = this->T._x(x, y, z, t);
            break;
        case (IndexBase + 1):
            dT = this->T._y(x, y, z, t);
            break;
        case (IndexBase + 2):
            dT = this->T._z(x, y, z, t);
            break;
        default: return ::std::numeric_limits<Scalar>::signaling_NaN();
    }
    const Scalar dmu = beta * mu_r / T_r * ::std::pow(T / T_r, beta - 1) * dT;
    return dmu;
}

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::Q_rho(
        T1 x, T2 y, T3 z, T4 t) const
{
    /* Computations retrieving primitive solution details */
    const Scalar rho   = this->rho   (x, y, z, t);  // shadow
    const Scalar rho_t = this->rho._t(x, y, z, t);
    const Scalar rho_x = this->rho._x(x, y, z, t);
    const Scalar rho_y = this->rho._y(x, y, z, t);
    const Scalar rho_z = this->rho._z(x, y, z, t);

    const Scalar u     = this->u   (x, y, z, t);    // shadow
    const Scalar u_x   = this->u._x(x, y, z, t);

    const Scalar v     = this->v   (x, y, z, t);    // shadow
    const Scalar v_y   = this->v._y(x, y, z, t);

    const Scalar w     = this->w   (x, y, z, t);    // shadow
    const Scalar w_z   = this->w._z(x, y, z, t);

    /* Computations stemming from the compressible, Newtonian fluid model */
    const Scalar rhou_x  = rho_x * u + rho * u_x;
    const Scalar rhov_y  = rho_y * v + rho * v_y;
    const Scalar rhow_z  = rho_z * w + rho * w_z;
    const Scalar Q_rho   = rho_t  + rhou_x + rhov_y + rhow_z;
    return Q_rho;
}

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::Q_rhou(
        T1 x, T2 y, T3 z, T4 t) const
{
    /* Computations retrieving primitive solution details */
    const Scalar rho   = this->rho   (x, y, z, t);  // shadow
    const Scalar rho_t = this->rho._t(x, y, z, t);
    const Scalar rho_x = this->rho._x(x, y, z, t);
    const Scalar rho_y = this->rho._y(x, y, z, t);
    const Scalar rho_z = this->rho._z(x, y, z, t);

    const Scalar u    = this->u    (x, y, z, t);    // shadow
    const Scalar u_t  = this->u._t (x, y, z, t);
    const Scalar u_x  = this->u._x (x, y, z, t);
    const Scalar u_xx = this->u._xx(x, y, z, t);
    const Scalar u_y  = this->u._y (x, y, z, t);
    const Scalar u_yy = this->u._yy(x, y, z, t);
    const Scalar u_z  = this->u._z (x, y, z, t);
    const Scalar u_zz = this->u._zz(x, y, z, t);

    const Scalar v    = this->v    (x, y, z, t);    // shadow
    const Scalar v_x  = this->v._x (x, y, z, t);
    const Scalar v_xy = this->v._xy(x, y, z, t);
    const Scalar v_y  = this->v._y (x, y, z, t);

    const Scalar w    = this->w    (x, y, z, t);    // shadow
    const Scalar w_x  = this->w._x (x, y, z, t);
    const Scalar w_xz = this->w._xz(x, y, z, t);
    const Scalar w_z  = this->w._z (x, y, z, t);

    const Scalar T    = this->T    (x, y, z, t);    // shadow
    const Scalar T_x  = this->T._x (x, y, z, t);
    const Scalar T_y  = this->T._y (x, y, z, t);
    const Scalar T_z  = this->T._z (x, y, z, t);

    /* Computations stemming from the constitutive relationships */
    const Scalar p_x      = rho_x * R * T + rho * R * T_x;
    const Scalar mu       = mu_r * ::std::pow(T / T_r, beta);
    const Scalar mu_x     = beta * mu_r / T_r * ::std::pow(T / T_r, beta - 1)
                          * T_x;
    const Scalar mu_y     = beta * mu_r / T_r * ::std::pow(T / T_r, beta - 1)
                          * T_y;
    const Scalar mu_z     = beta * mu_r / T_r * ::std::pow(T / T_r, beta - 1)
                          * T_z;
    const Scalar lambda   = lambda_r / mu_r * mu;
    const Scalar lambda_x = lambda_r / mu_r * mu_x;

    /* Computations stemming from the compressible, Newtonian fluid model */
    const Scalar rhou_t  = rho_t * u + rho * u_t;

    const Scalar rhouu_x = (rho_x * u * u) + (rho * u_x * u) + (rho * u * u_x);
    const Scalar rhouv_y = (rho_y * u * v) + (rho * u_y * v) + (rho * u * v_y);
    const Scalar rhouw_z = (rho_z * u * w) + (rho * u_z * w) + (rho * u * w_z);

    const Scalar tauxx_x = mu_x * (u_x  + u_x )
                         + lambda_x * (u_x  + v_y  + w_z )
                         + mu   * (u_xx + u_xx)
                         + lambda   * (u_xx + v_xy + w_xz);

    const Scalar tauxy_y = mu_y * (u_y + v_x) + mu * (u_yy + v_xy);
    const Scalar tauxz_z = mu_z * (u_z + w_x) + mu * (u_zz + w_xz);

    const Scalar Q_rhou = rhou_t + rhouu_x + rhouv_y + rhouw_z
                        + p_x - tauxx_x - tauxy_y - tauxz_z;

    return Q_rhou;
}

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::Q_rhov(
        T1 x, T2 y, T3 z, T4 t) const
{
    /* Computations retrieving primitive solution details */
    const Scalar rho   = this->rho   (x, y, z, t);  // shadow
    const Scalar rho_t = this->rho._t(x, y, z, t);
    const Scalar rho_x = this->rho._x(x, y, z, t);
    const Scalar rho_y = this->rho._y(x, y, z, t);
    const Scalar rho_z = this->rho._z(x, y, z, t);

    const Scalar u    = this->u    (x, y, z, t);    // shadow
    const Scalar u_x  = this->u._x (x, y, z, t);
    const Scalar u_xy = this->u._xy(x, y, z, t);
    const Scalar u_y  = this->u._y (x, y, z, t);

    const Scalar v    = this->v    (x, y, z, t);    // shadow
    const Scalar v_t  = this->v._t (x, y, z, t);
    const Scalar v_x  = this->v._x (x, y, z, t);
    const Scalar v_xx = this->v._xx(x, y, z, t);
    const Scalar v_y  = this->v._y (x, y, z, t);
    const Scalar v_yy = this->v._yy(x, y, z, t);
    const Scalar v_z  = this->v._z (x, y, z, t);
    const Scalar v_zz = this->v._zz(x, y, z, t);

    const Scalar w    = this->w    (x, y, z, t);    // shadow
    const Scalar w_y  = this->w._y (x, y, z, t);
    const Scalar w_yz = this->w._yz(x, y, z, t);
    const Scalar w_z  = this->w._z (x, y, z, t);

    const Scalar T    = this->T    (x, y, z, t);    // shadow
    const Scalar T_x  = this->T._x (x, y, z, t);
    const Scalar T_y  = this->T._y (x, y, z, t);
    const Scalar T_z  = this->T._z (x, y, z, t);

    /* Computations stemming from the constitutive relationships */
    const Scalar p_y      = rho_y * R * T + rho * R * T_y;
    const Scalar mu       = mu_r * ::std::pow(T / T_r, beta);
    const Scalar mu_x     = beta * mu_r / T_r * ::std::pow(T / T_r, beta - 1)
                          * T_x;
    const Scalar mu_y     = beta * mu_r / T_r * ::std::pow(T / T_r, beta - 1)
                          * T_y;
    const Scalar mu_z     = beta * mu_r / T_r * ::std::pow(T / T_r, beta - 1)
                          * T_z;
    const Scalar lambda   = lambda_r / mu_r * mu;
    const Scalar lambda_y = lambda_r / mu_r * mu_y;

    /* Computations stemming from the compressible, Newtonian fluid model */
    const Scalar rhov_t  = rho_t * v + rho * v_t;

    const Scalar rhouv_x = (rho_x * u * v) + (rho * u_x * v) + (rho * u * v_x);
    const Scalar rhovv_y = (rho_y * v * v) + (rho * v_y * v) + (rho * v * v_y);
    const Scalar rhovw_z = (rho_z * v * w) + (rho * v_z * w) + (rho * v * w_z);

    const Scalar tauyy_y = mu_y * (v_y  + v_y )
                         + lambda_y * (u_x  + v_y  + w_z )
                         + mu   * (v_yy + v_yy)
                         + lambda   * (u_xy + v_yy + w_yz);

    const Scalar tauxy_x = mu_x * (u_y + v_x) + mu * (u_xy + v_xx);
    const Scalar tauyz_z = mu_z * (v_z + w_y) + mu * (v_zz + w_yz);

    const Scalar Q_rhov = rhov_t + rhouv_x + rhovv_y + rhovw_z
                        + p_y - tauxy_x - tauyy_y - tauyz_z;

    return Q_rhov;
}

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::Q_rhow(
        T1 x, T2 y, T3 z, T4 t) const
{
    /* Computations retrieving primitive solution details */
    const Scalar rho   = this->rho   (x, y, z, t);  // shadow
    const Scalar rho_t = this->rho._t(x, y, z, t);
    const Scalar rho_x = this->rho._x(x, y, z, t);
    const Scalar rho_y = this->rho._y(x, y, z, t);
    const Scalar rho_z = this->rho._z(x, y, z, t);

    const Scalar u    = this->u    (x, y, z, t);    // shadow
    const Scalar u_x  = this->u._x (x, y, z, t);
    const Scalar u_xz = this->u._xz(x, y, z, t);
    const Scalar u_z  = this->u._z (x, y, z, t);

    const Scalar v    = this->v    (x, y, z, t);    // shadow
    const Scalar v_y  = this->v._y (x, y, z, t);
    const Scalar v_yz = this->v._yz(x, y, z, t);
    const Scalar v_z  = this->v._z (x, y, z, t);

    const Scalar w    = this->w    (x, y, z, t);    // shadow
    const Scalar w_t  = this->w._t (x, y, z, t);
    const Scalar w_x  = this->w._x (x, y, z, t);
    const Scalar w_xx = this->w._xx(x, y, z, t);
    const Scalar w_y  = this->w._y (x, y, z, t);
    const Scalar w_yy = this->w._yy(x, y, z, t);
    const Scalar w_z  = this->w._z (x, y, z, t);
    const Scalar w_zz = this->w._zz(x, y, z, t);

    const Scalar T    = this->T    (x, y, z, t);    // shadow
    const Scalar T_x  = this->T._x (x, y, z, t);
    const Scalar T_y  = this->T._y (x, y, z, t);
    const Scalar T_z  = this->T._z (x, y, z, t);

    /* Computations stemming from the constitutive relationships */
    const Scalar p_z      = rho_z * R * T + rho * R * T_z;
    const Scalar mu       = mu_r * ::std::pow(T / T_r, beta);
    const Scalar mu_x     = beta * mu_r / T_r * ::std::pow(T / T_r, beta - 1)
                          * T_x;
    const Scalar mu_y     = beta * mu_r / T_r * ::std::pow(T / T_r, beta - 1)
                          * T_y;
    const Scalar mu_z     = beta * mu_r / T_r * ::std::pow(T / T_r, beta - 1)
                          * T_z;
    const Scalar lambda   = lambda_r / mu_r * mu;
    const Scalar lambda_z = lambda_r / mu_r * mu_z;

    /* Computations stemming from the compressible, Newtonian fluid model */
    const Scalar rhow_t  = rho_t * w + rho * w_t;

    const Scalar rhouw_x = (rho_x * u * w) + (rho * u_x * w) + (rho * u * w_x);
    const Scalar rhovw_y = (rho_y * v * w) + (rho * v_y * w) + (rho * v * w_y);
    const Scalar rhoww_z = (rho_z * w * w) + (rho * w_z * w) + (rho * w * w_z);

    const Scalar tauzz_z = mu_z * (w_z  + w_z )
                         + lambda_z * (u_x  + v_y  + w_z )
                         + mu   * (w_zz + w_zz)
                         + lambda   * (u_xz + v_yz + w_zz);

    const Scalar tauxz_x = mu_x * (u_z + w_x) + mu * (u_xz + w_xx);
    const Scalar tauyz_y = mu_y * (v_z + w_y) + mu * (v_yz + w_yy);

    const Scalar Q_rhow = rhow_t + rhouw_x + rhovw_y + rhoww_z
                        + p_z - tauxz_x - tauyz_y - tauzz_z;

    return Q_rhow;
}

NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
template <typename T1, typename T2, typename T3, typename T4>
Scalar manufactured_solution<NSCTPL_MANUFACTURED_SOLUTION_TPNAMES>::Q_rhoe(
        T1 x, T2 y, T3 z, T4 t) const
{
    /* Computations retrieving primitive solution details */
    const Scalar rho   = this->rho   (x, y, z, t);  // shadow
    const Scalar rho_t = this->rho._t(x, y, z, t);
    const Scalar rho_x = this->rho._x(x, y, z, t);
    const Scalar rho_y = this->rho._y(x, y, z, t);
    const Scalar rho_z = this->rho._z(x, y, z, t);

    const Scalar u    = this->u    (x, y, z, t);    // shadow
    const Scalar u_t  = this->u._t (x, y, z, t);
    const Scalar u_x  = this->u._x (x, y, z, t);
    const Scalar u_xx = this->u._xx(x, y, z, t);
    const Scalar u_xy = this->u._xy(x, y, z, t);
    const Scalar u_xz = this->u._xz(x, y, z, t);
    const Scalar u_y  = this->u._y (x, y, z, t);
    const Scalar u_yy = this->u._yy(x, y, z, t);
    const Scalar u_z  = this->u._z (x, y, z, t);
    const Scalar u_zz = this->u._zz(x, y, z, t);

    const Scalar v    = this->v    (x, y, z, t);    // shadow
    const Scalar v_t  = this->v._t (x, y, z, t);
    const Scalar v_x  = this->v._x (x, y, z, t);
    const Scalar v_xx = this->v._xx(x, y, z, t);
    const Scalar v_xy = this->v._xy(x, y, z, t);
    const Scalar v_y  = this->v._y (x, y, z, t);
    const Scalar v_yy = this->v._yy(x, y, z, t);
    const Scalar v_yz = this->v._yz(x, y, z, t);
    const Scalar v_z  = this->v._z (x, y, z, t);
    const Scalar v_zz = this->v._zz(x, y, z, t);

    const Scalar w    = this->w    (x, y, z, t);    // shadow
    const Scalar w_t  = this->w._t (x, y, z, t);
    const Scalar w_x  = this->w._x (x, y, z, t);
    const Scalar w_xx = this->w._xx(x, y, z, t);
    const Scalar w_xz = this->w._xz(x, y, z, t);
    const Scalar w_y  = this->w._y (x, y, z, t);
    const Scalar w_yy = this->w._yy(x, y, z, t);
    const Scalar w_yz = this->w._yz(x, y, z, t);
    const Scalar w_z  = this->w._z (x, y, z, t);
    const Scalar w_zz = this->w._zz(x, y, z, t);

    const Scalar T    = this->T    (x, y, z, t);    // shadow
    const Scalar T_t  = this->T._t (x, y, z, t);
    const Scalar T_x  = this->T._x (x, y, z, t);
    const Scalar T_xx = this->T._xx(x, y, z, t);
    const Scalar T_y  = this->T._y (x, y, z, t);
    const Scalar T_yy = this->T._yy(x, y, z, t);
    const Scalar T_z  = this->T._z (x, y, z, t);
    const Scalar T_zz = this->T._zz(x, y, z, t);

    /* Computations stemming from the constitutive relationships */
    const Scalar e        = R * T   / (gamma - 1) + (u*u   + v*v   + w*w  ) / 2;
    const Scalar e_x      = R * T_x / (gamma - 1) + (u*u_x + v*v_x + w*w_x);
    const Scalar e_y      = R * T_y / (gamma - 1) + (u*u_y + v*v_y + w*w_y);
    const Scalar e_z      = R * T_z / (gamma - 1) + (u*u_z + v*v_z + w*w_z);
    const Scalar e_t      = R * T_t / (gamma - 1) + (u*u_t + v*v_t + w*w_t);
    const Scalar p        = rho * R * T;
    const Scalar p_x      = rho_x * R * T + rho * R * T_x;
    const Scalar p_y      = rho_y * R * T + rho * R * T_y;
    const Scalar p_z      = rho_z * R * T + rho * R * T_z;
    const Scalar mu       = mu_r * ::std::pow(T / T_r, beta);
    const Scalar mu_x     = beta * mu_r / T_r * ::std::pow(T / T_r, beta - 1)
                          * T_x;
    const Scalar mu_y     = beta * mu_r / T_r * ::std::pow(T / T_r, beta - 1)
                          * T_y;
    const Scalar mu_z     = beta * mu_r / T_r * ::std::pow(T / T_r, beta - 1)
                          * T_z;
    const Scalar lambda   = lambda_r / mu_r * mu;
    const Scalar lambda_x = lambda_r / mu_r * mu_x;
    const Scalar lambda_y = lambda_r / mu_r * mu_y;
    const Scalar lambda_z = lambda_r / mu_r * mu_z;
    const Scalar qx_x     = - kappa_r / mu_r * (mu_x * T_x + mu * T_xx);
    const Scalar qy_y     = - kappa_r / mu_r * (mu_y * T_y + mu * T_yy);
    const Scalar qz_z     = - kappa_r / mu_r * (mu_z * T_z + mu * T_zz);

    /* Computations stemming from the compressible, Newtonian fluid model */
    const Scalar rhoe_t  = rho_t * e + rho * e_t;

    const Scalar rhoue_x = (rho_x * u * e) + (rho * u_x * e) + (rho * u * e_x);
    const Scalar rhove_y = (rho_y * v * e) + (rho * v_y * e) + (rho * v * e_y);
    const Scalar rhowe_z = (rho_z * w * e) + (rho * w_z * e) + (rho * w * e_z);

    const Scalar tauxx = mu * (u_x + u_x) + lambda * (u_x + v_y + w_z);
    const Scalar tauyy = mu * (v_y + v_y) + lambda * (u_x + v_y + w_z);
    const Scalar tauzz = mu * (w_z + w_z) + lambda * (u_x + v_y + w_z);
    const Scalar tauxy = mu * (u_y + v_x);
    const Scalar tauxz = mu * (u_z + w_x);
    const Scalar tauyz = mu * (v_z + w_y);

    const Scalar tauxx_x = mu_x * (u_x  + u_x )
                         + lambda_x * (u_x  + v_y  + w_z )
                         + mu   * (u_xx + u_xx)
                         + lambda   * (u_xx + v_xy + w_xz);
    const Scalar tauyy_y = mu_y * (v_y  + v_y )
                         + lambda_y * (u_x  + v_y  + w_z )
                         + mu   * (v_yy + v_yy)
                         + lambda   * (u_xy + v_yy + w_yz);
    const Scalar tauzz_z = mu_z * (w_z  + w_z )
                         + lambda_z * (u_x  + v_y  + w_z )
                         + mu   * (w_zz + w_zz)
                         + lambda   * (u_xz + v_yz + w_zz);

    const Scalar tauxy_x = mu_x * (u_y + v_x) + mu * (u_xy + v_xx);
    const Scalar tauxy_y = mu_y * (u_y + v_x) + mu * (u_yy + v_xy);
    const Scalar tauxz_x = mu_x * (u_z + w_x) + mu * (u_xz + w_xx);
    const Scalar tauxz_z = mu_z * (u_z + w_x) + mu * (u_zz + w_xz);
    const Scalar tauyz_y = mu_y * (v_z + w_y) + mu * (v_yz + w_yy);
    const Scalar tauyz_z = mu_z * (v_z + w_y) + mu * (v_zz + w_yz);

    const Scalar pu_x = p_x * u + p * u_x;
    const Scalar pv_y = p_y * v + p * v_y;
    const Scalar pw_z = p_z * w + p * w_z;
    const Scalar utauxx_x = u_x * tauxx + u * tauxx_x;
    const Scalar vtauxy_x = v_x * tauxy + v * tauxy_x;
    const Scalar wtauxz_x = w_x * tauxz + w * tauxz_x;
    const Scalar utauxy_y = u_y * tauxy + u * tauxy_y;
    const Scalar vtauyy_y = v_y * tauyy + v * tauyy_y;
    const Scalar wtauyz_y = w_y * tauyz + w * tauyz_y;
    const Scalar utauxz_z = u_z * tauxz + u * tauxz_z;
    const Scalar vtauyz_z = v_z * tauyz + v * tauyz_z;
    const Scalar wtauzz_z = w_z * tauzz + w * tauzz_z;

    const Scalar Q_rhoe = rhoe_t + rhoue_x + rhove_y + rhowe_z
                        + pu_x + pv_y + pw_z
                        + qx_x + qy_y + qz_z
                        - utauxx_x - vtauxy_x - wtauxz_x
                        - utauxy_y - vtauyy_y - wtauyz_y
                        - utauxz_z - vtauyz_z - wtauzz_z;

    return Q_rhoe;
}

#undef NSCTPL_MANUFACTURED_SOLUTION_TPDECLARATION
#undef NSCTPL_MANUFACTURED_SOLUTION_TPNAMES

// ---------------------------------------------------------------------------
// Utilities: zero
// ---------------------------------------------------------------------------

namespace { // anonymous

template <typename Scalar> void zero_helper(const std::string&, Scalar& value)
{
    value = 0;
}

}

template <class T> void zero(T& t)
{
    void (*z)(const std::string&, typename T::scalar_type&) = &zero_helper;
    t.foreach_parameter(z);
}

// ---------------------------------------------------------------------------
// Utilities: isothermal channel recommended parameters
// ---------------------------------------------------------------------------

template <class T> void isothermal_channel(T& t)
{
    typedef typename T::scalar_type scalar;
    const scalar pi = 4 * atan(scalar(1));

    zero(t);   // All parameters not set here must be zeros

    t.gamma    = scalar(14)/scalar(10);
    t.R        = scalar(287);
    t.beta     = scalar(2)/scalar(3);
    t.mu_r     = scalar(1852)/scalar(100000000);
    t.T_r      = scalar(300);
    t.kappa_r  = (t.gamma*t.R*t.mu_r) / ((t.gamma - 1)*scalar(7)/scalar(10));
    t.lambda_r = - scalar(2)/scalar(3)*t.mu_r;
    t.Lx       = scalar(4)*pi;
    t.Ly       = scalar(2);
    t.Lz       = scalar(4)*pi/scalar(3);

    t.rho.a_0  = scalar(1);
    t.rho.a_xy = scalar(1) / scalar(11);
    t.rho.b_xy = scalar(3);
    t.rho.d_xy = scalar(3);
    t.rho.f_xy = scalar(3);
    t.rho.g_xy = pi / scalar(4);
    t.rho.a_y  = scalar(1) / scalar(7);
    t.rho.b_y  = scalar(1) / scalar(2);
    t.rho.f_y  = scalar(1);
    t.rho.g_y  = pi / scalar(4) - scalar(1) / scalar(20);
    t.rho.a_yz = scalar(1) / scalar(31);
    t.rho.b_yz = scalar(2);
    t.rho.d_yz = scalar(2);
    t.rho.f_yz = scalar(2);
    t.rho.g_yz = pi / scalar(4) + scalar(1) / scalar(20);

    t.u.a_xy = scalar(53) / scalar(37);
    t.u.b_xy = scalar(3);
    t.u.c_xy = - pi / scalar(2);
    t.u.d_xy = scalar(3);
    t.u.e_xy = - pi / scalar(2);
    t.u.f_xy = scalar(3);
    t.u.g_xy = pi / scalar(4);
    t.u.a_y  = scalar(53);
    t.u.b_y  = scalar(1) / scalar(2);
    t.u.c_y  = - pi / scalar(2);
    t.u.f_y  = scalar(1);
    t.u.g_y  = pi / scalar(4) - scalar(1) / scalar(20);
    t.u.a_yz = scalar(53) / scalar(41);
    t.u.b_yz = scalar(2);
    t.u.c_yz = - pi / scalar(2);
    t.u.d_yz = scalar(2);
    t.u.e_yz = - pi / scalar(2);
    t.u.f_yz = scalar(2);
    t.u.g_yz = pi / scalar(4) + scalar(1) / scalar(20);

    t.v.a_xy = scalar(3);
    t.v.b_xy = scalar(3);
    t.v.c_xy = - pi / scalar(2);
    t.v.d_xy = scalar(3);
    t.v.e_xy = - pi / scalar(2);
    t.v.f_xy = scalar(3);
    t.v.g_xy = pi / scalar(4);
    t.v.a_y  = scalar(2);
    t.v.b_y  = scalar(1) / scalar(2);
    t.v.c_y  = - pi / scalar(2);
    t.v.f_y  = scalar(1);
    t.v.g_y  = pi / scalar(4) - scalar(1) / scalar(20);
    t.v.a_yz = scalar(5);
    t.v.b_yz = scalar(2);
    t.v.c_yz = -pi / scalar(2);
    t.v.d_yz = scalar(2);
    t.v.e_yz = -pi / scalar(2);
    t.v.f_yz = scalar(2);
    t.v.g_yz = pi / scalar(4) + scalar(1) / scalar(20);

    t.w.a_xy = scalar(11);
    t.w.b_xy = scalar(3);
    t.w.c_xy = -pi / scalar(2);
    t.w.d_xy = scalar(3);
    t.w.e_xy = -pi / scalar(2);
    t.w.f_xy = scalar(3);
    t.w.g_xy = pi / scalar(4);
    t.w.a_y  = scalar(7);
    t.w.b_y  = scalar(1) / scalar(2);
    t.w.c_y  = - pi / scalar(2);
    t.w.f_y  = scalar(1);
    t.w.g_y  = pi / scalar(4) - scalar(1) / scalar(20);
    t.w.a_yz = scalar(13);
    t.w.b_yz = scalar(2);
    t.w.c_yz = - pi / scalar(2);
    t.w.d_yz = scalar(2);
    t.w.e_yz = - pi / scalar(2);
    t.w.f_yz = scalar(2);
    t.w.g_yz = pi / scalar(4) + scalar(1) / scalar(20);

    t.T.a_0  = scalar(300);
    t.T.a_xy = scalar(300) / scalar(17);
    t.T.b_xy = scalar(3);
    t.T.c_xy = - pi / scalar(2);
    t.T.d_xy = scalar(3);
    t.T.e_xy = - pi / scalar(2);
    t.T.f_xy = scalar(3);
    t.T.g_xy = pi / scalar(4);
    t.T.a_y  = scalar(300) / scalar(13);
    t.T.b_y  = scalar(1) / scalar(2);
    t.T.c_y  = - pi / scalar(2);
    t.T.f_y  = scalar(1);
    t.T.g_y  = pi / scalar(4) - scalar(1) / scalar(20);
    t.T.a_yz = scalar(300) / scalar(37);
    t.T.b_yz = scalar(2);
    t.T.c_yz = - pi / scalar(2);
    t.T.d_yz = scalar(2);
    t.T.e_yz = - pi / scalar(2);
    t.T.f_yz = scalar(2);
    t.T.g_yz = pi / scalar(4) + scalar(1) / scalar(20);
}

// ---------------------------------------------------------------------------
// Utilities: isothermal flat plate recommended parameters
// ---------------------------------------------------------------------------

template <class T> void isothermal_flat_plate(T& t)
{
    typedef typename T::scalar_type scalar;

    // Implemented as a delta relative to the isothermal channel
    isothermal_channel(t);
    t.rho.b_y = scalar(1) / scalar(4);
    t.u.b_y   = scalar(1) / scalar(4);
    t.v.b_y   = scalar(1) / scalar(4);
    t.w.b_y   = scalar(1) / scalar(4);
    t.T.b_y   = scalar(1) / scalar(4);
}

} // end namespace nsctpl

#endif /* NSCTPL_HPP */
