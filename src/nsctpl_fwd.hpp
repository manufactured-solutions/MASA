//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

// Using this file will require #include-ing <cmath>, <limits>, <sstream>, and
// <string>.  Those headers are not #include-d here to allow #include-ing this
// file inside any namespace, as is the consistent use of ::std instead of
// merely std.

#ifndef NSCTPL_FWD_HPP
#define NSCTPL_FWD_HPP

namespace nsctpl {

/**
 * Class providing operations for evaluating an analytical solution at a
 * given location and time.  The solution is of the form
 * \verbatim
 *       a_0                                         *cos(f_0 *t + g_0 )
 *     + a_x  * cos(b_x *x + c_x )                   *cos(f_x *t + g_x )
 *     + a_xy * cos(b_xy*x + c_xy)*cos(d_xy*y + e_xy)*cos(f_xy*t + g_xy)
 *     + a_xz * cos(b_xz*x + c_xz)*cos(d_xz*z + e_xz)*cos(f_xz*t + g_xz)
 *     + a_y  * cos(b_y *y + c_y )                   *cos(f_y *t + g_y )
 *     + a_yz * cos(b_yz*y + c_yz)*cos(d_yz*z + e_yz)*cos(f_yz*t + g_yz)
 *     + a_z  * cos(b_z *z + c_z )                   *cos(f_z *t + g_z )
 * \endverbatim
 * The member method names are non-traditional but permitted by the language
 * standards and well-aligned with mathematical notation.  Location and
 * time arguments are templated to allow extended precision intermediate
 * computations (followed by truncation) when possible.
 *
 * The foreach_parameter() member allows invoking an operation on each solution
 * parameter to aid parameter registration, initialization, or output.
 * Parameters may have an infix name added for use with foreach_parameter().
 * For example, setting the name 'phi' will cause foreach_parameter() to report
 * names like 'a_phix'.
 */
template <typename Scalar>
class primitive {

public:

    typedef Scalar scalar_type;  //!< Scalar type employed

     // X macro per http://drdobbs.com/blogs/cpp/228700289 which will
     // apply a macro on all solution parameter prefix/suffix pairs.
#define NSCTPL_FOR_EACH_PRIMITIVE_PARAMETER(apply)                             \
    apply(a_,0)                                                                \
    apply(a_,x) apply(a_,xy) apply(a_,xz) apply(a_,y) apply(a_,yz) apply(a_,z) \
    apply(b_,x) apply(b_,xy) apply(b_,xz) apply(b_,y) apply(b_,yz) apply(b_,z) \
    apply(c_,x) apply(c_,xy) apply(c_,xz) apply(c_,y) apply(c_,yz) apply(c_,z) \
                apply(d_,xy) apply(d_,xz)             apply(d_,yz)             \
                apply(e_,xy) apply(e_,xz)             apply(e_,yz)             \
    apply(f_,0)                                                                \
    apply(f_,x) apply(f_,xy) apply(f_,xz) apply(f_,y) apply(f_,yz) apply(f_,z) \
    apply(g_,0)                                                                \
    apply(g_,x) apply(g_,xy) apply(g_,xz) apply(g_,y) apply(g_,yz) apply(g_,z)

    // Declare all solution parameters as members, e.g. 'a_xy'
#define NSCTPL_APPLY(pre,suf) Scalar pre##suf;
    NSCTPL_FOR_EACH_PRIMITIVE_PARAMETER(NSCTPL_APPLY)
#undef NSCTPL_APPLY

    //! The name used infix in foreach_parameter' names.  For example,
    //! name = 'phi' implies parameters names like 'a_phix'.
    ::std::string name;

    const Scalar* Lx;           //!< Domain extent in x direction
    const Scalar* Ly;           //!< Domain extent in y direction
    const Scalar* Lz;           //!< Domain extent in z direction

    //! Construct an instance using \c name in the reported parameter names.
    //! The domain sizes are referenced from some external Lx, Ly, and Lz.
    //! All parameters set to zero at construction time.
    explicit primitive(const ::std::string &name,
                       const Scalar& Lx,
                       const Scalar& Ly,
                       const Scalar& Lz)
#define NSCTPL_APPLY(pre,suf) pre##suf(0),
        : NSCTPL_FOR_EACH_PRIMITIVE_PARAMETER(NSCTPL_APPLY)  // has final comma
#undef NSCTPL_APPLY
          name(name), Lx(&Lx), Ly(&Ly), Lz(&Lz)
    {}

    //! Construct an instance using \c name in the reported parameter names.
    //! Domain sizes Lx, Ly, and Lz \b must be set prior to invoking any
    //! member methods.  All parameters set to zero at construction time.
    explicit primitive(const ::std::string &name = "")
#define NSCTPL_APPLY(pre,suf) pre##suf(0),
        : NSCTPL_FOR_EACH_PRIMITIVE_PARAMETER(NSCTPL_APPLY)  // has final comma
#undef NSCTPL_APPLY
          name(name), Lx(NULL), Ly(NULL), Lz(NULL)
    {}

#define NSCTPL_APPLY_STRINGIFY(s) #s
#define NSCTPL_APPLY(pre,suf)              \
        os.clear(); os.str("");            \
        os << NSCTPL_APPLY_STRINGIFY(pre)  \
           << this->name                   \
           << NSCTPL_APPLY_STRINGIFY(suf); \
        f(os.str(), this->pre##suf);

    //! Invoke the binary function f on each parameter name and its value.
    template <typename BinaryFunction>
    void foreach_parameter(BinaryFunction& f) const {
        ::std::ostringstream os;
        NSCTPL_FOR_EACH_PRIMITIVE_PARAMETER(NSCTPL_APPLY)
    }

    //! Invoke the binary function f on each parameter name and its value.
    template <typename BinaryFunction>
    void foreach_parameter(BinaryFunction& f) {
        ::std::ostringstream os;
        NSCTPL_FOR_EACH_PRIMITIVE_PARAMETER(NSCTPL_APPLY)
    }

#undef NSCTPL_APPLY
#undef NSCTPL_APPLY_STRINGIFY

    //! Evaluate the solution
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar operator()(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to time
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _t(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to \c x
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _x(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's second derivative with respect to \c x
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _xx(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to \c x and \c y
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _xy(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to \c x and \c z
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _xz(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to \c y
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _y(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's second derivative with respect to \c y
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _yy(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to \c y and \c z
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _yz(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to \c z
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _z(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's second derivative with respect to \c z
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _zz(T1 x, T2 y, T3 z, T4 t) const;

private:

    static const Scalar twopi;  //!< \f$2\pi\f$ to \c Scalar precision

}; // end class

// It is handy to template manufactured_solution (just below) on the
// primitive functional forms chosen for rho, u, v, w, and T.  However
// SWIG, even version 2.0.3, has trouble with template template parameters:
// https://sourceforge.net/tracker/?func=detail&atid=101645&aid=1861407&group_id=1645
// Therefore only employ the template template parameter when SWIG is not used.

/**
 * Template that, given a primitive function for \c rho, \c u, \c v, \c w, and
 * \c T along with a floating point type, evaluates a manufactured solution and
 * the required forcing for the transient, compressible Navier--Stokes
 * equations with a power law viscosity.  Location and time arguments are
 * templated to allow extended precision intermediate computations (followed by
 * truncation) when possible.
 */
#ifndef SWIG
template <typename Scalar,
          int IndexBase = 0,
          template <typename> class Primitive = primitive>
#else
template <typename Scalar, int IndexBase = 0>
#endif
class manufactured_solution {

public:

    typedef Scalar scalar_type;                //!< Scalar type employed
    static const int index_base;               //!< Indexing base for gradients
#ifndef SWIG
    typedef Primitive<Scalar> primitive_type;  //!< Primitive solution type
#else
    typedef primitive<Scalar> primitive_type;  //!< Primitive solution type
#endif

    //! Scenario parameters
    //!@{
    Scalar gamma;            //!< Constant ratio of specific heats
    Scalar R;                //!< Gas constant
    Scalar beta;             //!< Viscosity power law exponent
    Scalar mu_r;             //!< Reference dynamic viscosity
    Scalar T_r;              //!< Reference temperature
    Scalar kappa_r;          //!< Reference thermal conductivity
    Scalar lambda_r;         //!< Reference bulk viscosity
    //!@}

    //! Analytic solutions (which contain additional parameters)
    //!@{
    primitive_type rho;      //!< Analytic solution for rho
    primitive_type u;        //!< Analytic solution for u
    primitive_type v;        //!< Analytic solution for v
    primitive_type w;        //!< Analytic solution for w
    primitive_type T;        //!< Analytic solution for T
    //!@}

    //! Domain extents
    //!@{
    Scalar Lx;               //!< Domain extent in x direction
    Scalar Ly;               //!< Domain extent in y direction
    Scalar Lz;               //!< Domain extent in z direction
    //!@}

    //! Default constructor
    manufactured_solution()
        : gamma(0), R(0), beta(0), mu_r(0), T_r(0), kappa_r(0), lambda_r(0),
          rho("rho", Lx, Ly, Lz),
          u  ("u"  , Lx, Ly, Lz),
          v  ("v"  , Lx, Ly, Lz),
          w  ("w"  , Lx, Ly, Lz),
          T  ("T"  , Lx, Ly, Lz),
          Lx(1), Ly(1), Lz(1)
    {
    }

    //! Invoke the binary function f on each parameter name and its constant value.
    template <typename BinaryFunction>
    void foreach_parameter(BinaryFunction& f) const {
        f(::std::string("gamma"),    gamma   );
        f(::std::string("R"),        R       );
        f(::std::string("beta"),     beta    );
        f(::std::string("mu_r"),     mu_r    );
        f(::std::string("T_r"),      T_r     );
        f(::std::string("kappa_r"),  kappa_r );
        f(::std::string("lambda_r"), lambda_r);
        rho.foreach_parameter(f);
        u.foreach_parameter(f);
        v.foreach_parameter(f);
        w.foreach_parameter(f);
        T.foreach_parameter(f);
        f(::std::string("Lx"), Lx);
        f(::std::string("Ly"), Ly);
        f(::std::string("Lz"), Lz);
    }

    //! Invoke the binary function f on each parameter name and its mutable value.
    template <typename BinaryFunction>
    void foreach_parameter(BinaryFunction& f) {
        f(::std::string("gamma"),    gamma   );
        f(::std::string("R"),        R       );
        f(::std::string("beta"),     beta    );
        f(::std::string("mu_r"),     mu_r    );
        f(::std::string("T_r"),      T_r     );
        f(::std::string("kappa_r"),  kappa_r );
        f(::std::string("lambda_r"), lambda_r);
        rho.foreach_parameter(f);
        u.foreach_parameter(f);
        v.foreach_parameter(f);
        w.foreach_parameter(f);
        T.foreach_parameter(f);
        f(::std::string("Lx"), Lx);
        f(::std::string("Ly"), Ly);
        f(::std::string("Lz"), Lz);
    }

    // Analytically determined quantities
    // Note that Primitive members can be used directly.
    // For example, T(x,y,z,t) or T._xx(x,y,z,t)
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_rho (T1 x, T2 y, T3 z, T4 t, int index) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_u   (T1 x, T2 y, T3 z, T4 t, int index) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_v   (T1 x, T2 y, T3 z, T4 t, int index) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_w   (T1 x, T2 y, T3 z, T4 t, int index) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_T   (T1 x, T2 y, T3 z, T4 t, int index) const;

    // Quantities built from the analytical solutions
    // TODO Build up q_u, q_v, q_w, q_e, q_T, q_p

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar e(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar p(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar mu(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar rhou(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar rhov(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar rhow(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar rhoe(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_e(T1 x, T2 y, T3 z, T4 t, int index) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_p(T1 x, T2 y, T3 z, T4 t, int index) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_mu(T1 x, T2 y, T3 z, T4 t, int index) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar Q_rho(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar Q_rhou(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar Q_rhov(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar Q_rhow(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar Q_rhoe(T1 x, T2 y, T3 z, T4 t) const;

}; // end class


// Utilities take generic templated types to allow reuse by
// any variants satisfying the necessary APIs.

/** Zero all of an instance's parameters using foreach_parameter */
template <class T> void zero(T& t);

/** Set recommended isothermal channel problem parameters per write up */
template <class T> void isothermal_channel(T& t);

/** Set recommended isothermal channel problem parameters per write up */
template <class T> void isothermal_flat_plate(T& t);

} // end namespace nsctpl

#endif /* NSCTPL_FWD_HPP */
