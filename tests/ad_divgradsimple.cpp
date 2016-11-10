// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012,2013 The PECOS Development Team
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
// $Author: clark $
// $Id:
//
// ad_divgradsimple.cpp: Check that the divergence of the gradient is correct
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <iostream>
#include <config.h>

#ifdef HAVE_METAPHYSICL

#include <tests.h>
#include "ad_masa.h"

using namespace MASA;

// Declarations
template <typename Scalar>
Scalar divgradsimple_x(Scalar x1, Scalar y1, Scalar z1);
template <typename Scalar>
Scalar divgradsimple_y(Scalar x1, Scalar y1, Scalar z1);
template <typename Scalar>
Scalar divgradsimple_z(Scalar x1, Scalar y1, Scalar z1);

typedef ShadowNumber<double, long double> RawScalar;
const unsigned int NDIM = 3;
typedef DualNumber<RawScalar, NumberVector<NDIM, RawScalar> > FirstDerivType;
typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
typedef SecondDerivType ADType;

template<typename Scalar>
int run_test()
{
  // CCP: I hacked this threshold to be very lenient of errors
  const Scalar thresh = 5000 * std::numeric_limits<Scalar>::epsilon();
  // parameters
  Scalar x, y, z;

  //problem size
  int nx = 200;  // number of points
  int lx=10;     // length
  Scalar dx=Scalar(lx)/Scalar(nx);
  Scalar dy=Scalar(lx)/Scalar(nx);
  Scalar dz=Scalar(lx)/Scalar(nx);

  // errors
  Scalar uerror, verror, werror;

  // evaluate MMS (1D)
  for(int i=0;i<nx;i++)
  {
    x = i*dx;
    for (int j=0;j<nx;j++)
    {
        y = i*dy;
        for (int k=0;k<nx;k++)
        {
            z = i*dz;
            Scalar ad_x_soln = divgradsimple_x<Scalar>(x,y,z);
            Scalar ad_y_soln = divgradsimple_y<Scalar>(x,y,z);
            Scalar ad_z_soln = divgradsimple_z<Scalar>(x,y,z);

            Scalar exact_x_soln = 6.0*x + 4.0;
            Scalar exact_y_soln = 6.0*y + 4.0;
            Scalar exact_z_soln = 6.0*z + 4.0;

	// test the result is roughly zero
	// choose between abs and rel error
#ifdef MASA_STRICT_REGRESSION

            uerror = std::abs(ad_x_soln - exact_x_soln);
            verror = std::abs(ad_y_soln - exact_y_soln);
            werror = std::abs(ad_z_soln - exact_z_soln);

#else

            uerror = std::abs(ad_x_soln - exact_x_soln)/std::abs(exact_x_soln);
            verror = std::abs(ad_y_soln - exact_y_soln)/std::abs(exact_y_soln);
            werror = std::abs(ad_z_soln - exact_z_soln)/std::abs(exact_z_soln);

#endif	
	        threshcheck(uerror,thresh);
	        threshcheck(verror,thresh);
	        threshcheck(werror,thresh);
        }
    }
  }

  // tests passed
  return 0;

} // end run_test

template <typename Scalar>
Scalar divgradsimple_x(Scalar x1, Scalar y1, Scalar z1)
{

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberVector<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  //
  NumberVector<NDIM, ADScalar> U;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  // Assign a simple manufactured solution to U
  U[0]       = x*x*x + y*y   + z*z;
  U[1]       = x*x   + y*y*y + z*z;
  U[2]       = x*x   + y*y   + z*z*z;

 NumberVector<NDIM, Scalar> divgrad = 
    raw_value(divergence(gradient(U)));

  return divgrad[0];
}

template <typename Scalar>
Scalar divgradsimple_y(Scalar x1, Scalar y1, Scalar z1)
{

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberVector<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  //
  NumberVector<NDIM, ADScalar> U;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  // Assign a simple manufactured solution to U
  U[0]       = x*x*x + y*y   + z*z;
  U[1]       = x*x   + y*y*y + z*z;
  U[2]       = x*x   + y*y   + z*z*z;

 NumberVector<NDIM, Scalar> divgrad = 
    raw_value(divergence(gradient(U)));

  return divgrad[1];
}

template <typename Scalar>
Scalar divgradsimple_z(Scalar x1, Scalar y1, Scalar z1)
{

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberVector<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  //
  NumberVector<NDIM, ADScalar> U;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  // Assign a simple manufactured solution to U
  U[0]       = x*x*x + y*y   + z*z;
  U[1]       = x*x   + y*y*y + z*z;
  U[2]       = x*x   + y*y   + z*z*z;

 NumberVector<NDIM, Scalar> divgrad = 
    raw_value(divergence(gradient(U)));

  return divgrad[2]; 
}

int main()
{
  int err=0;

  err += run_test<double>();

  return err;
}

#else // HAVE_METAPHYSICL

int main(void)
{
  return 77; // Autotools code for "skip test"
}

#endif // HAVE_METAPHYSICL
