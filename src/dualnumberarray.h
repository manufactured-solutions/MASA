// -*-c++-*-
//
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
//
// dualnumberarray.h: 
//
// $Id: $
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef __dualnumberarray_h__
#define __dualnumberarray_h__

template <typename T>
struct DivergenceType
{
  typedef typename T::derivatives_type::value_type divergence_type;
};

template <>
template <std::size_t size, typename T>
struct DivergenceType<NumberArray<size, T> >
{
  typedef typename T::derivatives_type::value_type divergence_type;
};

// For a vector of values a[i] each of which has a defined gradient,
// the divergence is the sum of derivative_wrt_xi(a[i])
template <std::size_t size, typename T>
inline
typename DivergenceType<T>::divergence_type
divergence(const NumberArray<size, T>& a)
{
  typename DivergenceType<T>::divergence_type returnval = 0;

  for (unsigned int i=0; i != size; ++i)
    returnval += a[i].derivatives()[i];

  return returnval;
}

// For a tensor of values, we take the divergence with respect to the
// last index.  This way we only have to have NumberArray<DualNumber>
// data structure support, not every other combination.

template <std::size_t size, typename T>
inline
NumberArray<size, typename DivergenceType<T>::divergence_type>
divergence(const NumberArray<size, NumberArray<size, T> >& a)
{
  NumberArray<size, typename DivergenceType<T>::divergence_type> returnval;

  for (unsigned int i=0; i != size; ++i)
    returnval[i] = divergence(a[i]);

  return returnval;
}

// For a vector of values, the gradient is going to be a tensor
template <std::size_t size, typename T>
inline
NumberArray<size, typename T::derivatives_type>
gradient(const NumberArray<size, T>& a)
{
  NumberArray<size, typename T::derivatives_type> returnval;

  for (unsigned int i=0; i != size; ++i)
    returnval[i] = a[i].derivatives();

  return returnval;
}

// The compiler can't choose between NumberArray-vs-T and
// T-vs-DualNumber CompareTypes specializations, so we have to choose
// for it.

template<std::size_t size, typename T, typename T2, typename D>
struct CompareTypes<NumberArray<size, T2>, DualNumber<T, D> > {
  typedef NumberArray<size, typename CompareTypes<T2, DualNumber<T, D> >::supertype> supertype;
};

template<std::size_t size, typename T, typename T2, typename D>
struct CompareTypes<DualNumber<T, D>, NumberArray<size, T2> > {
  typedef NumberArray<size, typename CompareTypes<T2, DualNumber<T, D> >::supertype> supertype;
};



#endif // __dualnumberarray_h__
