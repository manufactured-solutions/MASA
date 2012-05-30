#ifndef __dualnumberarray_h__
#define __dualnumberarray_h__


#include "dualnumber.h"
#include "numberarray.h"


template <std::size_t size, typename T>
struct DerivativeType<NumberArray<size, T> >
{
  typedef NumberArray<size, typename DerivativeType<T>::type> type;
};


template <std::size_t size, typename T>
struct DerivativesOf<NumberArray<size, T> >
{
  static
  typename DerivativeType<NumberArray<size, T> >::type
  derivative(const NumberArray<size, T>& a, unsigned int derivativeindex)
  {
    typename DerivativeType<NumberArray<size, T> >::type returnval;
    for (unsigned int i=0; i != size; ++i)
      returnval[i] = DerivativesOf<T>::derivative(a[i], derivativeindex);
  
    return returnval;
  }
};


template <std::size_t size, typename T, unsigned int derivativeindex>
struct DerivativeOf<NumberArray<size, T>, derivativeindex>
{
  static
  typename DerivativeType<NumberArray<size, T> >::type
  derivative(const NumberArray<size, T>& a)
  {
    typename DerivativeType<NumberArray<size, T> >::type returnval;
    for (unsigned int i=0; i != size; ++i)
      returnval[i] = DerivativeOf<T,derivativeindex>::derivative(a[i]);
  
    return returnval;
  }
};


// For a vector of values a[i] each of which has a defined gradient,
// the divergence is the sum of derivative_wrt_xi(a[i])

// For a tensor of values, we take the divergence with respect to the
// first index.
template <std::size_t size, typename T>
inline
typename DerivativeType<T>::type
divergence(const NumberArray<size, T>& a)
{
  typename DerivativeType<T>::type returnval = 0;

  for (unsigned int i=0; i != size; ++i)
    returnval += DerivativesOf<T>::derivative(a[i], i);

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
    returnval[i] = gradient(a[i]);

  return returnval;
}

// DualNumber is subordinate to NumberArray

#define DualNumberArray_comparisons(templatename) \
template<typename T, typename D, std::size_t size, typename T2, bool reverseorder> \
struct templatename<NumberArray<size, T2>, DualNumber<T, D>, reverseorder> { \
  typedef NumberArray<size, typename Symmetric##templatename<DualNumber<T, D>, T2, reverseorder>::supertype> supertype; \
}

DualNumberArray_comparisons(CompareTypes);
DualNumberArray_comparisons(PlusType);
DualNumberArray_comparisons(MinusType);
DualNumberArray_comparisons(MultipliesType);
DualNumberArray_comparisons(DividesType);

#endif // __dualnumberarray_h__
