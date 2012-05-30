
#ifndef __numberarray_h__
#define __numberarray_h__

#include <algorithm>
#include <ostream>

#include "compare_types.h"
#include "raw_type.h"

template <std::size_t size, typename T>
class NumberArray
{
public:
  typedef T value_type;

  template <unsigned int i>
  struct entry_type {
    typedef value_type type;
  };

  template <typename T2>
  struct rebind {
    typedef NumberArray<size, T2> other;
  };

  NumberArray() {}

  NumberArray(const T& val)
    { std::fill(_data, _data+size, val); }

  NumberArray(const T* vals)
    { std::copy(vals, vals+size, _data); }

  template <typename T2>
  NumberArray(NumberArray<size, T2> src)
    { if (size) std::copy(&src[0], &src[0]+size, _data); }

  template <typename T2>
  NumberArray(const T2& val)
    { std::fill(_data, _data+size, T(val)); }

  T& operator[](unsigned int i)
    { return _data[i]; }

  const T& operator[](unsigned int i) const
    { return _data[i]; }

  template <unsigned int i>
  typename entry_type<i>::type& get()
    { return _data[i]; }

  template <unsigned int i>
  const typename entry_type<i>::type& get() const
    { return _data[i]; }

  NumberArray<size,T> operator- () const {
    NumberArray<size,T> returnval;
    for (unsigned int i=0; i != size; ++i) returnval[i] = -_data[i];
    return returnval;
  }

  template <typename T2>
  NumberArray<size,T>& operator+= (const NumberArray<size,T2>& a)
    { for (unsigned int i=0; i != size; ++i) _data[i] += a[i]; return *this; }

  template <typename T2>
  NumberArray<size,T>& operator+= (const T2& a)
    { for (unsigned int i=0; i != size; ++i) _data[i] += a; return *this; }

  template <typename T2>
  NumberArray<size,T>& operator-= (const NumberArray<size,T2>& a)
    { for (unsigned int i=0; i != size; ++i) _data[i] -= a[i]; return *this; }

  template <typename T2>
  NumberArray<size,T>& operator-= (const T2& a)
    { for (unsigned int i=0; i != size; ++i) _data[i] -= a; return *this; }

  template <typename T2>
  NumberArray<size,T>& operator*= (const NumberArray<size,T2>& a)
    { for (unsigned int i=0; i != size; ++i) _data[i] *= a[i]; return *this; }

  template <typename T2>
  NumberArray<size,T>& operator*= (const T2& a)
    { for (unsigned int i=0; i != size; ++i) _data[i] *= a; return *this; }

  template <typename T2>
  NumberArray<size,T>& operator/= (const NumberArray<size,T2>& a)
    { for (unsigned int i=0; i != size; ++i) _data[i] /= a[i]; return *this; }

  template <typename T2>
  NumberArray<size,T>& operator/= (const T2& a)
    { for (unsigned int i=0; i != size; ++i) _data[i] /= a; return *this; }

  template <typename T2>
  typename SymmetricMultipliesType<T,T2>::supertype
  dot (const NumberArray<size,T2>& a)
  {
    typename SymmetricMultipliesType<T,T2>::supertype returnval = 0;
    for (unsigned int i=0; i != size; ++i)
      returnval += _data[i] * a[i];
    return returnval;
  }

  template <typename T2>
  NumberArray<size, NumberArray<size, typename SymmetricMultipliesType<T,T2>::supertype> >
  outerproduct (const NumberArray<size,T2>& a)
  {
    NumberArray<size, NumberArray<size, typename SymmetricMultipliesType<T,T2>::supertype> > returnval;

    for (unsigned int i=0; i != size; ++i)
      for (unsigned int j=0; j != size; ++j)
        returnval[i][j] = _data[i] * a[j];

    return returnval;
  }

  static NumberArray<size, NumberArray<size, T> > identity()
  {
    NumberArray<size, NumberArray<size, T> > returnval(0);
  
    for (unsigned int i=0; i != size; ++i)
      returnval[i][i] = 1;

    return returnval;
  }

private:
  T _data[size];
};



//
// Non-member functions
//

template <unsigned int size,
          unsigned int index1=0, typename Data1=void,
          unsigned int index2=0, typename Data2=void,
          unsigned int index3=0, typename Data3=void,
          unsigned int index4=0, typename Data4=void,
          unsigned int index5=0, typename Data5=void,
          unsigned int index6=0, typename Data6=void,
          unsigned int index7=0, typename Data7=void,
          unsigned int index8=0, typename Data8=void>
struct NumberArrayOf
{
  typedef
  typename CompareTypes<Data1,
    typename CompareTypes<Data2,
      typename CompareTypes<Data3,
        typename CompareTypes<Data4,
          typename CompareTypes<Data5,
            typename CompareTypes<Data6,
              typename CompareTypes<Data7,Data8>::supertype
            >::supertype
          >::supertype
        >::supertype
      >::supertype
    >::supertype
  >::supertype supertype;

  typedef NumberArray<size, supertype> type;
};



template <std::size_t size, unsigned int index, typename T>
struct NumberArrayUnitVector
{
  typedef NumberArray<size, T> type;

  static const type value() {
    type returnval = 0;
    returnval[index] = 1;
    return returnval;
  }
};


template <std::size_t size, typename T>
struct NumberArrayFullVector
{
  typedef NumberArray<size,T> type;

  static const type value() {
    type returnval;
    for (unsigned int i=0; i != size; ++i)
      returnval[i] = 1;
    return returnval;
  }
};



template <std::size_t size, typename T>
inline
NumberArray<size, NumberArray<size, T> >
transpose(const NumberArray<size, NumberArray<size, T> >& a)
{
  NumberArray<size, NumberArray<size, T> > returnval = a;

  for (unsigned int i=0; i != size; ++i)
    for (unsigned int j=i+1; j != size; ++j)
      std::swap(returnval[i][j], returnval[j][i]);

  return returnval;
}



#define NumberArray_op_ab(opname, atype, btype, newtype) \
template <std::size_t size, typename T, typename T2> \
inline \
typename newtype::supertype \
operator opname (const atype& a, const btype& b) \
{ \
  typedef typename newtype::supertype TS; \
  TS returnval(a); \
  returnval opname##= b; \
  return returnval; \
}

#define NumberArray_op(opname, typecomparison) \
NumberArray_op_ab(opname, NumberArray<size MacroComma T>, NumberArray<size MacroComma T2>, \
                  typecomparison##Type<NumberArray<size MacroComma T> MacroComma NumberArray<size MacroComma T2> >) \
NumberArray_op_ab(opname,                             T , NumberArray<size MacroComma T2>, \
                  typecomparison##Type<NumberArray<size MacroComma T2> MacroComma T MacroComma true>) \
NumberArray_op_ab(opname, NumberArray<size MacroComma T>,                             T2 , \
                  typecomparison##Type<NumberArray<size MacroComma T> MacroComma T2>)

NumberArray_op(+,Plus)
NumberArray_op(-,Minus)
NumberArray_op(*,Multiplies)
NumberArray_op(/,Divides)

namespace std {

#define NumberArray_std_unary(funcname) \
template <std::size_t size, typename T> \
inline \
NumberArray<size, T> \
funcname (const NumberArray<size, T>& a) \
{ \
  NumberArray<size, T> returnval; \
 \
  for (unsigned int i=0; i != size; ++i) \
    returnval[i] = std::funcname(a[i]); \
 \
  return returnval; \
}


#define NumberArray_std_binary_abab(funcname, atype, btype, abtypes, aarg, barg) \
template <std::size_t size, typename T, typename T2> \
inline \
typename CompareTypes<abtypes>::supertype \
funcname (const atype& a, const btype& b) \
{ \
  typedef typename CompareTypes<abtypes>::supertype TS; \
  TS returnval; \
 \
  for (unsigned int i=0; i != size; ++i) \
    returnval[i] = std::funcname(aarg, barg); \
 \
  return returnval; \
}

#define NumberArray_std_binary(funcname) \
NumberArray_std_binary_abab(funcname, NumberArray<size MacroComma T>, NumberArray<size MacroComma T2>, \
                            NumberArray<size MacroComma T> MacroComma NumberArray<size MacroComma T2>, a[i], b[i]) \
NumberArray_std_binary_abab(funcname,                             T , NumberArray<size MacroComma T2>, \
                            NumberArray<size MacroComma T2> MacroComma T,                              a,    b[i]) \
NumberArray_std_binary_abab(funcname, NumberArray<size MacroComma T>,                             T2 , \
                            NumberArray<size MacroComma T> MacroComma T2,                              a[i],    b)

NumberArray_std_binary(pow)
NumberArray_std_unary(exp)
NumberArray_std_unary(log)
NumberArray_std_unary(log10)
NumberArray_std_unary(sin)
NumberArray_std_unary(cos)
NumberArray_std_unary(tan)
NumberArray_std_unary(asin)
NumberArray_std_unary(acos)
NumberArray_std_unary(atan)
NumberArray_std_binary(atan2)
NumberArray_std_unary(sinh)
NumberArray_std_unary(cosh)
NumberArray_std_unary(tanh)
NumberArray_std_unary(sqrt)
NumberArray_std_unary(abs)
NumberArray_std_binary(max)
NumberArray_std_binary(min)
NumberArray_std_unary(ceil)
NumberArray_std_unary(floor)
NumberArray_std_binary(fmod)


template <std::size_t size, typename T>
class numeric_limits<NumberArray<size, T> > : 
  public raw_numeric_limits<NumberArray<size, T>, T> {};

} // namespace std

#define NumberArray_operator_binary_abab(opname, atype, btype, aarg, barg) \
template <std::size_t size, typename T, typename T2> \
inline \
NumberArray<size, bool> \
operator opname (const atype& a, const btype& b) \
{ \
  NumberArray<size, bool> returnval; \
 \
  for (unsigned int i=0; i != size; ++i) \
    returnval[i] = (aarg opname barg); \
 \
  return returnval; \
}

#define NumberArray_operator_binary(opname) \
NumberArray_operator_binary_abab(opname, NumberArray<size MacroComma T>, NumberArray<size MacroComma T2>, a[i], b[i]) \
NumberArray_operator_binary_abab(opname,                             T , NumberArray<size MacroComma T2>, a,    b[i]) \
NumberArray_operator_binary_abab(opname, NumberArray<size MacroComma T>,                             T2 , a[i], b)

NumberArray_operator_binary(<)
NumberArray_operator_binary(<=)
NumberArray_operator_binary(>)
NumberArray_operator_binary(>=)
NumberArray_operator_binary(==)
NumberArray_operator_binary(!=)

template <std::size_t size, typename T>
inline
std::ostream&      
operator<< (std::ostream& output, const NumberArray<size,T>& a)
{
  output << '{';
  if (size)
    output << a[0];
  for (unsigned int i=1; i<size; ++i)
    output << ',' << a[i];
  output << '}';
  return output;
}


// CompareTypes, RawType specializations

#define NumberArray_comparisons(templatename) \
template<std::size_t size, typename T, bool reverseorder> \
struct templatename<NumberArray<size,T>, NumberArray<size,T>, reverseorder> { \
  typedef NumberArray<size, T> supertype; \
}; \
 \
template<std::size_t size, typename T, typename T2, bool reverseorder> \
struct templatename<NumberArray<size,T>, NumberArray<size,T2>, reverseorder> { \
  typedef NumberArray<size, typename Symmetric##templatename<T, T2, reverseorder>::supertype> supertype; \
}; \
 \
template<std::size_t size, typename T, typename T2, bool reverseorder> \
struct templatename<NumberArray<size, T>, T2, reverseorder, \
                    typename boostcopy::enable_if<BuiltinTraits<T2> >::type> { \
  typedef NumberArray<size, typename Symmetric##templatename<T, T2, reverseorder>::supertype> supertype; \
}

NumberArray_comparisons(CompareTypes);
NumberArray_comparisons(PlusType);
NumberArray_comparisons(MinusType);
NumberArray_comparisons(MultipliesType);
NumberArray_comparisons(DividesType);

template <std::size_t size, typename T>
struct RawType<NumberArray<size, T> >
{
  typedef NumberArray<size, typename RawType<T>::value_type> value_type;

  static value_type value(const NumberArray<size, T>& a)
    {
      value_type returnval;
      for (unsigned int i=0; i != size; ++i)
        returnval[i] = RawType<T>::value(a[i]);
      return returnval;
    }
};

#endif // __numberarray_h__
