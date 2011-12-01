
#ifndef __shadownumber_h__
#define __shadownumber_h__

#include <algorithm>
#include <iostream>

#include "compare_types.h"
#include "raw_type.h"

template <typename T, typename S>
class ShadowNumber
{
public:
  typedef T value_type;

  typedef S shadow_type;

  ShadowNumber() {}

  template <typename T2, typename S2>
  ShadowNumber(const T2& v, const S2& s) : _val(v), _shadow(s) {}

  template <typename T2>
  ShadowNumber(const T2& val) : _val(val), _shadow(val) {}

  template <typename T2, typename S2>
  ShadowNumber(ShadowNumber<T2, S2>& other) : _val(other._val), _shadow(other._shadow) {}

  T& value() { return _val; }

  const T& value() const { return _val; }

  S& shadow() { return _shadow; }

  const S& shadow() const { return _shadow; }

  ShadowNumber<T,S> operator- () { return ShadowNumber<T,S> (-_val, -_shadow); }

  template <typename T2, typename S2>
  ShadowNumber<T,S>& operator+= (const ShadowNumber<T2,S2>& a)
    { _val += a.value(); _shadow += a.shadow(); return *this; }

  template <typename T2>
  ShadowNumber<T,S>& operator+= (const T2& a)
    { _val += a; _shadow += a; return *this; }

  template <typename T2, typename S2>
  ShadowNumber<T,S>& operator-= (const ShadowNumber<T2,S2>& a)
    { _val -= a.value(); _shadow -= a.shadow(); return *this; }

  template <typename T2>
  ShadowNumber<T,S>& operator-= (const T2& a)
    { _val -= a; _shadow -= a; return *this; }

  template <typename T2, typename S2>
  ShadowNumber<T,S>& operator*= (const ShadowNumber<T2,S2>& a)
    { _val *= a.value(); _shadow *= a.shadow(); return *this; }

  template <typename T2>
  ShadowNumber<T,S>& operator*= (const T2& a)
    { _val *= a; _shadow *= a; return *this; }

  template <typename T2, typename S2>
  ShadowNumber<T,S>& operator/= (const ShadowNumber<T2,S2>& a)
    { _val /= a.value(); _shadow /= a.shadow(); return *this; }

  template <typename T2>
  ShadowNumber<T,S>& operator/= (const T2& a)
    { _val /= a; _shadow /= a; return *this; }

private:
  T _val;
  S _shadow;
};

//
// Non-member functions
//

#define ShadowNumber_op(opname) \
template <typename T, typename S, typename T2, typename S2> \
inline \
ShadowNumber<typename CompareTypes<T,T2>::supertype, \
             typename CompareTypes<S,S2>::supertype> \
operator opname (const ShadowNumber<T,S>& a, const ShadowNumber<T2,S2>& b) \
{ \
  typedef typename CompareTypes<T,T2>::supertype TS; \
  typedef typename CompareTypes<S,S2>::supertype SS; \
  ShadowNumber<TS, SS> returnval(a); \
  returnval opname##= b; \
  return returnval; \
} \
 \
template <typename T, typename S, typename T2> \
inline \
typename boostcopy::enable_if_c<ScalarTraits<T2>::value, \
ShadowNumber<typename CompareTypes<T,T2>::supertype, S> \
>::type \
operator opname (const ShadowNumber<T,S>& a, const T2& b) \
{ \
  typedef typename CompareTypes<T,T2>::supertype TS; \
  ShadowNumber<TS, S> returnval(a); \
  returnval opname##= b; \
  return returnval; \
 \
} \
template <typename T, typename T2, typename S> \
inline \
typename boostcopy::enable_if_c<ScalarTraits<T>::value, \
ShadowNumber<typename CompareTypes<T,T2>::supertype, S> \
>::type \
operator opname (const T& a, const ShadowNumber<T2,S>& b) \
{ \
  typedef typename CompareTypes<T,T2>::supertype TS; \
  ShadowNumber<TS, S> returnval(a); \
  returnval opname##= b; \
  return returnval; \
}

ShadowNumber_op(+)
ShadowNumber_op(-)
ShadowNumber_op(*)
ShadowNumber_op(/)

namespace std {

// Some forward declarations necessary for recursive DualNumbers

template <typename T, typename S>
inline ShadowNumber<T,S> cos   (const ShadowNumber<T,S>& a);

template <typename T, typename S>
inline ShadowNumber<T,S> cosh  (const ShadowNumber<T,S>& a);

// Now just combined declaration/definitions

#define ShadowNumber_std_unary(funcname) \
template <typename T, typename S> \
inline \
ShadowNumber<T, S> \
funcname (const ShadowNumber<T, S>& a) \
{ \
  return ShadowNumber<T,S>(std::funcname(a.value()), std::funcname(a.shadow())); \
}


#define ShadowNumber_std_binary(funcname) \
template <typename T, typename S, typename T2, typename S2> \
inline \
ShadowNumber<typename CompareTypes<T,T2>::supertype, \
             typename CompareTypes<S,S2>::supertype> \
funcname (const ShadowNumber<T,S>& a, const ShadowNumber<T2,S2>& b) \
{ \
  typedef typename CompareTypes<T,T2>::supertype TS; \
  typedef typename CompareTypes<S,S2>::supertype SS; \
  return ShadowNumber<TS, SS> (std::funcname(a.value(), b.value()), \
                               std::funcname(a.shadow(), b.shadow())); \
} \
 \
template <typename T, typename S, typename T2> \
inline \
ShadowNumber<typename CompareTypes<T,T2>::supertype, S> \
funcname (const ShadowNumber<T,S>& a, const T2& b) \
{ \
  typedef typename CompareTypes<T,T2>::supertype TS; \
  return ShadowNumber<TS, S> (std::funcname(a.value(), b), \
                              std::funcname(a.shadow(), b)); \
} \
 \
template <typename T, typename T2, typename S> \
inline \
ShadowNumber<typename CompareTypes<T,T2>::supertype, S> \
funcname (const T& a, const ShadowNumber<T2,S>& b) \
{ \
  typedef typename CompareTypes<T,T2>::supertype TS; \
  return ShadowNumber<TS, S> (std::funcname(a, b.value()), \
                              std::funcname(a, b.shadow())); \
}


ShadowNumber_std_binary(pow)
ShadowNumber_std_unary(exp)
ShadowNumber_std_unary(log)
ShadowNumber_std_unary(log10)
ShadowNumber_std_unary(sin)
ShadowNumber_std_unary(cos)
ShadowNumber_std_unary(tan)
ShadowNumber_std_unary(asin)
ShadowNumber_std_unary(acos)
ShadowNumber_std_unary(atan)
ShadowNumber_std_binary(atan2)
ShadowNumber_std_unary(sinh)
ShadowNumber_std_unary(cosh)
ShadowNumber_std_unary(tanh)
ShadowNumber_std_unary(sqrt)
ShadowNumber_std_unary(abs)
ShadowNumber_std_binary(max)
ShadowNumber_std_binary(min)
ShadowNumber_std_unary(ceil)
ShadowNumber_std_unary(floor)
ShadowNumber_std_binary(fmod)

template <typename T, typename S>
class numeric_limits<ShadowNumber<T, S> > :
  public raw_numeric_limits<ShadowNumber<T, S>, T> {};

} // namespace std

#define ShadowNumber_operator_binary(opname) \
template <typename T, typename S, typename T2, typename S2> \
inline \
ShadowNumber<bool, bool> \
operator opname (const ShadowNumber<T,S>& a, const ShadowNumber<T2,S2>& b) \
{ \
  return ShadowNumber<bool, bool> (a.value() opname b.value(), a.shadow() opname b.shadow()); \
} \
 \
template <typename T, typename S, typename T2> \
inline \
ShadowNumber<bool, bool> \
operator opname (const ShadowNumber<T,S>& a, const T2& b) \
{ \
  return ShadowNumber<bool, bool> (a.value() opname b, a.shadow() opname b); \
} \
 \
template <typename T, typename T2, typename S> \
inline \
ShadowNumber<bool, bool> \
operator opname (const T& a, const ShadowNumber<T2,S>& b) \
{ \
  return ShadowNumber<bool, bool> (a opname b.value(), a opname b.shadow()); \
}


ShadowNumber_operator_binary(<)
ShadowNumber_operator_binary(<=)
ShadowNumber_operator_binary(>)
ShadowNumber_operator_binary(>=)
ShadowNumber_operator_binary(==)
ShadowNumber_operator_binary(!=)

template <typename T, typename S>
inline
std::ostream&      
operator<< (std::ostream& output, const ShadowNumber<T,S>& a)
{
  return output << '(' << a.value() << ',' << a.shadow() << ')';
}


// ScalarTraits, RawType, CompareTypes specializations

template <typename T, typename S>
struct ScalarTraits<ShadowNumber<T, S> >
{
  static const bool value = ScalarTraits<T>::value;
};

template<typename T, typename S>
struct CompareTypes<ShadowNumber<T,S>, ShadowNumber<T,S> > {
  typedef ShadowNumber<T, S> supertype;
};

template<typename T, typename S, typename T2, typename S2>
struct CompareTypes<ShadowNumber<T,S>, ShadowNumber<T2,S2> > {
  typedef ShadowNumber<typename CompareTypes<T, T2>::supertype,
                       typename CompareTypes<S, S2>::supertype> supertype;
};

template<typename T, typename S, typename T2>
struct CompareTypes<ShadowNumber<T, S>, T2> {
  typedef ShadowNumber<typename CompareTypes<T, T2>::supertype, S> supertype;
};

template<typename T, typename T2, typename S>
struct CompareTypes<T, ShadowNumber<T2, S> > {
  typedef ShadowNumber<typename CompareTypes<T, T2>::supertype, S> supertype;
};



template <typename T, typename S>
struct RawType<ShadowNumber<T, S> >
{
  typedef typename RawType<T>::value_type value_type;

  static value_type value(const ShadowNumber<T, S>& a) {
    const S max_value = std::max(S(a.value()), a.shadow());
    if (max_value) {
      const S relative_error = (a.value() - a.shadow()) / max_value;
      if (relative_error > 10*std::numeric_limits<T>::epsilon())
        std::cerr << "Shadow relative error = " << relative_error << std::endl;
    }
    return a.value();
  }
};

#endif // __shadownumber_h__
