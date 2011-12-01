
#ifndef __dualnumber_h__
#define __dualnumber_h__

#include <ostream>

#include "compare_types.h"
#include "raw_type.h"

template <typename T, typename D=T>
class DualNumber
{
public:
  typedef T value_type;

  typedef D derivatives_type;

  DualNumber();

  template <typename T2>
  DualNumber(const T2& val);

  template <typename T2, typename D2>
  DualNumber(const T2& val, const D2& deriv);

  T& value() { return _val; }

  const T& value() const { return _val; }

  D& derivatives() { return _deriv; }

  const D& derivatives() const { return _deriv; }

  DualNumber<T,D> operator- () { return DualNumber<T,D>(-_val, -_deriv); }

  template <typename T2, typename D2>
  DualNumber<T,D>& operator+= (const DualNumber<T2,D2>& a);

  template <typename T2>
  DualNumber<T,D>& operator+= (const T2& a);

  template <typename T2, typename D2>
  DualNumber<T,D>& operator-= (const DualNumber<T2,D2>& a);

  template <typename T2>
  DualNumber<T,D>& operator-= (const T2& a);

  template <typename T2, typename D2>
  DualNumber<T,D>& operator*= (const DualNumber<T2,D2>& a);

  template <typename T2>
  DualNumber<T,D>& operator*= (const T2& a);

  template <typename T2, typename D2>
  DualNumber<T,D>& operator/= (const DualNumber<T2,D2>& a);

  template <typename T2>
  DualNumber<T,D>& operator/= (const T2& a);


private:
  T _val;
  D _deriv;
};



// Helper class to handle partial specialization for DualNumber
// constructors

template <typename T, typename D>
struct DualNumberConstructor
{
  static T value(const DualNumber<T,D>& v) { return v.value(); }

  template <typename T2>
  static T value(const T2& v) { return v; }

  template <typename T2, typename D2>
  static T value(const T2& v, const D2&) { return v; }

  template <typename T2, typename D2>
  static T value(const DualNumber<T2,D2>& v) { return v.value(); }

  template <typename T2>
  static D deriv(const T2&) { return 0.; }

  template <typename T2, typename D2>
  static D deriv(const DualNumber<T2,D2>& v) { return v.derivatives(); }

  template <typename T2, typename D2>
  static D deriv(const T2&, const D2& d) { return d; }
};

template <typename T, typename D, typename DD>
struct DualNumberConstructor<DualNumber<T,D>, DD>
{
  static T value(const DualNumber<DualNumber<T,D>, DD>& v) { return v.value(); }

  template <typename T2>
  static DualNumber<T,D> value(const T2& v) { return v; }

  template <typename T2, typename D2>
  static DualNumber<T,D> value(const T2& v, const D2& d) { return DualNumber<T,D>(v,d); }

  template <typename D2>
  static DualNumber<T,D> value(const DualNumber<T,D>& v, const D2&) { return v; }

  template <typename T2>
  static DD deriv(const T2&) { return 0.; }

  template <typename T2, typename D2>
  static DD deriv(const DualNumber<T2,D2>& v) { return v.derivatives(); }

  template <typename T2, typename D2>
  static DD deriv(const T2&, const D2& d) { return d; }
};


//
// Member function definitions
//

template <typename T, typename D>
inline
DualNumber<T,D>::DualNumber() :
  _val(0.), _deriv(0.) {}

template <typename T, typename D>
template <typename T2>
inline
DualNumber<T,D>::DualNumber(const T2& val) :
  _val  (DualNumberConstructor<T,D>::value(val)),
  _deriv(DualNumberConstructor<T,D>::deriv(val)) {}

template <typename T, typename D>
template <typename T2, typename D2>
inline
DualNumber<T,D>::DualNumber(const T2& val,
                            const D2& deriv) :
  _val  (DualNumberConstructor<T,D>::value(val,deriv)),
  _deriv(DualNumberConstructor<T,D>::deriv(val,deriv)) {}

#define DualNumber_op(opname, simplecalc, dualcalc) \
template <typename T, typename D> \
template <typename T2> \
inline \
DualNumber<T,D>& \
DualNumber<T,D>::operator opname##= (const T2& in) \
{ \
  simplecalc; \
  this->value() opname##= in; \
  return *this; \
} \
 \
template <typename T, typename D> \
template <typename T2, typename D2> \
inline \
DualNumber<T,D>& \
DualNumber<T,D>::operator opname##= (const DualNumber<T2,D2>& in) \
{ \
  dualcalc; \
  this->value() opname##= in.value(); \
  return *this; \
} \
 \
template <typename T, typename D, typename T2, typename D2> \
inline \
DualNumber<typename CompareTypes<T,T2>::supertype,  \
           typename CompareTypes<D,D2>::supertype>  \
operator opname (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b) \
{ \
  typedef typename CompareTypes<T,T2>::supertype TS; \
  typedef typename CompareTypes<D,D2>::supertype DS; \
  DualNumber<TS,DS> returnval = a; \
  returnval opname##= b; \
  return returnval; \
} \
 \
template <typename T, typename T2, typename D> \
inline \
typename boostcopy::enable_if_c<ScalarTraits<T>::value, \
DualNumber<typename CompareTypes<T,T2>::supertype, D> \
>::type \
operator opname (const T& a, const DualNumber<T2,D>& b) \
{ \
  typedef typename CompareTypes<T,T2>::supertype TS; \
  DualNumber<TS,D> returnval = a; \
  returnval opname##= b; \
  return returnval; \
} \
 \
template <typename T, typename D, typename T2> \
inline \
typename boostcopy::enable_if_c<ScalarTraits<T2>::value, \
DualNumber<typename CompareTypes<T,T2>::supertype, D> \
>::type \
operator opname (const DualNumber<T,D>& a, const T2& b) \
{ \
  typedef typename CompareTypes<T,T2>::supertype TS; \
  DualNumber<TS,D> returnval = a; \
  returnval opname##= b; \
  return returnval; \
}



DualNumber_op(+, , this->derivatives() += in.derivatives())

DualNumber_op(-, , this->derivatives() -= in.derivatives())

DualNumber_op(*, this->derivatives() *= in,
  this->derivatives() *= in.value();
  this->derivatives() += this->value() * in.derivatives();)

DualNumber_op(/, this->derivatives() /= in,
  this->derivatives() /= in.value();
  this->derivatives() -= this->value()/(in.value()*in.value()) * in.derivatives();
)



namespace std {

// Some forward declarations necessary for recursive DualNumbers

template <typename T, typename D>
inline DualNumber<T,D> cos   (const DualNumber<T,D>& a);

template <typename T, typename D>
inline DualNumber<T,D> cosh  (const DualNumber<T,D>& a);

// Now just combined declaration/definitions

#define DualNumber_std_unary(funcname, derivative, precalc) \
template <typename T, typename D> \
inline \
DualNumber<T,D> funcname   (const DualNumber<T,D>& in) \
{ \
  T funcval = std::funcname(in.value()); \
  precalc; \
  return DualNumber<T,D>(funcval, (derivative) * in.derivatives()); \
}

DualNumber_std_unary(sqrt, 1 / (2 * funcval),)
DualNumber_std_unary(exp, funcval,)
DualNumber_std_unary(log, 1 / in.value(),)
DualNumber_std_unary(log10, 1 / in.value() * (1/std::log(T(10.))),)
DualNumber_std_unary(sin, std::cos(in.value()),)
DualNumber_std_unary(cos, -std::sin(in.value()),)
DualNumber_std_unary(tan, sec_in * sec_in, T sec_in = 1 / std::cos(in.value()))
DualNumber_std_unary(asin, 1 / std::sqrt(1 - in.value()*in.value()),)
DualNumber_std_unary(acos, -1 / std::sqrt(1 - in.value()*in.value()),)
DualNumber_std_unary(atan, 1 / (1 + in.value()*in.value()),)
DualNumber_std_unary(sinh, std::cosh(in.value()),)
DualNumber_std_unary(cosh, std::sinh(in.value()),)
DualNumber_std_unary(tanh, sech_in * sech_in, T sech_in = 1 / std::cos(in.value()))
DualNumber_std_unary(abs, (in.value() > 0) - (in.value() < 0),) // std < and > return 0 or 1
DualNumber_std_unary(ceil, 0,)
DualNumber_std_unary(floor, 0,)

#define DualNumber_std_binary(funcname, derivative) \
template <typename T, typename D, typename T2, typename D2> \
inline \
DualNumber<typename CompareTypes<T,T2>::supertype, \
           typename CompareTypes<D,D2>::supertype> \
funcname (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b) \
{ \
  typedef typename CompareTypes<T,T2>::supertype TS; \
  typedef typename CompareTypes<D,D2>::supertype DS; \
 \
  TS funcval = std::funcname(a.value(), b.value()); \
  return DualNumber<TS,DS>(funcval, derivative); \
} \
 \
template <typename T, typename T2, typename D> \
inline \
typename boostcopy::enable_if_c<ScalarTraits<T>::value, \
DualNumber<typename CompareTypes<T,T2>::supertype, D> \
>::type \
funcname (const T& a, const DualNumber<T2,D>& b) \
{ \
  typedef typename CompareTypes<T,T2>::supertype TS; \
  DualNumber<TS, D> newa(a); \
  return std::funcname(newa, b); \
} \
 \
template <typename T, typename T2, typename D> \
inline \
typename boostcopy::enable_if_c<ScalarTraits<T2>::value, \
DualNumber<typename CompareTypes<T,T2>::supertype, D> \
>::type \
funcname (const DualNumber<T,D>& a, const T2& b) \
{ \
  typedef typename CompareTypes<T,T2>::supertype TS; \
  DualNumber<TS, D> newb(b); \
  return std::funcname(a, newb); \
}

DualNumber_std_binary(pow, 
  funcval * (b.value() * a.derivatives() / a.value() + b.derivatives() * std::log(a.value())))
DualNumber_std_binary(atan2,
  (b.value() * a.derivatives() - a.value() * b.derivatives()) /
  (b.value() * b.value() + a.value() * a.value()))
DualNumber_std_binary(max,
  (a.value() > b.value()) ?  a : b)
DualNumber_std_binary(min,
  (a.value() > b.value()) ?  b : a)
DualNumber_std_binary(fmod, a.derivatives())

template <typename T, typename D>
class numeric_limits<DualNumber<T, D> > : 
  public raw_numeric_limits<DualNumber<T, D>, T> {};

} // namespace std

#define DualNumber_compare(opname) \
template <typename T, typename D, typename T2, typename D2> \
inline \
bool \
operator opname  (const DualNumber<T,D>& a, const DualNumber<T2,D2>& b) \
{ \
  return (a.value() opname b.value()); \
} \
 \
template <typename T, typename T2, typename D2> \
inline \
typename boostcopy::enable_if_c<ScalarTraits<T>::value, \
bool \
>::type \
operator opname  (const T& a, const DualNumber<T2,D2>& b) \
{ \
  return (a opname b.value()); \
} \
 \
template <typename T, typename T2, typename D> \
inline \
typename boostcopy::enable_if_c<ScalarTraits<T2>::value, \
bool \
>::type \
operator opname  (const DualNumber<T,D>& a, const T2& b) \
{ \
  return (a.value() opname b.value()); \
}

DualNumber_compare(>)
DualNumber_compare(>=)
DualNumber_compare(<)
DualNumber_compare(<=)
DualNumber_compare(==)
DualNumber_compare(!=)

template <typename T, typename D>
inline
std::ostream&      
operator<< (std::ostream& output, const DualNumber<T,D>& a)
{
  return output << '(' << a.value() << ',' << a.derivatives() << ')';
}


// ScalarTraits, RawType, CompareTypes specializations

template <typename T, typename D>
struct ScalarTraits<DualNumber<T, D> >
{
  static const bool value = ScalarTraits<T>::value;
};

template <typename T, typename D>
struct RawType<DualNumber<T, D> >
{
  typedef typename RawType<T>::value_type value_type;

  static value_type value(const DualNumber<T, D>& a) { return raw_value(a.value()); }
};

template<typename T, typename T2, typename D>
struct CompareTypes<T, DualNumber<T2, D> > {
  typedef DualNumber<typename CompareTypes<T, T2>::supertype,
                     typename CompareTypes<D, T2>::supertype> supertype;
};

template<typename T, typename T2, typename D>
struct CompareTypes<DualNumber<T, D>, T2> {
  typedef DualNumber<typename CompareTypes<T, T2>::supertype,
                     typename CompareTypes<D, T2>::supertype> supertype;
};

template<typename T, typename D, typename T2, typename D2>
struct CompareTypes<DualNumber<T, D>, DualNumber<T2, D2> > {
  typedef DualNumber<typename CompareTypes<T, T2>::supertype,
                     typename CompareTypes<D, D2>::supertype> supertype;
};

template<typename T, typename D>
struct CompareTypes<DualNumber<T, D>, DualNumber<T, D> > {
  typedef DualNumber<T, D> supertype;
};

#endif // __dualnumber_h__
