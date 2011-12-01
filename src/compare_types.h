
#ifndef __compare_types_h__
#define __compare_types_h__

// System includes
#include <complex>
#include <limits>

// Copy of boost's enable_if_c

namespace boostcopy {
  template <bool B, class T = void>
    struct enable_if_c {
      typedef T type;
    };

  template <class T>
    struct enable_if_c<false, T> {};
}



// Complete list of scalar classes, needed for disambiguation
template <typename T>
struct ScalarTraits {
      static const bool value = false;
};

#define ScalarTraits_true(type) \
template<> \
struct ScalarTraits<type> { static const bool value = true; }

ScalarTraits_true(char);
ScalarTraits_true(short);
ScalarTraits_true(int);
ScalarTraits_true(long);
ScalarTraits_true(unsigned char);
ScalarTraits_true(unsigned short);
ScalarTraits_true(unsigned int);
ScalarTraits_true(unsigned long);
ScalarTraits_true(float);
ScalarTraits_true(double);
ScalarTraits_true(long double);

template<typename T>
struct ScalarTraits<std::complex<T> > { static const bool value = ScalarTraits<T>::value; };

// Operators using different but compatible types need a return value
// based on whichever type the other can be upconverted into.  For
// instance, the proper return type for
// TypeVector<float>::operator*(double) is TypeVector<double>.  In
// general, an operation using types S and T should return a value
// based on CompareTypes<S,T>::supertype

template<typename S, typename T>
struct CompareTypes {
  typedef void supertype;
};

template<typename T>
struct CompareTypes<T, T> {
  typedef T supertype;
};

template<typename T>
struct CompareTypes<T, std::complex<T> > {
  typedef std::complex<T> supertype;
};

template<typename T>
struct CompareTypes<std::complex<T>, T> {
  typedef std::complex<T> supertype;
};

// There's got to be some magic template way to do these better - but the best
// thing on the net requires a bunch of Alexandrescu's code and doesn't work
// with older compilers

#define CompareTypes_super(a,b,super) \
	template<> \
	struct CompareTypes<a, b> { \
	  typedef super supertype; \
	}

#define SUPERTYPE(mysub,mysuper) \
        CompareTypes_super(mysub, mysuper, mysuper); \
        CompareTypes_super(mysuper, mysub, mysuper); \
        CompareTypes_super(std::complex<mysub>, mysuper, std::complex<mysuper>); \
        CompareTypes_super(mysuper, std::complex<mysub>, std::complex<mysuper>); \
        CompareTypes_super(mysub, std::complex<mysuper>, std::complex<mysuper>); \
        CompareTypes_super(std::complex<mysuper>, mysub, std::complex<mysuper>); \
        CompareTypes_super(std::complex<mysub>, std::complex<mysuper>, std::complex<mysuper>); \
        CompareTypes_super(std::complex<mysuper>, std::complex<mysub>, std::complex<mysuper>)

SUPERTYPE(unsigned char, short);
SUPERTYPE(unsigned char, int);
SUPERTYPE(unsigned char, float);
SUPERTYPE(unsigned char, double);
SUPERTYPE(unsigned char, long double);
SUPERTYPE(unsigned short, int);
SUPERTYPE(unsigned short, float);
SUPERTYPE(unsigned short, double);
SUPERTYPE(unsigned short, long double);
SUPERTYPE(unsigned int, float);
SUPERTYPE(unsigned int, double);
SUPERTYPE(unsigned int, long double);
SUPERTYPE(char, short);
SUPERTYPE(char, int);
SUPERTYPE(char, float);
SUPERTYPE(char, double);
SUPERTYPE(char, long double);
SUPERTYPE(short, int);
SUPERTYPE(short, float);
SUPERTYPE(short, double);
SUPERTYPE(short, long double);
SUPERTYPE(int, float);
SUPERTYPE(int, double);
SUPERTYPE(int, long double);
SUPERTYPE(float, double);
SUPERTYPE(float, long double);
SUPERTYPE(double, long double);

// gcc can't tell which of the following is the most specialized?  Weak.
/*
template<typename S, typename T>
struct CompareTypes<std::complex<S>, std::complex<T> > {
  typedef std::complex<typename CompareTypes<S, T>::supertype> supertype;
};

template<typename S, typename T>
struct CompareTypes<std::complex<S>, T> {
  typedef std::complex<typename CompareTypes<S, T>::supertype> supertype;
};

template<typename S, typename T>
struct CompareTypes<S, std::complex<T> > {
  typedef std::complex<typename CompareTypes<S, T>::supertype> supertype;
};
*/


/*  
 *
 *  Raw Types
 *
 *
 */

template <typename T>
struct RawType
{
  typedef T value_type;

  static value_type value(const T& a) { return a; }
};

// Make the user syntax slightly nicer
template <typename T>
inline
typename RawType<T>::value_type
raw_value(const T& a) { return RawType<T>::value(a); }


template <typename NewType, typename OldType>
class raw_numeric_limits
{
public:
  static const bool is_specialized = true;
  static OldType min() throw() { return NewType(std::numeric_limits<OldType>::min()); }
  static OldType max() throw() { return NewType(std::numeric_limits<OldType>::max()); }
  static const int  digits = std::numeric_limits<OldType>::digits;
  static const int  digits10 = std::numeric_limits<OldType>::digits10;
  static const bool is_signed = std::numeric_limits<OldType>::is_signed;
  static const bool is_integer = std::numeric_limits<OldType>::is_integer;
  static const bool is_exact = std::numeric_limits<OldType>::is_exact;
  static const int radix = std::numeric_limits<OldType>::radix;
  static OldType epsilon() throw() {return NewType(std::numeric_limits<OldType>::epsilon()); }
  static OldType round_error() throw() {return NewType(std::numeric_limits<OldType>::round_error()); }

  static const int  min_exponent = std::numeric_limits<OldType>::min_exponent;
  static const int  min_exponent10 = std::numeric_limits<OldType>::min_exponent10;
  static const int  max_exponent = std::numeric_limits<OldType>::max_exponent;
  static const int  max_exponent10 = std::numeric_limits<OldType>::max_exponent10;

  static const bool has_infinity = std::numeric_limits<OldType>::has_infinity;
  static const bool has_quiet_NaN = std::numeric_limits<OldType>::has_quiet_NaN;
  static const bool has_signaling_NaN = std::numeric_limits<OldType>::has_signaling_NaN;
  static const std::float_denorm_style has_denorm = std::numeric_limits<OldType>::has_denorm;
  static const bool has_denorm_loss = std::numeric_limits<OldType>::has_denorm_loss;
  static OldType infinity() throw() {return NewType(std::numeric_limits<OldType>::infinity()); }
  static OldType quiet_NaN() throw() {return NewType(std::numeric_limits<OldType>::quiet_NaN()); }
  static OldType signaling_NaN() throw() {return NewType(std::numeric_limits<OldType>::signaling_NaN()); }
  static OldType denorm_min() throw() {return NewType(std::numeric_limits<OldType>::denorm_min()); }

  static const bool is_iec559 = std::numeric_limits<OldType>::is_iec559;
  static const bool is_bounded = std::numeric_limits<OldType>::is_bounded;
  static const bool is_modulo = std::numeric_limits<OldType>::is_modulo;

  static const bool traps = std::numeric_limits<OldType>::traps;
  static const bool tinyness_before = std::numeric_limits<OldType>::tinyness_before;
  static const std::float_round_style round_style = std::numeric_limits<OldType>::round_style;
};


#endif // __compare_types_h__
