#ifndef __raw_type_h__
#define __raw_type_h__

#include <limits>


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


#endif // __raw_type_h__
