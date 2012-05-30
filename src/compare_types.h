
// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#ifndef __compare_types_h__
#define __compare_types_h__

// System includes
#include <complex>

// Compile-time assertions:
// Use of ctassert<E>, where E is a constant expression,
// will cause a compile-time error unless E evaulates to
// a nonzero integral value.

template <bool t>
struct ctassert {
  enum { N = 1 - 2 * int(!t) };
    // 1 if t is true, -1 if t is false.
  static char A[N];

  static void apply() {}
};

template <typename T1=void,
          typename T2=void,
          typename T3=void,
          typename T4=void,
          typename T5=void,
          typename T6=void,
          typename T7=void,
          typename T8=void,
          typename T9=void>
struct ctprint
{
  typedef typename T9::compiler_should_print_types force_an_error;
  force_an_error f;
  static void apply() {}
};

template <>
struct ctprint<>
{
  static void apply() {}
};

template <bool t>
char ctassert<t>::A[N];

// Copy of boost's enable_if_c

namespace boostcopy {
  template <bool B, class T = void>
    struct enable_if_c {
      typedef T type;
    };

  template <class T>
    struct enable_if_c<false, T> {};

  template <class Cond, class T = void>
    struct enable_if : public enable_if_c<Cond::value, T> {};

  template <class B, class T>
    struct enable_if_class {
      typedef T type;
    };

}

// We can pass commas for template arguments through one level of
// macro expansion by "escaping" them this way:
#define MacroComma ,


// List of scalar and builtin classes, useful for disambiguation
template <typename T>
struct ScalarTraits {
      static const bool value = false;
};

template <typename T>
struct BuiltinTraits {
      static const bool value = false;
};

#define ScalarBuiltin_true(type) \
template<> \
struct ScalarTraits<type> { static const bool value = true; }; \
template<> \
struct BuiltinTraits<type> { static const bool value = true; }

ScalarBuiltin_true(char);
ScalarBuiltin_true(short);
ScalarBuiltin_true(int);
ScalarBuiltin_true(long);
ScalarBuiltin_true(unsigned char);
ScalarBuiltin_true(unsigned short);
ScalarBuiltin_true(unsigned int);
ScalarBuiltin_true(unsigned long);
ScalarBuiltin_true(float);
ScalarBuiltin_true(double);
ScalarBuiltin_true(long double);

template<typename T>
struct ScalarTraits<std::complex<T> > { static const bool value = ScalarTraits<T>::value; };

template<typename T>
struct BuiltinTraits<std::complex<T> > { static const bool value = BuiltinTraits<T>::value; };



// Operators using different but compatible types need a return value
// based on whichever type the other can be upconverted into.  For
// instance, the proper return type for
// TypeVector<float>::operator*(double) is TypeVector<double>.  In
// general, an operation using types S and T should return a value
// based on CompareTypes<S,T>::supertype

template<typename S, typename T, bool reverseorder=false, typename Enable=void>
struct CompareTypes {

// All specializations need to define supertype.  But compilers give
// better error messages for forgot-to-specialize-CompareTypes bugs if
// we leave supertype undefined in the unspecialized template
// definition.

// typedef void supertype;

//   typedef S nosupertype;
};


// A tricky SFINAE idiom for testing whether a particular CompareTypes
// combination has been properly specialized:

template<typename T>
struct DefinesSupertype
{
private:
    typedef char                      yes;
    typedef struct { char array[2]; } no;

    template<typename C> static yes test(typename C::supertype*);
    template<typename C> static no  test(...);
public:
    static const bool value = (sizeof(test<T>(0)) == sizeof(yes));
};


// PlusType, MinusType, MultipliesType, and DividesType are usually just the
// same as CompareTypes, but user types may want to specialize further for
// efficiency.
template<typename S, typename T, bool reverseorder=false, typename Enable=void>
struct PlusType {
};

template<typename S, typename T, bool reverseorder=false, typename Enable=void>
struct MinusType {
};

template<typename S, typename T, bool reverseorder=false, typename Enable=void>
struct MultipliesType {
};

template<typename S, typename T, bool reverseorder=false, typename Enable=void>
struct DividesType {
};


// typenames may need a MacroComma...
#define CompareTypes_default_Type(functor, typenames, typename1, typename2, enabletype) \
template<typenames bool reverseorder> \
struct functor##Type<typename1, typename2, reverseorder, enabletype> \
{ \
   typedef typename CompareTypes<typename1,typename2,reverseorder>::supertype supertype; \
}

#define CompareTypes_default_Types(typenames,typename1,typename2, enabletype) \
CompareTypes_default_Type(Plus,typenames,typename1,typename2, enabletype); \
CompareTypes_default_Type(Minus,typenames,typename1,typename2, enabletype); \
CompareTypes_default_Type(Multiplies,typenames,typename1,typename2, enabletype); \
CompareTypes_default_Type(Divides,typenames,typename1,typename2, enabletype) \


template<bool reverseorder>
struct CompareTypes<void, void, reverseorder> {
  typedef void supertype;
};

template<typename T, bool reverseorder>
struct CompareTypes<T, void, reverseorder> {
  typedef T supertype;
};

template<typename T, bool reverseorder>
struct CompareTypes<T, T, reverseorder> {
  typedef T supertype;
};

template<typename T, bool reverseorder>
struct CompareTypes<T, std::complex<T>, reverseorder> {
  typedef std::complex<T> supertype;
};

template<typename T, bool reverseorder>
struct CompareTypes<std::complex<T>, T, reverseorder> {
  typedef std::complex<T> supertype;
};

// There's got to be some magic template way to do these better - but the best
// thing on the net requires a bunch of Alexandrescu's code and doesn't work
// with older compilers

#define CompareTypes_super(a,b,super) \
	template<bool reverseorder> \
	struct CompareTypes<a, b, reverseorder> { \
	  typedef super supertype; \
	}; \
        CompareTypes_default_Types(,a,b,void)

#define CompareTypes_all(mysub,mysuper) \
        CompareTypes_super(mysub, mysuper, mysuper); \
        CompareTypes_super(mysuper, mysub, mysuper); \
        CompareTypes_super(std::complex<mysub>, mysuper, std::complex<mysuper>); \
        CompareTypes_super(mysuper, std::complex<mysub>, std::complex<mysuper>); \
        CompareTypes_super(mysub, std::complex<mysuper>, std::complex<mysuper>); \
        CompareTypes_super(std::complex<mysuper>, mysub, std::complex<mysuper>); \
        CompareTypes_super(std::complex<mysub>, std::complex<mysuper>, std::complex<mysuper>); \
        CompareTypes_super(std::complex<mysuper>, std::complex<mysub>, std::complex<mysuper>)

#define CompareTypes_single(mytype) \
        CompareTypes_super(mytype, mytype, mytype)

CompareTypes_single(unsigned char);
CompareTypes_single(unsigned short);
CompareTypes_single(unsigned int);
CompareTypes_single(char);
CompareTypes_single(short);
CompareTypes_single(int);
CompareTypes_single(float);
CompareTypes_single(double);
CompareTypes_single(long double);

CompareTypes_all(unsigned char, short);
CompareTypes_all(unsigned char, int);
CompareTypes_all(unsigned char, float);
CompareTypes_all(unsigned char, double);
CompareTypes_all(unsigned char, long double);
CompareTypes_all(unsigned short, int);
CompareTypes_all(unsigned short, float);
CompareTypes_all(unsigned short, double);
CompareTypes_all(unsigned short, long double);
CompareTypes_all(unsigned int, float);
CompareTypes_all(unsigned int, double);
CompareTypes_all(unsigned int, long double);
CompareTypes_all(char, short);
CompareTypes_all(char, int);
CompareTypes_all(char, float);
CompareTypes_all(char, double);
CompareTypes_all(char, long double);
CompareTypes_all(short, int);
CompareTypes_all(short, float);
CompareTypes_all(short, double);
CompareTypes_all(short, long double);
CompareTypes_all(int, float);
CompareTypes_all(int, double);
CompareTypes_all(int, long double);
CompareTypes_all(float, double);
CompareTypes_all(float, long double);
CompareTypes_all(double, long double);

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


// We can define CompareTypes template specializations with user types
// asymmetrically, to assist in disambiguation of templated functions
// or classes.  But sometimes we just need the supertype:

// FIXME: this won't work yet for cases where CompareTypes depends on
// reverseorder

#define Symmetric_definition(templatename) \
template<typename S, typename T, bool reverseorder=false, \
	 bool Simple=DefinesSupertype<templatename<S,T> >::value> \
struct Symmetric##templatename { \
  typedef typename templatename<S,T,reverseorder>::supertype supertype; \
}; \
 \
template<typename S, typename T, bool reverseorder> \
struct Symmetric##templatename<S, T, reverseorder, false> \
{ \
  typedef typename templatename<T,S,!reverseorder>::supertype supertype; \
}



Symmetric_definition(CompareTypes);
Symmetric_definition(PlusType);
Symmetric_definition(MinusType);
Symmetric_definition(MultipliesType);
Symmetric_definition(DividesType);

#endif // __compare_types_h__
