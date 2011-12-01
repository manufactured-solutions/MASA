#ifndef __dualshadowarray_h__
#define __dualshadowarray_h__


#include "shadownumber.h"
#include "dualnumberarray.h"

// The compiler can't choose between NumberArray-vs-T and
// T-vs-ShadowNumber CompareTypes specializations...

template<std::size_t size, typename T, typename T2, typename S>
struct CompareTypes<NumberArray<size, T2>, ShadowNumber<T, S> > {
  typedef NumberArray<size, typename CompareTypes<T2, ShadowNumber<T, S> >::supertype> supertype;
};

template<std::size_t size, typename T, typename T2, typename S>
struct CompareTypes<ShadowNumber<T, S>, NumberArray<size, T2> > {
  typedef NumberArray<size, typename CompareTypes<T2, ShadowNumber<T, S> >::supertype> supertype;
};

template<typename T, typename D, typename T2, typename S>
struct CompareTypes<DualNumber<T, D>, ShadowNumber<T2, S> > {
  typedef DualNumber<typename CompareTypes<T, ShadowNumber<T2, S> >::supertype,
                     typename CompareTypes<D, ShadowNumber<T2, S> >::supertype> supertype;
};

template<typename T, typename D, typename T2, typename S>
struct CompareTypes<ShadowNumber<T2, S>, DualNumber<T, D> > {
  typedef DualNumber<typename CompareTypes<T, ShadowNumber<T2, S> >::supertype,
                     typename CompareTypes<D, ShadowNumber<T2, S> >::supertype> supertype;
};



#endif // __dualshadowarray_h__
