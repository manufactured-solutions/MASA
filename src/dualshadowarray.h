#ifndef __dualshadowarray_h__
#define __dualshadowarray_h__


// Order of declarations is important here?
#include "dualshadow.h"
#include "dualnumberarray.h"

// ShadowNumber is subordinate to NumberArray:

#define DualShadowArray_comparisons(templatename) \
template<std::size_t size, typename T, typename T2, typename S, bool reverseorder> \
struct templatename<NumberArray<size, T2>, ShadowNumber<T, S>, reverseorder> { \
  typedef NumberArray<size, typename Symmetric##templatename<T2, ShadowNumber<T, S>, reverseorder>::supertype> supertype; \
}

DualShadowArray_comparisons(CompareTypes);
DualShadowArray_comparisons(PlusType);
DualShadowArray_comparisons(MinusType);
DualShadowArray_comparisons(MultipliesType);
DualShadowArray_comparisons(DividesType);

#endif // __dualshadowarray_h__
