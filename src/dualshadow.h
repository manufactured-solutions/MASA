#ifndef __dualshadow_h__
#define __dualshadow_h__

// Order of declarations is important here?
#include "shadownumber.h"
#include "dualnumber.h"

// ShadowNumber is subordinate to DualNumber:

#define DualShadow_comparisons(templatename) \
template<typename T, typename D, typename T2, typename S, bool reverseorder> \
struct templatename<DualNumber<T, D>, ShadowNumber<T2, S>, reverseorder> { \
  typedef DualNumber<typename Symmetric##templatename<T, ShadowNumber<T2, S>, reverseorder>::supertype, \
                     typename Symmetric##templatename<D, ShadowNumber<T2, S>, reverseorder>::supertype> supertype; \
}

DualShadow_comparisons(CompareTypes);
DualShadow_comparisons(PlusType);
DualShadow_comparisons(MinusType);
DualShadow_comparisons(MultipliesType);
DualShadow_comparisons(DividesType);

#endif // __dualshadow_h__

