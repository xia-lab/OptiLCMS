#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
// extern void continuousPtsAboveThreshold(void *, void *, void *, void *, void *, void *);
// extern void continuousPtsAboveThresholdIdx(void *, void *, void *, void *, void *, void *);
// extern void DescendMin(void *, void *, void *, void *, void *);
// extern void FindEqualGreaterM(void *, void *, void *, void *, void *);
// extern void RectUnique(void *, void *, void *, void *, void *, void *, void *);
// extern void WhichColMax(void *, void *, void *, void *);
// extern void DescendZero(void *, void *, void *, void *, void *);
// extern void ColMax(void *, void *, void *, void *);
// 
// 
// static const R_CMethodDef CEntries[] = {
//   {"continuousPtsAboveThreshold",          (DL_FUNC) &continuousPtsAboveThreshold,           6},
//   {"continuousPtsAboveThresholdIdx",          (DL_FUNC) &continuousPtsAboveThresholdIdx,           6},
//   {"DescendMin",          (DL_FUNC) &DescendMin,           5},
//   {"FindEqualGreaterM",          (DL_FUNC) &FindEqualGreaterM,           5},
//   {"RectUnique",          (DL_FUNC) &RectUnique,           7},
//   {"WhichColMax",          (DL_FUNC) &WhichColMax,           4},
//   {"DescendZero",          (DL_FUNC) &DescendZero,           5},
//   {"ColMax",          (DL_FUNC) &ColMax,           4},
//   {NULL, NULL, 0}
// };


// void R_init_OptiLCMS(DllInfo *info) {
//   R_RegisterCCallable("OptiLCMS", "continuousPtsAboveThreshold", (DL_FUNC) &continuousPtsAboveThreshold);
//   R_RegisterCCallable("OptiLCMS", "DescendZero", (DL_FUNC) &DescendZero);
// }

//void R_init_OptiLCMS(DllInfo *info) {
//  R_RegisterCCallable("OptiLCMS", "add",  (DL_FUNC) &CEntries);
//}