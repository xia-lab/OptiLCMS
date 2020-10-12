#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void continuousPtsAboveThreshold(void *, void *, void *, void *, void *, void *);
extern void continuousPtsAboveThresholdIdx(void *, void *, void *, void *, void *, void *);
extern void DescendMin(void *, void *, void *, void *, void *);
extern void FindEqualGreaterM(void *, void *, void *, void *, void *);
extern void RectUnique(void *, void *, void *, void *, void *, void *, void *);
extern void WhichColMax(void *, void *, void *, void *);
extern void DescendZero(void *, void *, void *, void *, void *);
extern void ColMax(void *, void *, void *, void *);


/* .Call calls */
extern SEXP findmzROI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP getMZ(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP getEIC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP binYonX(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP binYonX_multi(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impute_with_linear_interpolation(SEXP, SEXP);
extern SEXP impute_with_linear_interpolation_base(SEXP, SEXP, SEXP);
extern SEXP breaks_on_nBins(SEXP, SEXP, SEXP, SEXP);
extern SEXP R_set_obiwarp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
  {"continuousPtsAboveThreshold",          (DL_FUNC) &continuousPtsAboveThreshold,           6},
  {"continuousPtsAboveThresholdIdx",          (DL_FUNC) &continuousPtsAboveThresholdIdx,           6},
  {"DescendMin",          (DL_FUNC) &DescendMin,           5},
  {"FindEqualGreaterM",          (DL_FUNC) &FindEqualGreaterM,           5},
  {"RectUnique",          (DL_FUNC) &RectUnique,           7},
  {"WhichColMax",          (DL_FUNC) &WhichColMax,           4},
  {"DescendZero",          (DL_FUNC) &DescendZero,           5},
  {"ColMax",          (DL_FUNC) &ColMax,           4},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"findmzROI",  (DL_FUNC) &findmzROI,  10},
    {"getMZ",  (DL_FUNC) &getMZ,  6},
    {"getEIC",  (DL_FUNC) &getEIC,  6},
    {"binYonX",  (DL_FUNC) &binYonX,  14},
    {"binYonX_multi",  (DL_FUNC) &binYonX_multi,  14},
    {"impute_with_linear_interpolation",  (DL_FUNC) &impute_with_linear_interpolation,  2},
    {"impute_with_linear_interpolation_base",  (DL_FUNC) &impute_with_linear_interpolation_base,  3},
    {"breaks_on_nBins",  (DL_FUNC) &breaks_on_nBins,  4},
    {"R_set_obiwarp", (DL_FUNC) &R_set_obiwarp, 18},
    {NULL, NULL, 0}
};

 void R_init_OptiLCMS(DllInfo *dll)
 {
   R_registerRoutines(dll, CEntries, CallEntries,NULL,NULL);
   R_useDynamicSymbols(dll, FALSE);
 }

//void R_init_OptiLCMS(DllInfo *info) {
//  R_RegisterCCallable("OptiLCMS", "add",  (DL_FUNC) &CEntries);
//}
