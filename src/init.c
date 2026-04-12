#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* .C calls */
extern void compute_similarity(void *, void *, void *, void *, void *, void *, void *);
extern void regForest_pattern(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void regForest_predict(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void regForest_represent(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void regForest_similarity(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void regRF_time_series(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"compute_similarity",   (DL_FUNC) &compute_similarity,    7},
    {"regForest_pattern",    (DL_FUNC) &regForest_pattern,    21},
    {"regForest_predict",    (DL_FUNC) &regForest_predict,    21},
    {"regForest_represent",  (DL_FUNC) &regForest_represent,  19},
    {"regForest_similarity", (DL_FUNC) &regForest_similarity, 20},
    {"regRF_time_series",    (DL_FUNC) &regRF_time_series,    32},
    {NULL, NULL, 0}
};

void R_init_LPStimeSeries(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
