#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "clue.h"

static R_NativePrimitiveArgType solve_LSAP_t[3] = { REALSXP, INTSXP, INTSXP };
static R_NativePrimitiveArgType clue_ultrametric_penalty_function_t[3] = { REALSXP, INTSXP, REALSXP };
static R_NativePrimitiveArgType clue_ultrametric_penalty_gradient_t[3] = { REALSXP, INTSXP, REALSXP };
static R_NativePrimitiveArgType clue_dissimilarity_count_inversions_t[4] = { REALSXP, REALSXP, INTSXP, REALSXP };

#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}

static const R_CMethodDef cMethods[] = {
    CDEF(solve_LSAP),
    CDEF(clue_ultrametric_penalty_function),
    CDEF(clue_ultrametric_penalty_gradient),
    CDEF(clue_dissimilarity_count_inversions),
    {NULL, NULL, 0}
};

void R_init_clue(DllInfo *dll)
{
    R_useDynamicSymbols(dll, FALSE);
    R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
}
