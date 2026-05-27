/* init.c — R routine registration for the ChainBinomial package */
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* Forward declarations */
SEXP r_gen_population(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP r_transtat(SEXP, SEXP, SEXP);
SEXP r_simulate_single(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP r_estimate_single(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallMethods[] = {
    {"r_gen_population",  (DL_FUNC)&r_gen_population,   8},
    {"r_transtat",        (DL_FUNC)&r_transtat,          3},
    {"r_simulate_single", (DL_FUNC)&r_simulate_single, 10},
    {"r_estimate_single", (DL_FUNC)&r_estimate_single, 11},
    {NULL, NULL, 0}
};

void R_init_ChainBinomial(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallMethods, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}
