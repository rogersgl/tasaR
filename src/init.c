#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP pandaseq_merge_tokens(SEXP args_in);

static const R_CallMethodDef CallEntries[] = {
    {"pandaseq_merge_tokens", (DL_FUNC) &pandaseq_merge_tokens, 1},
    {NULL, NULL, 0}
};

void R_init_tasaR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}