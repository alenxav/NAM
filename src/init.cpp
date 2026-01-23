#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Forward declaration of the function from pegs.cpp
extern "C" SEXP PEGS_lapack(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// Define the table of functions available to .Call()
static const R_CallMethodDef CallEntries[] = {
  // { "name_in_R", (DL_FUNC) &function_pointer, number_of_arguments }
  {"PEGS_lapack", (DL_FUNC) &PEGS_lapack, 8},
  
  // This terminator is required
  {NULL, NULL, 0}
};

// This function is called by R when the package is loaded.
// The name MUST be R_init_packagename
void R_init_pegs(DllInfo *dll) {
  // Register the C++ functions with R
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  
  // Disable symbol search, as we have registered all routines
  R_useDynamicSymbols(dll, FALSE);
}