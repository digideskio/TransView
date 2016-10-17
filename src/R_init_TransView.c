
#include <R_ext/Rdynload.h>
#include "slice_dc.h"
#include "construct_dc.h"

static const R_CallMethodDef callMethods[] =
{
	{"slice_dc", (DL_FUNC)&slice_dc, 7},
    {"construct_dc", (DL_FUNC)&construct_dc, 3},
    {"approx_window", (DL_FUNC)&approx_window, 3},
    {NULL,NULL, 0}
};

void R_init_construct_dc(DllInfo *dll)
{
    R_registerRoutines(dll,NULL,callMethods,NULL,NULL);
}
