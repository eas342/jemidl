#include <stdio.h>
#include <stdarg.h>
#include "idl_export.h"

IDL_VPTR fzroots2(int argc, IDL_VPTR *argv, char *argk)
{
  typedef struct {
    IDL_KW_RESULT_FIRST_FIELD;
    int force_type;
    IDL_LONG do_double;
    double eps;
    IDL_LONG no_polish;
    IDL_VPTR tc_input;
  } KW_RESULT;
  static IDL_KW_PAR kw_pars[] = {
    {"DOUBLE", IDL_TYP_LONG, 1, 0, 
     IDL_KW_OFFSETOF(force_type), IDL_KW_OFFSETOF(do_double) },
    {"EPS", IDL_TYP_DOUBLE, 1, 0, 0, IDL_KW_OFFSETOF(eps) },
    {"NO_POLISH", IDL_TYP_LONG, 1, IDL_KW_ZERO,
     0, IDL_KW_OFFSETOF(no_polish) },
    {"TC_INPUT", 0, 1, IDL_KW_OUT|IDL_KW_ZERO,
     0, IDL_KW_OFFSETOF(tc_input) },
    { NULL }
  };
  
  KW_RESULT kw;
  IDL_VPTR result;
  IDL_VPTR c_raw;
  IDL_VPTR c_tc;
  IDL_MEMINT m;
  void *outdata;
  IDL_ARRAY_DIM dim;
  int rtype;
  static IDL_ALLTYPES zero;

  kw.eps = 2.0e-6;
  (void) IDL_KWProcessByOffset(argc, argv, argk,
			       kw_pars, &c_raw, 1, &kw);
  
  IDL_ENSURE_ARRAY(c_raw);
  IDL_ENSURE_SIMPLE(c_raw);
  if (c_raw->value.arr->n_dim != 1)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"Input argument must be a column vector.");
  m = c_raw->value.arr->dim[0];
  if (--m <= 0)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
		"Input array does not have enough elements");
  if (kw.tc_input)
    IDL_StoreScalar(kw.tc_input, IDL_TYP_LONG, &zero);

  if (kw.force_type) {
    rtype = kw.do_double ? IDL_TYP_DCOMPLEX : IDL_TYP_COMPLEX;
  } else {
    rtype = ((c_raw->type == IDL_TYP_DOUBLE)
	     || (c_raw->type == IDL_TYP_DCOMPLEX))
      ? IDL_TYP_DCOMPLEX : IDL_TYP_COMPLEX;
  }
  dim[0] = m;
  outdata = (void *) 
    IDL_MakeTempArray(rtype, 1, dim, IDL_ARR_INI_ZERO, &result);
  
  if (c_raw->type == rtype) {
    c_tc = c_raw;
  } else {
    c_tc = IDL_BasicTypeConversion(1, &c_raw, rtype);
  }
  
  // here would be the zroots call...

  if (kw.tc_input) IDL_VarCopy(c_tc, kw.tc_input);
  else if (c_raw != c_tc) IDL_Deltmp(c_tc);

  IDL_KW_FREE;
  return result;
}
