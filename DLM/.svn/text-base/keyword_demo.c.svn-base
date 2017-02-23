#include <stdio.h>
#include "idl_export.h"

void keyword_demo(int argc, IDL_VPTR *argv, char *argk) 
{
  typedef struct {
    IDL_KW_RESULT_FIRST_FIELD;
    IDL_LONG arr_data[10];
    int arr_there;
    IDL_MEMINT arr_n;
  } KW_RESULT;
  static IDL_KW_ARR_DESC_R arr_d = { IDL_KW_OFFSETOF(arr_data), 3, 10,
				     IDL_KW_OFFSETOF(arr_n) };
  static IDL_KW_PAR kw_pars[] = {
    {"ARRAY", IDL_TYP_LONG, 1, IDL_KW_ARRAY, IDL_KW_OFFSETOF(arr_there), IDL_CHARA(arr_d) },
    { NULL }
    };

  KW_RESULT kw;
  
  (void) IDL_KWProcessByOffset(argc, argv, argk, kw_pars, 
			       (IDL_VPTR *) 0, 1, &kw);

  if (kw.arr_there) 
    printf("banana\n");
  
  IDL_KW_FREE;
}
