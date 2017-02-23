#include <stdio.h>
#include "idl_export.h"

IDL_VPTR drizzle1d(int argc, IDL_VPTR *argv)
{

  IDL_VPTR inbins;  //input
  IDL_VPTR influx;  //input
  IDL_VPTR inivar;  //input
  IDL_VPTR outbins; //input
  IDL_VPTR preflux; //input&output
  IDL_VPTR preivar; //input&output

  //data vectors
  double *influx_d;
  double *inbins_d;
  double *inivar_d;
  double *outbins_d;
  double *preflux_d;
  double *preivar_d;
  //n_elements variables
  long influx_n;
  long inbins_n;
  long inivar_n;
  long outbins_n;
  long preflux_n;
  long preivar_n;

  //grab positional arguments
  inbins = argv[0];
  influx = argv[1];
  inivar = argv[2];
  outbins = argv[3];
  preflux = argv[4];
  preivar = argv[5];
  
  //some error checking
  IDL_ENSURE_SIMPLE(inbins);
  IDL_ENSURE_SIMPLE(influx); 
  IDL_ENSURE_SIMPLE(inivar);
  IDL_ENSURE_SIMPLE(outbins);
  IDL_ENSURE_SIMPLE(preflux);
  IDL_ENSURE_SIMPLE(preivar);
  IDL_ENSURE_ARRAY(inbins);
  IDL_ENSURE_ARRAY(influx); 
  IDL_ENSURE_ARRAY(inivar);
  IDL_ENSURE_ARRAY(outbins);
  IDL_ENSURE_ARRAY(preflux);
  IDL_ENSURE_ARRAY(preivar);
  
  //grab data vectors...
  inbins_d = (double *) inbins->value.arr->data;
  influx_d = (double *) influx->value.arr->data;  
  inivar_d = (double *) inivar->value.arr->data;  
  outbins_d = (double *) outbins->value.arr->data;
  preflux_d = (double *) preflux->value.arr->data;  
  preivar_d = (double *) preivar->value.arr->data;  
  //and lengths...
  inbins_n = *inbins->value.arr->dim;
  influx_n = *influx->value.arr->dim;
  inivar_n = *inivar->value.arr->dim;
  outbins_n = *outbins->value.arr->dim;
  preflux_n = *preflux->value.arr->dim;
  preivar_n = *preivar->value.arr->dim;

  IDL_MEMINT out_n = preflux_n+preivar_n;
  
  IDL_VPTR out;
  double *outflux_d = (double *) 
    IDL_MakeTempArray(IDL_TYP_DOUBLE, 1, &out_n, IDL_ARR_INI_NOP, &out);
  double *outivar_d = outflux_d+preflux_n;

  int i;
  for (i=0; i<preflux_n;i++){
    outflux_d[i] = preflux_d[i];
    outivar_d[i] = preivar_d[i];
  }
  
  //actual algorithm should go here...
  long inbin, outbin;
  for(inbin=0; inbin < influx_n; inbin++) {
    for(outbin=0; outbin < outbins_n-1; outbin++ ) {
      if (outbins_d[outbin] >= inbins_d[inbin+1]) break;
      if (inbins_d[inbin] >= outbins_d[outbin+1]) continue;
      //case 1: input bin overhangs the output bin on the left
      if (inbins_d[inbin] <= outbins_d[outbin]) {
	double fraction = (inbins_d[inbin+1]-outbins_d[outbin])/
	  (inbins_d[inbin+1]-inbins_d[inbin]);
	double newbinivar = fraction * inivar_d[inbin] + outivar_d[outbin];
	if (newbinivar == 0.) continue;
	outflux_d[outbin] = (influx_d[inbin]*fraction*inivar_d[inbin]
			     + outflux_d[outbin]*outivar_d[outbin]) / newbinivar;
	outivar_d[outbin] = newbinivar;
      } 
      //case 2: input bin overhangs the output bin on the right
      else if (inbins_d[inbin+1] >= outbins_d[outbin+1]) {
	double fraction = (outbins_d[outbin+1]-inbins_d[inbin])/
	  (inbins_d[inbin+1] - inbins_d[inbin]);
	double newbinivar = fraction*inivar_d[inbin]+outivar_d[outbin];
	if (newbinivar == 0.) continue;
	outflux_d[outbin] = (influx_d[inbin]*fraction*inivar_d[inbin]
			     + outflux_d[outbin]*outivar_d[outbin]) / newbinivar;
	outivar_d[outbin] = newbinivar;
      } 
      //case 3: input bin is totally encompassed by the output bin
      else {
	double newbinivar = inivar_d[inbin]+outivar_d[outbin];
	if (newbinivar == 0.) continue;
	outflux_d[outbin] = (influx_d[inbin]*inivar_d[inbin]
			     + outflux_d[outbin]*outivar_d[outbin])/newbinivar;
	outivar_d[outbin] = newbinivar;
      }
    }
  }
  return(out);
}

int IDL_Load(void)
{
  static IDL_SYSFUN_DEF2 function_addr[] = {
    { drizzle1d, "DRIZZLE1D", 6, 6, 0, 0},
  };
  return IDL_SysRtnAdd(function_addr, TRUE,
		       IDL_CARRAY_ELTS(function_addr));
}
