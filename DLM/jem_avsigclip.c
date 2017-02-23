#include <stdio.h>
#include "idl_export.h"

float avsigclip(int, float *, float *, float, float, int);
float weighted_ave(int, float *, float *);

IDL_VPTR jem_avsigclip(int argc, IDL_VPTR *argv)
{
  //variables passed as arguments
  IDL_VPTR data;
  IDL_VPTR weight;
  IDL_VPTR pulllo;
  IDL_VPTR pullhi;
  IDL_VPTR maxiter;

  //C variables to store arguments
  float   *data_d;
  float   *weight_d;
  float   pulllo_d;
  float   pullhi_d;
  IDL_INT maxiter_d;

  //parse input arguments
  data    = argv[0];
  weight  = argv[1];
  pulllo  = argv[2];
  pullhi  = argv[3];
  maxiter = argv[4];

  //some simple checking
  IDL_ENSURE_SIMPLE(data);
  IDL_ENSURE_SIMPLE(weight);
  IDL_ENSURE_SIMPLE(pulllo);
  IDL_ENSURE_SIMPLE(pullhi);
  IDL_ENSURE_SIMPLE(maxiter);
  IDL_ENSURE_ARRAY(data);
  IDL_ENSURE_ARRAY(weight);

  //check types
  if (data->type    != IDL_TYP_FLOAT)
    data    = IDL_CvtFlt(1, argv);
  if (weight->type  != IDL_TYP_FLOAT)
    weight  = IDL_CvtFlt(1, argv+1);
  if (pulllo->type  != IDL_TYP_FLOAT)
    pulllo  = IDL_CvtFlt(1, argv+2);
  if (pullhi->type  != IDL_TYP_FLOAT)
    pullhi  = IDL_CvtFlt(1, argv+3);
  if (maxiter->type != IDL_TYP_INT)
    maxiter = IDL_CvtFix(1, argv+4);
  
  //point to actual array values
  data_d   = (float *) data->value.arr->data;
  weight_d = (float *) weight->value.arr->data;

  //handle non-arrays
  pullhi_d  = pullhi->value.f;
  pulllo_d  = -pulllo->value.f;
  maxiter_d = maxiter->value.i;

  //dimensions of arrays
  IDL_MEMINT *data_dim = data->value.arr->dim;
  IDL_MEMINT out_dim[2];
  out_dim[0] = data_dim[0];
  out_dim[1] = data_dim[1];
  
  //create output array
  IDL_VPTR out;
  float *out_d = (float *) IDL_MakeTempArray(IDL_TYP_FLOAT, 2, out_dim, IDL_ARR_INI_NOP, &out);

  //actual algorithm
  int nx, ny, nexp;
  int ix, iy, iexp;
  float *tempdata, *tempweight; 
  nx = out_dim[0];
  ny = out_dim[1];
  nexp = data_dim[2];
  tempdata = malloc(nexp*sizeof(float));
  tempweight = malloc(nexp*sizeof(float));
  for (ix=0; ix<nx; ix++){
    for (iy=0; iy<ny; iy++){
      for (iexp=0; iexp<nexp; iexp++){
	tempdata[iexp]=data_d[ix+iy*(nx)+iexp*(nx*ny)];
	tempweight[iexp]=weight_d[ix+iy*(nx)+iexp*(nx*ny)];
      }
      out_d[ix+iy*nx]=avsigclip(nexp, tempdata, tempweight, pulllo_d, pullhi_d, maxiter_d);
    }
  }
  
  free(tempdata);
  free(tempweight);

  IDL_DELTMP(data);
  IDL_DELTMP(weight);
  IDL_DELTMP(pulllo);
  IDL_DELTMP(pullhi);
  IDL_DELTMP(maxiter);

  return out;
}

float avsigclip(int size, float *data, float *weight, float pulllo, float pullhi, int maxiter)
{
  int iiter, i;
  int nGood;
  float *dGood, *wGood;
  float pull;
  float ave;

  dGood = malloc(size * sizeof(float));
  wGood = malloc(size * sizeof(float));
  ave = weighted_ave(size, data, weight);

  for (iiter=0; iiter<maxiter; iiter++){
    nGood = 0;
    for (i=0; i<size; i++) {
      pull = (data[i]-ave)*sqrt(weight[i]);
      if (pull > pulllo && pull < pullhi) {
	dGood[nGood] = data[i];
	wGood[nGood++] = weight[i];
	//printf("good!\n");
      } //else printf("rejected!\n");
    }
    if (nGood == 0) {
      iiter = maxiter;
    } else {
      ave = weighted_ave(nGood, dGood, wGood);
    }
  }
  free(dGood);
  free(wGood);
  return ave;
}

float weighted_ave(int size, float *array, float *weight)
{
  int i;
  float total=0.;
  float weighttotal=0.;
  for(i=0; i<size; i++){
    total += array[i]*weight[i];
    weighttotal += weight[i];
  }
  return total/weighttotal;
}

int IDL_Load(void)
{
  static IDL_SYSFUN_DEF2 function_addr[] = {
    { jem_avsigclip, "JEM_AVSIGCLIP", 5, 5, 0, 0}, 
  };
  return IDL_SysRtnAdd(function_addr, TRUE, 
		       IDL_CARRAY_ELTS(function_addr));
}
