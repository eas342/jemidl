#include <stdio.h>
#include "idl_export.h"

#define ABS(x) (((x) >=0) ? (x) : -(x))

double getpeak(double array[], int NX, int NY, int *ipk, int *jpk);

IDL_VPTR clean(int argc, IDL_VPTR *argv)
{
  //variables passed as arguments
  IDL_VPTR image;
  IDL_VPTR psf;
  IDL_VPTR loopgain;
  IDL_VPTR niter;

  //C variables to store arguments
  double *image_d;
  double *psf_d;
  double loopgain_d;
  IDL_LONG niter_d;

  //get input arguments
  image = argv[0];
  psf = argv[1];
  loopgain = argv[2];
  niter = argv[3];

  //some simple checking
  IDL_ENSURE_SIMPLE(image);
  IDL_ENSURE_SIMPLE(psf);
  IDL_ENSURE_SIMPLE(loopgain);
  IDL_ENSURE_SIMPLE(niter);
  IDL_ENSURE_ARRAY(image);
  IDL_ENSURE_ARRAY(psf);

  //convert to type double
  if (image->type != IDL_TYP_DOUBLE)
    image = IDL_CvtDbl(1,argv);
  if (psf->type != IDL_TYP_DOUBLE)
    psf = IDL_CvtDbl(1,argv+1);
  if (loopgain->type != IDL_TYP_DOUBLE)
    loopgain = IDL_CvtDbl(1,argv+2);
  if (niter->type != IDL_TYP_LONG)
    niter = IDL_CvtLng(1,argv+3);

  //point to actual array values
  image_d = (double *) image->value.arr->data;
  psf_d   = (double *) psf->value.arr->data;
  //grab scalar values
  loopgain_d = loopgain->value.d;
  niter_d = niter->value.l;

  //dimensions of arrays
  IDL_MEMINT *image_dim = image->value.arr->dim;
  IDL_MEMINT *psf_dim   = psf->value.arr->dim;
  IDL_MEMINT *outimage_dim = (IDL_MEMINT *) malloc(IDL_MAX_ARRAY_DIM*sizeof(IDL_MEMINT));
  outimage_dim[0] = image_dim[0];
  outimage_dim[1] = image_dim[1];
  outimage_dim[2] = 2;

  //create output array
  IDL_VPTR outimage;
  double *outimage_d = (double *) IDL_MakeTempArray(IDL_TYP_DOUBLE, 3, outimage_dim, IDL_ARR_INI_NOP, &outimage);

  //actual code here...

  long iter;
  int i,j;
  int NX,NY,NXp,NYp,ipk,jpk,iii,jjj,ipsf0,jpsf0;
  double *workimage,*gaussimage;
  double vmax,sflux;

  NX = image_dim[0];
  NY = image_dim[1];
  NXp = psf_dim[0];
  NYp = psf_dim[1];

  //  printf("NX:%i NY:%i NXp:%i NYp:%i\n", NX, NY, NXp, NYp);

  if ((workimage = (double *)malloc(NX*NY*sizeof(double))) == NULL )
    {fprintf(stderr,"Not enough mem for image.\n"); exit(1);}
  if ((gaussimage = (double *)malloc(NX*NY*sizeof(double))) == NULL )
    {fprintf(stderr,"Not enough mem for image.\n"); exit(1);}

  //workimage will be the original image and then the residual after each iteration...
  for(i=0;i<NX*NY;i++){
    workimage[i]=image_d[i];
    gaussimage[i]=0.;
  }

  double psfmax=getpeak(psf_d,NXp,NYp,&ipsf0,&jpsf0);


  //scale the psf to total(psf)=1.
  double psftotal=0.;
  for(i=0;i<NXp*NYp;i++)
    psftotal += psf_d[i];
  for(i=0;i<NXp*NYp;i++)
    psf_d[i] /= psftotal;

  //main CLEAN loop
  for (iter=0; iter<niter_d; iter++) {
    vmax = getpeak(workimage,NX,NY,&ipk,&jpk);
    sflux = loopgain_d*vmax;
    
    //subtract scaled PSF and add to gaussimage
    for(j=0; j<NYp; j++) {
      for(i=0; i<NXp; i++) {
	iii=ipk-ipsf0+i;
	jjj=jpk-jpsf0+j;
	if (jjj<0 || jjj>=NY || iii<0 || iii>=NX)
	  continue;
	workimage[jjj*NX + iii] -= (sflux*psf_d[j*NXp+i]);
      }
    }
    gaussimage[jpk*NX + ipk] += sflux;
  }

  //return gaussimage and workimage
  for (i=0;i<NX*NY;i++)
    outimage_d[i] = gaussimage[i];
  for (i=NX*NY; i<2*NX*NY;i++)
    outimage_d[i] = workimage[i-NX*NY];

  free(workimage);
  free(gaussimage);
  free(outimage_dim);

  //clean up and return results.
  IDL_DELTMP(image);
  IDL_DELTMP(psf);
  IDL_DELTMP(loopgain);
  IDL_DELTMP(niter);

  return(outimage);
}

double getpeak(double array[], int NX, int NY, int *ipk, int *jpk)
{
  double vmax, zz;
  int i,j;

  vmax=0.0;
  for(j=0;j<NY;j++) {
    for(i=0;i<NX;i++) {
      if((zz = array[j*NX+i]) > vmax) {
	vmax = zz;
	*ipk = i;
	*jpk = j;
      }
    }
  }
  return vmax;
}

int IDL_Load(void)
{
  static IDL_SYSFUN_DEF2 function_addr[] = {
    { clean, "CLEAN3", 4, 4, 0, 0}, 
  };
  return IDL_SysRtnAdd(function_addr, TRUE, 
		       IDL_CARRAY_ELTS(function_addr));
}
