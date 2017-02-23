#include <stdio.h>
#include "idl_export.h"

#define ABS(x) (((x) >=0) ? (x) : -(x))
#define defGAIN 0.003
#define defSIGCUT 5.0
#define defTHRESHCUT 7.5
#define defNITER 500000
#define defPRINTITER 1000

double getpeak(double array[], int NX, int NY, int *ipk, int *jpk);
void avevar(double array[], int size, double *ave, double *var);

IDL_VPTR clean(int argc, IDL_VPTR *argv)
{
  //variables passed as arguments
  IDL_VPTR image;
  IDL_VPTR psf;
  IDL_VPTR gauss;
  IDL_VPTR thresh;

  //C variables to store arguments
  double *image_d;
  double *psf_d;
  double *gausspsf_d;
  double gauss_d;
  double thresh_d;

  //get input arguments
  image = argv[0];
  psf = argv[1];
  gauss = argv[2];
  thresh = argv[3];

  //some simple checking
  IDL_ENSURE_SIMPLE(image);
  IDL_ENSURE_SIMPLE(psf);
  IDL_ENSURE_SIMPLE(gauss);
  IDL_ENSURE_SIMPLE(thresh);
  IDL_ENSURE_ARRAY(image);
  IDL_ENSURE_ARRAY(psf);

  //convert to type double
  if (image->type != IDL_TYP_DOUBLE)
    image = IDL_CvtDbl(1,argv);
  if (psf->type != IDL_TYP_DOUBLE)
    psf = IDL_CvtDbl(1,argv+1);
  if (gauss->type != IDL_TYP_DOUBLE)
    gauss = IDL_CvtDbl(1,argv+2);
  if (thresh->type != IDL_TYP_DOUBLE)
    thresh = IDL_CvtDbl(1,argv+3);

  //point to actual array values
  image_d = (double *) image->value.arr->data;
  psf_d   = (double *) psf->value.arr->data;
  //size of gaussian with which to convolve CLEAN map before adding residual...
  gauss_d = gauss->value.d;
  thresh_d = thresh->value.d;

  //dimensions of arrays
  IDL_MEMINT *image_dim = image->value.arr->dim;
  IDL_MEMINT *psf_dim   = psf->value.arr->dim;

  //create output array
  IDL_VPTR outimage;
  double *outimage_d = (double *) IDL_MakeTempArray(IDL_TYP_DOUBLE, 2, image_dim, IDL_ARR_INI_NOP, &outimage);

  //actual code here...

  long niter=defNITER, iter;
  int i,j;
  int NX,NY,NXp,NYp,ipk,jpk,iii,jjj,ipsf0,jpsf0;
  double *workimage,*psfimage,*gaussimage;
  double gain=defGAIN;
  double vmax,sflux;
  double ave, var;

  NX = image_dim[0];
  NY = image_dim[1];
  NXp = psf_dim[0];
  NYp = psf_dim[1];

  if ((workimage = (double *)malloc(NX*NY*sizeof(double))) == NULL )
    {fprintf(stderr,"Not enough mem for image.\n"); exit(1);}
  if ((psfimage = (double *)malloc(NX*NY*sizeof(double))) == NULL )
    {fprintf(stderr,"Not enough mem for image.\n"); exit(1);}
  if ((gaussimage = (double *)malloc(NX*NY*sizeof(double))) == NULL )
    {fprintf(stderr,"Not enough mem for image.\n"); exit(1);}
  if ((gausspsf_d = (double*)malloc(NXp*NYp*sizeof(double))) == NULL )
    {fprintf(stderr,"Not enough mem for guasspsf_d.\n"); exit(1);}

  //workimage will be the original image and then the residual after each iteration...
  for(i=0;i<NX*NY;i++){
    workimage[i]=image_d[i];
    psfimage[i]=0.;
    gaussimage[i]=0.;
  }

  double psfmax=getpeak(psf_d,NXp,NYp,&ipsf0,&jpsf0);

  //create gausspsf_d
  double sigma = gauss_d/2.35482;
  for(i=0;i<NXp;i++){
    for(j=0;j<NYp;j++){
      //now average over this pixel...
      int ii=0,jj=0;
      gausspsf_d[j*NXp+i]=0.;
      for(ii=-2;ii<=2;ii++) {
	for(jj=-2;jj<=2;jj++) {
	  double x=i-ipsf0+(ii*0.2);
	  double y=j-jpsf0+(jj*0.2);
	  double r2=x*x+y*y;
	  gausspsf_d[j*NXp+i] += exp(-0.5*r2/sigma/sigma);
	}
      }
    }
  }

  //scale the psf to total(psf)=1.
  double psftotal=0.;
  double gausstotal=0.;
  for(i=0;i<NXp*NYp;i++){
    psftotal += psf_d[i];
    gausstotal += gausspsf_d[i];
  }
  for(i=0;i<NXp*NYp;i++){
    psf_d[i] /= psftotal;
    gausspsf_d[i] /= gausstotal;
  }

  //main CLEAN loop
  for (iter=0; iter<niter; iter++) {
    vmax = getpeak(workimage,NX,NY,&ipk,&jpk);
    sflux = gain*vmax;
    
    //subtract scaled PSF and add to gaussimage
    for(j=0; j<NYp; j++) {
      for(i=0; i<NXp; i++) {
	iii=ipk-ipsf0+i;
	jjj=jpk-jpsf0+j;
	if (jjj<0 || jjj>=NY || iii<0 || iii>=NX)
	  continue;
	workimage[jjj*NX + iii] -= (sflux*psf_d[j*NXp+i]);
	gaussimage[jjj*NX + iii] += (sflux*gausspsf_d[j*NXp+i]);
      }
    }
    if (iter%defPRINTITER==0) {
      avevar(workimage,NX*NY,&ave,&var);
      //printf("%i %6.3f %6.3f\n",iter, vmax/sqrt(var), vmax/thresh_d);
      if (vmax<sqrt(var)*defSIGCUT && vmax<thresh_d*defTHRESHCUT) break;
    }
  }

  //ADD BACK RESIDUAL
  for (i=0;i<NX*NY;i++) {
    outimage_d[i] = workimage[i]+gaussimage[i];
  }

  free(workimage);
  free(psfimage);
  free(gaussimage);
  free(gausspsf_d);

  //clean up and return results.
  IDL_DELTMP(image);
  IDL_DELTMP(psf);
  IDL_DELTMP(gauss);

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

void avevar(double array[], int size, double *ave, double *var)
{
  double s, ep;
  int j;
  *ave=0.0;
  for (j=0;j<size;j++) *ave += array[j];
  *ave /= size;
  *var=ep=0.0;
  for (j=0;j<size;j++){
    s=array[j]-*ave;
    ep += s;
    *var += s*s;
  }
  *var=(*var-ep*ep/size)/(size-1);
}


int IDL_Load(void)
{
  static IDL_SYSFUN_DEF2 function_addr[] = {
    { clean, "CLEAN2", 4, 4, 0, 0}, 
  };
  return IDL_SysRtnAdd(function_addr, TRUE, 
		       IDL_CARRAY_ELTS(function_addr));
}
