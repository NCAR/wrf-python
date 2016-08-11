#include <stdio.h>
#include <strings.h>
#include <math.h>
#include "wrapper.h"

#define ERRLEN 512

extern void NGCALLF(dcomputetk,DCOMPUTETK)(double *,double *,double *,int *);
extern void NGCALLF(dcomputetd,DCOMPUTETD)(double *,double *,double *,int *);
extern void NGCALLF(dcomputerh,DCOMPUTERH)(double *,double *,double *,
                                           double *,int *);
extern void NGCALLF(dcomputeseaprs,DCOMPUTESEAPRS)(int *,int *,int *,
                                                   double *,double *,
                                                   double *,double *,
                                                   double *,double *,
                                                   double *,double *,
												   int *, char *, int);

extern void NGCALLF(dinterp3dz,DINTERP3DZ)(double *,double *,double *,
                                           double *,int *,int *, int*,
                                           double *);

extern void NGCALLF(dinterp2dxy,DINTERP2DXY)(double *,double *,double *,
                                             int *,int *,int *, int*);

extern void NGCALLF(dinterp1d,DINTERP1D)(double *,double *,double *,double *,
                                         int *, int *, double *);

extern void NGCALLF(dfilter2d,DFILTER2D)(double *, double *, int *, int *, 
                                         int *, double *);

extern void NGCALLF(filter2d,FILTER2D)(float *, float *, int *, int *, 
                                       int *, float *);

extern void NGCALLF(dgetijlatlong,DGETIJLATLONG)(double *, double *, double *,
                                                 double *, int *, int *,
                                                 int *, int *, int *);


extern void NGCALLF(dcomputeuvmet,DCOMPUTEUVMET)(double *, double *, double *,
                                                 double *, double *, double *,
                                                 double *, double *, double *,
                                                 double *, int *, int *, 
                                                 int *, int *, int *,
                                                 logical *,double *, double*,
                                                 double *);

extern void NGCALLF(dcomputeiclw,DCOMPUTEICLW)(double *, double *, double *, 
                                               int *, int *, int *);

extern void NGCALLF(dbint3d,DBINT3D)(double *,double *,double *, double *,
                                     int *, int *, int *, int *,
                                     int *, int *, int *);


extern void NGCALLF(dcomputepv,DCOMPUTEPV)(double *, double *, double *, 
                                           double *, double *, double *, 
                                           double *, double *, double *, 
                                           double *, double *, int *, int *, 
                                           int *, int *, int *);

extern void NGCALLF(dcomputeabsvort,DCOMPUTEABSVORT)(double *, double *,
                                                     double *, double *,
                                                     double *, double *,
                                                     double *, double *,
                                                     double *, int *, int *,
                                                     int *, int *, int *);

extern void NGCALLF(calcdbz,CALCDBZ)(double *, double *, double *, double *,
                                     double *, double *, double *, int *, 
                                     int *, int *, int *, int *, int *);

extern void NGCALLF(dlltoij,DLLTOIJ)(int *, double *, double *, double *, 
                                     double *, double *, double *, double *, 
                                     double *, double *, double *, double *, 
                                     double *, double *, double *, double *, 
                                     double *, int *, char *, int);

extern void NGCALLF(dijtoll,DIJTOLL)(int *, double *, double *, double *, 
                                     double *, double *, double *, double *, 
                                     double *, double *, double *, double *, 
                                     double *, double *, double *, double *, 
                                     double *, int *, char *, int);

extern void NGCALLF(deqthecalc,DEQTHECALC)(double *, double *, double *,
                                           double *, int *, int *, int *);

extern void NGCALLF(dcapecalc3d,DCAPECALC3D)(double *prs, double *tmk, 
                                             double *qvp, double *ght,
                                             double *ter, double *sfp, 
                                             double *cape, double *cin, 
                                             double *cmsg,
                                             int *miy, int *mjx, int *mkzh, 
                                             int *i3dflag, int *ter_follow,
                                             char *psafile,
											 int *errstat, char *errmsg,
											 int, int);

extern void NGCALLF(dcalrelhl,DCALRELHL)(double *u, double *v, double *z,
                                         double *ter, double *top, 
                                         double *sreh, int *miy, int *mjx,
                                         int *mkzh);

extern void NGCALLF(dcalcuh,DCALCUH)(int *, int *, int *, int *, double *, 
                                     double *, double *, double *, double *,
                                     double *, double *, double *, double *,
                                     double *, double *, double *);

extern void NGCALLF(plotgrids_var,PLOTGRIDS_VAR)(char *fname, float *plotvar, int);

extern void NGCALLF(plotfmt_open,PLOTFMT_OPEN)(char *cfilename, int *istatus,
                                               int);

extern void NGCALLF(plotfmt_close,PLOTFMT_CLOSE)();

extern void NGCALLF(plotfmt_rdhead,PLOTFMT_RDHEAD)(int *istatus,
                                                   float *rhead,
                                                   char *cfield, 
                                                   char *chdate,
                                                   char *cunits,
                                                   char *cmapsc,
                                                   char *cdesc,
                                                   int,int,int,int,int);

extern void NGCALLF(plotfmt_rddata,PLOTFMT_RDDATA)(int *istatus,
                                                   int *nx, int *ny,
                                                   float *slab);

extern void NGCALLF(wetbulbcalc,WETBULBCALC)(double *prs,double *tmk, 
                                             double *qvp, double *twb,
                                             int *nx, int *ny,int *nz,
					     char *,int);

extern void NGCALLF(omgcalc,OMGCALC)(double *qvp, double *tmk, double *www,
                                     double *prs, double *omg, 
                                     int *mx,int *my, int *mz);

extern void NGCALLF(virtual_temp,VIRTUAL_TEMP)(double *temp,double *ratmix,
                                               double *tv,int *,int *,int *);

extern NclDimRec *get_wrf_dim_info(int,int,int,ng_size_t*);

extern char *get_psa_file();

extern void var_zero(double *, ng_size_t);

extern void convert_to_hPa(double *, ng_size_t);

extern void flip_it(double *, double *, ng_size_t, ng_size_t);

NhlErrorTypes wrf_tk_W( void )
{
/*
 * Input array variables
 */
  void *p, *theta;
  double *tmp_p = NULL;
  double *tmp_theta = NULL;
  int ndims_p, ndims_theta;
  ng_size_t dsizes_p[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_theta[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_p, type_theta;
/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info;

/*
 * Output variable and attributes.
 */
  void *t;
  NclQuark *description, *units;
  char *cdescription, *cunits;
  double *tmp_t = NULL;
  int size_tt;
  NclBasicDataTypes type_t;
  NclObjClass type_obj_t;
/*
 * Various
 */
  ng_size_t i, nx, size_leftmost, index_p;
  int inx;

/*
 * Variables for returning the output array with attributes attached.
 */
  int att_id;
  ng_size_t dsizes[1];
  NclMultiDValData att_md, return_md;
  NclVar tmp_var;
  NclStackEntry return_data;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
  p = (void*)NclGetArgValue(
           0,
           2,
           &ndims_p,
           dsizes_p,
           NULL,
           NULL,
           &type_p,
           DONT_CARE);

  theta = (void*)NclGetArgValue(
           1,
           2,
           &ndims_theta,
           dsizes_theta,
           NULL,
           NULL,
           &type_theta,
           DONT_CARE);
/*
 * Error checking. Input variables must be same size.
 */
  if(ndims_p != ndims_theta) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_tk: The p and theta arrays must have the same number of dimensions");
    return(NhlFATAL);
  }
  for(i = 0; i < ndims_p; i++) {
    if(dsizes_p[i] != dsizes_theta[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_tk: p and theta must be the same dimensionality");
      return(NhlFATAL);
    }
  }

/*
 * Retrieve dimension names from the "theta" variable, if any.
 * These dimension names will later be attached to the output variable.
 */
  dim_info = get_wrf_dim_info(1,2,ndims_theta,dsizes_theta);

/*
 * Calculate size of leftmost dimensions.
 */
  size_leftmost = 1;
  for(i = 0; i < ndims_p-1; i++) size_leftmost *= dsizes_p[i];
  nx = dsizes_p[ndims_p-1];
  size_tt = size_leftmost * nx;

/*
 * Test dimension sizes.
 */
  if(nx > INT_MAX) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_tk: nx = %ld is greater than INT_MAX", nx);
    return(NhlFATAL);
  }
  inx = (int) nx;

/* 
 * Allocate space for coercing input arrays.  If the input p or theta
 * are already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 *
 * The output type defaults to float, unless any of the two input arrays
 * are double.
 */
  type_t     = NCL_float;
  type_obj_t = nclTypefloatClass;
  if(type_p != NCL_double) {
    tmp_p = (double *)calloc(nx,sizeof(double));
    if(tmp_p == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_tk: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_t     = NCL_double;
    type_obj_t = nclTypedoubleClass;
  }

  if(type_theta != NCL_double) {
    tmp_theta = (double *)calloc(nx,sizeof(double));
    if(tmp_theta == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_tk: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_t     = NCL_double;
    type_obj_t = nclTypedoubleClass;
  }

/*
 * Allocate space for output array.
 */ 
  if(type_t == NCL_double) {
    t = (double *)calloc(size_tt,sizeof(double));
    if(t == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_tk: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }
  else {
    t     = (float *)calloc(size_tt,sizeof(float));
    tmp_t = (double *)calloc(nx,sizeof(double));
    if(tmp_t == NULL || t == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_tk: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }
/*
 * Loop across leftmost dimensions and call the Fortran routine for each
 * one-dimensional subsection.
 */
  index_p = 0;
  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of p (tmp_p) to double if necessary.
 */
    if(type_p != NCL_double) {
      coerce_subset_input_double(p,tmp_p,index_p,type_p,nx,0,NULL,NULL);
    }
    else {
      tmp_p = &((double*)p)[index_p];
    }
/*
 * Coerce subsection of theta (tmp_theta) to double if ncessary.
 */
    if(type_theta != NCL_double) {
      coerce_subset_input_double(theta,tmp_theta,index_p,type_theta,nx,
                                 0,NULL,NULL);
    }
    else {
      tmp_theta = &((double*)theta)[index_p];
    }

/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_t == NCL_double) tmp_t = &((double*)t)[index_p];
/*
 * Call Fortran routine.
 */
      NGCALLF(dcomputetk,DCOMPUTETK)(tmp_t,tmp_p,tmp_theta,&inx);

/*
 * Coerce output back to float if necessary.
 */
    if(type_t == NCL_float) {
      coerce_output_float_only(t,tmp_t,nx,index_p);
    }

    index_p += nx;    /* Increment index */
  }
/*
 * Free up memory.
 */
  if(type_p     != NCL_double) NclFree(tmp_p);
  if(type_theta != NCL_double) NclFree(tmp_theta);
  if(type_t     != NCL_double) NclFree(tmp_t);

/*
 * Set up some attributes ("description" and "units") to return.
 * Note that if the input arrays are anything but 2D, the units
 * will be "Temperature", and "2m Temperature" otherwise.
 */
  if(ndims_p != 2) {
    cdescription = (char *)calloc(12,sizeof(char));
    strcpy(cdescription,"Temperature");
  }
  else {
    cdescription = (char *)calloc(15,sizeof(char));
    strcpy(cdescription,"2m Temperature");
  }
  cunits       = (char *)calloc(2,sizeof(char));
  strcpy(cunits,"K");
  description = (NclQuark*)NclMalloc(sizeof(NclQuark));
  units       = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *description = NrmStringToQuark(cdescription);
  *units       = NrmStringToQuark(cunits);
  free(cunits);
  free(cdescription);

/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)t,
                            NULL,
                            ndims_p,
                            dsizes_p,
                            TEMPORARY,
                            NULL,
                            type_obj_t
                            );
/*
 * Set up attributes to return.
 */
  att_id = _NclAttCreate(NULL,NULL,Ncl_Att,0,NULL);

  dsizes[0] = 1;
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)description,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "description",
             att_md,
             NULL
             );
    
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)units,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "units",
             att_md,
             NULL
             );
    
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          att_id,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );

  NclFree(dim_info);
/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);

}

NhlErrorTypes wrf_td_W( void )
{
/*
 * Input array variables
 */
  void *p, *qv;
  double *tmp_p, *tmp_qv;
  int ndims_p, ndims_qv;
  ng_size_t dsizes_p[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_qv[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_p, type_qv;
/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info;

/*
 * Output variable and attributes.
 */
  void *t;
  NclQuark *description, *units;
  char *cdescription, *cunits;
  double *tmp_t = NULL;
  int size_tt;
  NclBasicDataTypes type_t;
  NclObjClass type_obj_t;
/*
 * Various
 */
  ng_size_t i, np, nx, size_leftmost, index_p;
  int inx;

/*
 * Variables for returning the output array with attributes attached.
 */
  int att_id;
  ng_size_t dsizes[1];
  NclMultiDValData att_md, return_md;
  NclVar tmp_var;
  NclStackEntry return_data;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
  p = (void*)NclGetArgValue(
           0,
           2,
           &ndims_p,
           dsizes_p,
           NULL,
           NULL,
           &type_p,
           DONT_CARE);

  qv = (void*)NclGetArgValue(
           1,
           2,
           &ndims_qv,
           dsizes_qv,
           NULL,
           NULL,
           &type_qv,
           DONT_CARE);

/*
 * Error checking. Input variables must be same size.
 */
  if(ndims_p != ndims_qv) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_td: The p and qv arrays must have the same number of dimensions");
    return(NhlFATAL);
  }
  for(i = 0; i < ndims_p; i++) {
    if(dsizes_p[i] != dsizes_qv[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_td: p and qv must be the same dimensionality");
      return(NhlFATAL);
    }
  }

/*
 * Retrieve dimension names from the "qvapor" variable, if any.
 * These dimension names will later be attached to the output variable.
 */
  dim_info = get_wrf_dim_info(1,2,ndims_qv,dsizes_qv);

/*
 * Calculate size of leftmost dimensions.
 */
  size_leftmost = 1;
  for(i = 0; i < ndims_p-1; i++) size_leftmost *= dsizes_p[i];
  nx = dsizes_p[ndims_p-1];
  size_tt = size_leftmost * nx;

/*
 * Test dimension sizes.
 */
  if(nx > INT_MAX) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_td: nx = %ld is greater than INT_MAX", nx);
    return(NhlFATAL);
  }
  inx = (int) nx;

/* 
 * Allocate space for coercing input arrays.  If the input p or qv
 * are already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 *
 * The output type defaults to float, unless any of the two input arrays
 * are double.
 */
  type_t     = NCL_float;
  type_obj_t = nclTypefloatClass;
/*
 * Allocate space for tmp_p no matter what, because we have to
 * convert the values from hPa to Pa, and we don't want to do
 * this to the original array.
 */
  tmp_p = (double *)calloc(nx,sizeof(double));
  if(tmp_p == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_td: Unable to allocate memory for coercing 't' to double");
    return(NhlFATAL);
  }

  if(type_p == NCL_double) {
    type_t     = NCL_double;
    type_obj_t = nclTypedoubleClass;
  }

/*
 * Allocate space for tmp_qv no matter what, because we want to set
 * values of qv that are less than zero to zero, but we don't want
 * these values retained when the function is done.
 */
  tmp_qv = (double *)malloc(nx*sizeof(double));
  if(tmp_qv == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_td: Unable to allocate memory for coercing 'qv' to double");
    return(NhlFATAL);
  }
  if(type_qv == NCL_double) {
    type_t     = NCL_double;
    type_obj_t = nclTypedoubleClass;
  }

/*
 * Allocate space for output array.
 */ 
  if(type_t == NCL_double) {
    t = (double *)calloc(size_tt,sizeof(double));
    if(t == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_td: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }
  else {
    t     = (float *)calloc(size_tt,sizeof(float));
    tmp_t = (double *)calloc(nx,sizeof(double));
    if(tmp_t == NULL || t == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_td: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }
/*
 * Loop across leftmost dimensions and call the Fortran routine for each
 * one-dimensional subsection.
 */
  index_p = 0;
  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of p (tmp_p) to double if necessary. Otherwise,
 * just do a memcpy. Afterwards, convert the p values to Pa units,
 * as they are coming in as hPa.
 */
    if(type_p != NCL_double) {
      coerce_subset_input_double(p,tmp_p,index_p,type_p,nx,0,NULL,NULL);
    }
    else {
      (void *)memcpy((void*)((char*)tmp_p),
                     (void*)((char*)p + (index_p*sizeof(double))),
                     sizeof(double)*nx);
    }
    for(np = 0; np < nx; np++) tmp_p[np] *= 0.01;
/*
 * Coerce subsection of qv (tmp_qv) to double if ncessary. Otherwise,
 * just do a memcpy. Afterwards, set all values < 0 to 0.
 */
    if(type_qv != NCL_double) {
      coerce_subset_input_double(qv,tmp_qv,index_p,type_qv,nx,
                                 0,NULL,NULL);
    }
    else {
      (void *)memcpy((void*)((char*)tmp_qv),
                     (void*)((char*)qv + (index_p*sizeof(double))),
                     sizeof(double)*nx);
    }
    var_zero(tmp_qv, nx);    /* Set all values < 0 to 0. */

/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_t == NCL_double) tmp_t = &((double*)t)[index_p];
/*
 * Call Fortran routine.
 */
    NGCALLF(dcomputetd,DCOMPUTETD)(tmp_t,tmp_p,tmp_qv,&inx);

/*
 * Coerce output back to float if necessary.
 */
    if(type_t == NCL_float) {
      coerce_output_float_only(t,tmp_t,nx,index_p);
    }

    index_p += nx;    /* Increment index */
  }
/*
 * Free up memory.
 */
  NclFree(tmp_qv);
  NclFree(tmp_p);
  if(type_t  != NCL_double) NclFree(tmp_t);

/*
 * Set up some attributes ("description" and "units") to return.
 * Note that if the input arrays are anything but 2D, the units
 * will be "Temperature", and "2m Temperature" otherwise.
 */
  if(ndims_p != 2) {
    cdescription = (char *)calloc(21,sizeof(char));
    strcpy(cdescription,"Dewpoint Temperature");
  }
  else {
    cdescription = (char *)calloc(24,sizeof(char));
    strcpy(cdescription,"2m Dewpoint Temperature");
  }
  cunits       = (char *)calloc(2,sizeof(char));
  strcpy(cunits,"C");
  description = (NclQuark*)NclMalloc(sizeof(NclQuark));
  units       = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *description = NrmStringToQuark(cdescription);
  *units       = NrmStringToQuark(cunits);
  free(cdescription);
  free(cunits);

/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)t,
                            NULL,
                            ndims_p,
                            dsizes_p,
                            TEMPORARY,
                            NULL,
                            type_obj_t
                            );
/*
 * Set up attributes to return.
 */
  att_id = _NclAttCreate(NULL,NULL,Ncl_Att,0,NULL);

  dsizes[0] = 1;
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)description,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "description",
             att_md,
             NULL
             );
    
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)units,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "units",
             att_md,
             NULL
             );
    

  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          att_id,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );

  NclFree(dim_info);
/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);

}

NhlErrorTypes wrf_rh_W( void )
{
/*
 * Input array variables
 */
  void *qv, *p, *t;
  double *tmp_qv = NULL; 
  double *tmp_p = NULL; 
  double *tmp_t = NULL;
  int ndims_qv, ndims_p, ndims_t;
  ng_size_t dsizes_qv[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_p[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_t[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_qv, type_p, type_t;
/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info;

/*
 * Output variable and attributes.
 */
  void *rh;
  NclQuark *description, *units;
  char *cdescription, *cunits;
  double *tmp_rh = NULL;
  int size_rh;
  NclBasicDataTypes type_rh;
  NclObjClass type_obj_rh;
/*
 * Various
 */
  ng_size_t i, nx, size_leftmost, index_qv;
  int inx;

/*
 * Variables for returning the output array with attributes and/or
 * dimension names attached.
 */
  int att_id;
  ng_size_t dsizes[1];
  NclMultiDValData att_md, return_md;
  NclVar tmp_var;
  NclStackEntry return_data;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
  qv = (void*)NclGetArgValue(
           0,
           3,
           &ndims_qv,
           dsizes_qv,
           NULL,
           NULL,
           &type_qv,
           DONT_CARE);

  p = (void*)NclGetArgValue(
           1,
           3,
           &ndims_p,
           dsizes_p,
           NULL,
           NULL,
           &type_p,
           DONT_CARE);

  t = (void*)NclGetArgValue(
           2,
           3,
           &ndims_t,
           dsizes_t,
           NULL,
           NULL,
           &type_t,
           DONT_CARE);

/*
 * Error checking. Input variables must be same size.
 */
  if(ndims_qv != ndims_t || ndims_p != ndims_t) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_rh: The qv, p, and t arrays must have the same number of dimensions");
    return(NhlFATAL);
  }
  for(i = 0; i < ndims_qv; i++) {
    if(dsizes_qv[i] != dsizes_t[i] || dsizes_p[i] != dsizes_t[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_rh: qv, p, and t must be the same dimensionality");
      return(NhlFATAL);
    }
  }

/*
 * Retrieve dimension names from the "t" variable, if any.
 * These dimension names will later be attached to the output variable.
 */
  dim_info = get_wrf_dim_info(2,3,ndims_t,dsizes_t);

/*
 * Calculate size of leftmost dimensions.
 */
  size_leftmost = 1;
  for(i = 0; i < ndims_qv-1; i++) size_leftmost *= dsizes_qv[i];
  nx = dsizes_qv[ndims_qv-1];
  size_rh = size_leftmost * nx;

/*
 * Test dimension sizes.
 */
  if(nx > INT_MAX) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_rh: nx = %ld is greater than INT_MAX", nx);
    return(NhlFATAL);
  }
  inx = (int) nx;

/* 
 * Allocate space for coercing input arrays.  If the input p or t
 * are already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 *
 * The output type defaults to float, unless any of the two input arrays
 * are double.
 */
  type_rh     = NCL_float;
  type_obj_rh = nclTypefloatClass;
/*
 * Allocate space for tmp_qv no matter what, because we want to set
 * values of qv that are less than zero to zero, but we don't want
 * these values retained when the function is done.
 */
  tmp_qv = (double*)malloc(nx * sizeof(double));
  if(tmp_qv == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_rh: Unable to allocate memory for coercing 'qv' to double");
    return(NhlFATAL);
  }

  if(type_qv == NCL_double) {
    type_rh     = NCL_double;
    type_obj_rh = nclTypedoubleClass;
  }

  if(type_p != NCL_double) {
    tmp_p = (double *)calloc(nx,sizeof(double));
    if(tmp_p == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_rh: Unable to allocate memory for coercing 'p' to double");
      return(NhlFATAL);
    }
  }
  else {
    type_rh     = NCL_double;
    type_obj_rh = nclTypedoubleClass;
  }

  if(type_t != NCL_double) {
    tmp_t = (double *)calloc(nx,sizeof(double));
    if(tmp_t == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_rh: Unable to allocate memory for coercing 't' to double");
      return(NhlFATAL);
    }
  }
  else {
    type_rh     = NCL_double;
    type_obj_rh = nclTypedoubleClass;
  }

/*
 * Allocate space for output array.
 */ 
  if(type_rh == NCL_double) {
    rh = (double *)calloc(size_rh,sizeof(double));
    if(rh == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_rh: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }
  else {
    rh     = (float *)calloc(size_rh,sizeof(float));
    tmp_rh = (double *)calloc(nx,sizeof(double));
    if(tmp_rh == NULL || rh == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_rh: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }
/*
 * Loop across leftmost dimensions and call the Fortran routine for each
 * one-dimensional subsection.
 */
  index_qv = 0;
  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of qv (tmp_qv) to double if necessary. Otherwise,
 * just do a memcpy. Afterwards, set all values < 0 to 0.
 */
    if(type_qv != NCL_double) {
      coerce_subset_input_double(qv,tmp_qv,index_qv,type_qv,nx,0,NULL,NULL);
    }
    else {
      (void *)memcpy((void*)((char*)tmp_qv),
                     (void*)((char*)qv + (index_qv*sizeof(double))),
                     sizeof(double)*nx);
    }
    var_zero(tmp_qv, nx);    /* Set all values < 0 to 0. */
/*
 * Coerce subsection of p (tmp_p) to double if necessary.
 */
    if(type_p != NCL_double) {
      coerce_subset_input_double(p,tmp_p,index_qv,type_p,nx,0,NULL,NULL);
    }
    else {
      tmp_p = &((double*)p)[index_qv];
    }
/*
 * Coerce subsection of t (tmp_t) to double if ncessary.
 */
    if(type_t != NCL_double) {
      coerce_subset_input_double(t,tmp_t,index_qv,type_t,nx,
                                 0,NULL,NULL);
    }
    else {
      tmp_t = &((double*)t)[index_qv];
    }

/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_rh == NCL_double) tmp_rh = &((double*)rh)[index_qv];
/*
 * Call Fortran routine.
 */
    NGCALLF(dcomputerh,DCOMPUTERH)(tmp_qv,tmp_p,tmp_t,tmp_rh,&inx);

/*
 * Coerce output back to float if necessary.
 */
    if(type_rh == NCL_float) {
      coerce_output_float_only(rh,tmp_rh,nx,index_qv);
    }

    index_qv += nx;    /* Increment index */
  }
/*
 * Free up memory.
 */
  NclFree(tmp_qv);
  if(type_p  != NCL_double) NclFree(tmp_p);
  if(type_t  != NCL_double) NclFree(tmp_t);
  if(type_rh != NCL_double) NclFree(tmp_rh);

/*
 * Set up some attributes ("description" and "units") to return.
 */
  cdescription = (char *)calloc(18,sizeof(char));
  cunits       = (char *)calloc(2,sizeof(char));
  strcpy(cdescription,"Relative Humidity");
  strcpy(cunits,"%");
  description = (NclQuark*)NclMalloc(sizeof(NclQuark));
  units       = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *description = NrmStringToQuark(cdescription);
  *units       = NrmStringToQuark(cunits);
  free(cdescription);
  free(cunits);

/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)rh,
                            NULL,
                            ndims_qv,
                            dsizes_qv,
                            TEMPORARY,
                            NULL,
                            type_obj_rh
                            );
/*
 * Set up attributes to return.
 */
  att_id = _NclAttCreate(NULL,NULL,Ncl_Att,0,NULL);

  dsizes[0] = 1;
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)description,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "description",
             att_md,
             NULL
             );
    
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)units,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "units",
             att_md,
             NULL
             );
    
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          att_id,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );

  NclFree(dim_info);
/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);
}

NhlErrorTypes wrf_slp_W( void )
{
/*
 * Input array variables
 */
  void *z, *t, *p, *q;
  double *tmp_z = NULL; 
  double *tmp_t = NULL; 
  double *tmp_p = NULL; 
  double *tmp_q = NULL; 
  int ndims_z, ndims_t, ndims_p, ndims_q;
  ng_size_t dsizes_z[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_t[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_p[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_q[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_z, type_t, type_p, type_q;
/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info = NULL;
  NclDimRec *dim_info_t = NULL;

/*
 * Output variable.
 */
  void *slp;
  NclQuark *description, *units;
  char *cdescription, *cunits;
  double *tmp_slp = NULL;
  int ndims_slp;
  ng_size_t *dsizes_slp;
  ng_size_t size_slp;
  NclBasicDataTypes type_slp;
  NclObjClass type_obj_slp;
  int errstat;
  char* errmsg;
/*
 * Various
 */
  ng_size_t i, nx, ny, nz, nxy, nxyz, size_leftmost, index_nxy, index_nxyz;
  double *tmp_t_sea_level, *tmp_t_surf, *tmp_level;
  int inx, iny, inz;
/*
 * Variables for returning the output array with attributes and/or
 * dimension names attached.
 */
  int att_id;
  ng_size_t dsizes[1];
  NclMultiDValData att_md, return_md;
  NclVar tmp_var;
  NclStackEntry return_data;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
  z = (void*)NclGetArgValue(
           0,
           4,
           &ndims_z,
           dsizes_z,
           NULL,
           NULL,
           &type_z,
           DONT_CARE);

  t = (void*)NclGetArgValue(
           1,
           4,
           &ndims_t,
           dsizes_t,
           NULL,
           NULL,
           &type_t,
           DONT_CARE);

  p = (void*)NclGetArgValue(
           2,
           4,
           &ndims_p,
           dsizes_p,
           NULL,
           NULL,
           &type_p,
           DONT_CARE);

  q = (void*)NclGetArgValue(
           3,
           4,
           &ndims_q,
           dsizes_q,
           NULL,
           NULL,
           &type_q,
           DONT_CARE);

/*
 * Error checking. Input variables must be same size, and must have at least
 * 3 dimensions.
 */
  if(ndims_z != ndims_t || ndims_z != ndims_p || ndims_z != ndims_q) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_slp: The z, t, p, and q arrays must have the same number of dimensions");
    return(NhlFATAL);
  }
  if(ndims_z < 3) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_slp: The z, t, p, and q arrays must have at least 3 dimensions");
    return(NhlFATAL);
  }
  for(i = 0; i < ndims_z; i++) {
    if(dsizes_z[i] != dsizes_t[i] || dsizes_z[i] != dsizes_p[i] ||
       dsizes_z[i] != dsizes_q[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_slp: z, t, p, and q must be the same dimensionality");
      return(NhlFATAL);
    }
  }
/*
 * Allocate space to set dimension sizes.
 */
  ndims_slp  = ndims_z-1;
  dsizes_slp = (ng_size_t*)calloc(ndims_slp,sizeof(ng_size_t));  
  if( dsizes_slp == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_slp: Unable to allocate memory for holding dimension sizes");
    return(NhlFATAL);
  }
/*
 * Set sizes for output array and calculate size of leftmost dimensions.
 * The output array will have one less dimension than the four input arrays.
 */
  size_leftmost = 1;
  for(i = 0; i < ndims_z-3; i++) {
    dsizes_slp[i] = dsizes_z[i];
    size_leftmost *= dsizes_z[i];
  }
  nx = dsizes_z[ndims_z-1];
  ny = dsizes_z[ndims_z-2];
  nz = dsizes_z[ndims_z-3];
  dsizes_slp[ndims_slp-1] = nx;
  dsizes_slp[ndims_slp-2] = ny;
  nxy  = nx * ny;
  nxyz = nxy * nz;
  size_slp = size_leftmost * nxy;
/*
 * Test dimension sizes.
 */
  if((nx > INT_MAX) || (ny > INT_MAX) || (nz > INT_MAX)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_slp: nx, ny, and/or nz is greater than INT_MAX");
    return(NhlFATAL);
  }
  inx = (int) nx;
  iny = (int) ny;
  inz = (int) nz;

/*
 * Get dimension info to see if we have named dimensions.
 * This will be used for return variable.
 */
  dim_info_t = get_wrf_dim_info(1,4,ndims_t,dsizes_t);
  if(dim_info_t != NULL) {
    dim_info = malloc(sizeof(NclDimRec)*ndims_slp);
    if(dim_info == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_slp: Unable to allocate memory for holding dimension information");
      return(NhlFATAL);
    }
    for(i = 0; i < ndims_z-3; i++) {
      dim_info[i] = dim_info_t[i];
    }
    dim_info[ndims_slp-1] = dim_info_t[ndims_t-1];
    dim_info[ndims_slp-2] = dim_info_t[ndims_t-2];
  }

/* 
 * Allocate space for coercing input arrays.  If the input q, p, or t
 * are already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 *
 * The output type defaults to float, unless any of the two input arrays
 * are double.
 */
  type_slp     = NCL_float;
  type_obj_slp = nclTypefloatClass;

/*
 * Allocate space for tmp_q no matter what, because we want to set
 * values of q that are less than zero to zero, but we don't want
 * these values retained when the function is done.
 */
  tmp_q = (double*)malloc(nxyz * sizeof(double));
  if(tmp_q == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_slp: Unable to allocate memory for coercing 'q' to double");
    return(NhlFATAL);
  }

  if(type_q == NCL_double) {
    type_slp     = NCL_double;
    type_obj_slp = nclTypedoubleClass;
  }

  if(type_z != NCL_double) {
    tmp_z = (double *)calloc(nxyz,sizeof(double));
    if(tmp_z == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_slp: Unable to allocate memory for coercing 'z' to double");
      return(NhlFATAL);
    }
  }
  else {
    type_slp     = NCL_double;
    type_obj_slp = nclTypedoubleClass;
  }

  if(type_t != NCL_double) {
    tmp_t = (double *)calloc(nxyz,sizeof(double));
    if(tmp_t == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_slp: Unable to allocate memory for coercing 't' to double");
      return(NhlFATAL);
    }
  }
  else {
    type_slp     = NCL_double;
    type_obj_slp = nclTypedoubleClass;
  }

  if(type_p != NCL_double) {
    tmp_p = (double *)calloc(nxyz,sizeof(double));
    if(tmp_p == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_slp: Unable to allocate memory for coercing 'p' to double");
      return(NhlFATAL);
    }
  }
  else {
    type_slp     = NCL_double;
    type_obj_slp = nclTypedoubleClass;
  }

/*
 * Allocate space for work arrays.
 */ 
  tmp_t_sea_level = (double *)calloc(nxy,sizeof(double));
  tmp_t_surf      = (double *)calloc(nxy,sizeof(double));
  tmp_level       = (double *)calloc(nxy,sizeof(double));
  if(tmp_t_sea_level == NULL || tmp_t_surf == NULL || tmp_level == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_slp: Unable to allocate memory for temporary arrays");
      return(NhlFATAL);
  }

/*
 * Allocate space for output array.
 */ 
  if(type_slp == NCL_double) {
    slp = (double *)calloc(size_slp,sizeof(double));
    if(slp == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_slp: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }
  else {
    slp     = (float *)calloc(size_slp,sizeof(float));
    tmp_slp = (double *)calloc(nxy,sizeof(double));
    if(tmp_slp == NULL || slp == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_slp: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }

  /* Allocate space for errmsg*/
  errmsg = (char *) calloc(ERRLEN, sizeof(char))

/*
 * Loop across leftmost dimensions and call the Fortran routine
 * for each three-dimensional subsection.
 */
  index_nxy = index_nxyz = 0;
  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of q (tmp_q) to double if necessary.  Otherwise,
 * just do a memcpy. Afterwards, set all values < 0 to 0.
 */
    if(type_q != NCL_double) {
      coerce_subset_input_double(q,tmp_q,index_nxyz,type_q,nxyz,0,NULL,NULL);
    }
    else {
      (void *)memcpy((void*)((char*)tmp_q),
                     (void*)((char*)q + (index_nxyz*sizeof(double))),
                     sizeof(double)*nxyz);
    }
    var_zero(tmp_q, nxyz);   /* Set all values < 0 to 0. */
/*
 * Coerce subsection of z (tmp_z) to double if necessary.
 */
    if(type_z != NCL_double) {
      coerce_subset_input_double(z,tmp_z,index_nxyz,type_z,nxyz,0,NULL,NULL);
    }
    else {
      tmp_z = &((double*)z)[index_nxyz];
    }
/*
 * Coerce subsection of p (tmp_p) to double if necessary.
 */
    if(type_p != NCL_double) {
      coerce_subset_input_double(p,tmp_p,index_nxyz,type_p,nxyz,0,NULL,NULL);
    }
    else {
      tmp_p = &((double*)p)[index_nxyz];
    }
/*
 * Coerce subsection of t (tmp_t) to double if ncessary.
 */
    if(type_t != NCL_double) {
      coerce_subset_input_double(t,tmp_t,index_nxyz,type_t,nxyz,0,NULL,NULL);
    }
    else {
      tmp_t = &((double*)t)[index_nxyz];
    }

/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_slp == NCL_double) tmp_slp = &((double*)slp)[index_nxy];
/*
 * Call Fortran routine.
 */
    errstat = 0;
    errmsg = "";
    NGCALLF(dcomputeseaprs,DCOMPUTESEAPRS)(&inx,&iny,&inz,tmp_z,tmp_t,tmp_p,
                                           tmp_q,tmp_slp,tmp_t_sea_level,
                                           tmp_t_surf,tmp_level,&errstat,
										   errmsg, ERRLEN);

    /* Terminate if there was an error */
    if (errstat != 0) {
    	fprintf(stderr, errmsg);
    	exit(errstat);
    }
/*
 * Coerce output back to float if necessary.
 */
    if(type_slp == NCL_float) {
      coerce_output_float_only(slp,tmp_slp,nxy,index_nxy);
    }

    index_nxyz += nxyz;    /* Increment indices */
    index_nxy  += nxy;
  }
/*
 * Free up memory.
 */
  NclFree(tmp_q);
  if(type_p   != NCL_double) NclFree(tmp_p);
  if(type_t   != NCL_double) NclFree(tmp_t);
  if(type_z   != NCL_double) NclFree(tmp_z);
  if(type_slp != NCL_double) NclFree(tmp_slp);

  NclFree(tmp_t_sea_level);
  NclFree(tmp_t_surf);
  NclFree(tmp_level);

  NclFree(errmsg);

/*
 * Set up some attributes ("description" and "units") to return.
 */
  cdescription = (char *)calloc(19,sizeof(char));
  cunits       = (char *)calloc(4,sizeof(char));
  strcpy(cdescription,"Sea Level Pressure");
  strcpy(cunits,"hPa");
  description = (NclQuark*)NclMalloc(sizeof(NclQuark));
  units       = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *description = NrmStringToQuark(cdescription);
  *units       = NrmStringToQuark(cunits);
  free(cdescription);
  free(cunits);

/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)slp,
                            NULL,
                            ndims_slp,
                            dsizes_slp,
                            TEMPORARY,
                            NULL,
                            type_obj_slp
                            );
/*
 * Set up attributes to return.
 */
  att_id = _NclAttCreate(NULL,NULL,Ncl_Att,0,NULL);

  dsizes[0] = 1;
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)description,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "description",
             att_md,
             NULL
             );
    
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)units,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "units",
             att_md,
             NULL
             );
    
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          att_id,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );

  NclFree(dsizes_slp);
  if(dim_info != NULL) NclFree(dim_info);
  NclFree(dim_info_t);

/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);

}

NhlErrorTypes wrf_interp_3d_z_W( void )
{
/*
 * Input array variables
 */
  void *v3d, *z, *loc;
  double *tmp_v3d = NULL; 
  double *tmp_z = NULL; 
  double *tmp_loc = NULL; 
  int ndims_v3d, ndims_z;
  ng_size_t dsizes_v3d[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_z[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_v3d, type_z, type_loc;
/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info_v3d;
  NclDimRec *dim_info = NULL;

/*
 * Variables for retrieving attributes from "v3d".
 */
  NclAttList  *attr_list;
  NclAtt  attr_obj;
  NclStackEntry   stack_entry;
  NrmQuark *description, *units;
  char *cdesc = NULL;
  char *cunits = NULL;
  logical found_desc = False, found_units = False;
/*
 * Output variable.
 */
  void *v2d;
  double *tmp_v2d = NULL;
  int ndims_v2d;
  ng_size_t *dsizes_v2d;
  ng_size_t size_v2d;
  NclBasicDataTypes type_v2d;
  NclObjClass type_obj_v2d;
  NclScalar missing_v2d;

/*
 * Variables for returning the output array with attributes attached.
 */
  int att_id;
  ng_size_t dsizes[1];
  NclMultiDValData att_md, return_md;
  NclVar tmp_var;
  NclStackEntry return_data;
  NclQuark *qdesc, *qunits;

/*
 * Various
 */
  ng_size_t i, nx, ny, nz, nxy, nxyz, size_leftmost, index_v3d, index_v2d;
  int inx, iny, inz;
  double vmsg;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
  v3d = (void*)NclGetArgValue(
           0,
           3,
           &ndims_v3d,
           dsizes_v3d,
           NULL,
           NULL,
           &type_v3d,
           DONT_CARE);

  z = (void*)NclGetArgValue(
           1,
           3,
           &ndims_z,
           dsizes_z,
           NULL,
           NULL,
           &type_z,
           DONT_CARE);

  loc = (void*)NclGetArgValue(
           2,
           3,
           NULL,
           NULL,
           NULL,
           NULL,
           &type_loc,
           DONT_CARE);

/*
 * Error checking. First two input variables must be same size.
 */
  if(ndims_v3d < 3) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_3d_z: The v3d and z arrays must have at least 3 dimensions");
    return(NhlFATAL);
  }

  if(ndims_v3d != ndims_z) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_3d_z: The v3d and z arrays must have the same number of dimensions");
    return(NhlFATAL);
  }

  for(i = 0; i < ndims_v3d; i++) {
    if(dsizes_v3d[i] != dsizes_z[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_3d_z: v3d and z must be the same dimensionality");
      return(NhlFATAL);
    }
  }

/*
 * Check if v3d has any attributes, namely "description" or "units".
 * These attributes will be attached to the return variable v2d.
 */
  stack_entry = _NclGetArg(0, 3, DONT_CARE);
  switch (stack_entry.kind) {
  case NclStk_VAR:
    if (stack_entry.u.data_var->var.att_id != -1) {
      attr_obj = (NclAtt) _NclGetObj(stack_entry.u.data_var->var.att_id);
      if (attr_obj == NULL) {
        break;
      }
    }
    else {
/*
 * att_id == -1 ==> no optional args given.
 */
      break;
    }
/* 
 * Get optional arguments. If none are specified, then return
 * missing values.
 */
    if (attr_obj->att.n_atts == 0) {
      break;
    }
    else {
/*
 * Get list of attributes.
 */
      attr_list = attr_obj->att.att_list;
/*
 * Loop through attributes and check them.
 */
      while (attr_list != NULL) {
        if ((strcmp(attr_list->attname, "description")) == 0) {
          description = (NrmQuark *) attr_list->attvalue->multidval.val;
          cdesc       = NrmQuarkToString(*description);
          found_desc  = True;
        }
        if ((strcmp(attr_list->attname, "units")) == 0) {
          units  = (NrmQuark *) attr_list->attvalue->multidval.val;
          cunits = NrmQuarkToString(*units);
          found_units  = True;
        }
        attr_list = attr_list->next;
      }
    }
  default:
    break;
  }

/*
 * Calculate size of leftmost dimensions and set dimension sizes for 
 * output array.
 *
 * The output array will have one less dimension than v3d/z input arrays.
 */
  nx = dsizes_v3d[ndims_v3d-1];
  ny = dsizes_v3d[ndims_v3d-2];
  nz = dsizes_v3d[ndims_v3d-3];
  nxy  = nx * ny;
  nxyz = nxy * nz;

/*
 * Test dimension sizes.
 */
  if((nx > INT_MAX) || (ny > INT_MAX) || (nz > INT_MAX)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_3d_z: nx, ny, and/or nz is greater than INT_MAX");
    return(NhlFATAL);
  }
  inx = (int) nx;
  iny = (int) ny;
  inz = (int) nz;

  ndims_v2d = ndims_v3d-1;
  dsizes_v2d = (ng_size_t*)calloc(ndims_v2d,sizeof(ng_size_t));  
  if( dsizes_v2d == NULL ) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_3d_z: Unable to allocate memory for holding dimension sizes");
    return(NhlFATAL);
  }
  size_leftmost = 1;
  for(i = 0; i < ndims_v3d-3; i++) {
    dsizes_v2d[i] = dsizes_v3d[i];
    size_leftmost *= dsizes_v3d[i];
  }
  dsizes_v2d[ndims_v2d-2] = ny;
  dsizes_v2d[ndims_v2d-1] = nx;

  size_v2d = size_leftmost * nxy;

/*
 * Retrieve dimension names from the "v3d" variable, if any.
 * These dimension names will later be attached to the output variable.
 */
  dim_info_v3d = get_wrf_dim_info(0,3,ndims_v3d,dsizes_v3d);
  if(dim_info_v3d != NULL) {
    dim_info = malloc(sizeof(NclDimRec)*ndims_v2d);
    if(dim_info == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_3d_z: Unable to allocate memory for holding dimension information");
      return(NhlFATAL);
    }
    for(i = 0; i < ndims_v3d-3; i++) {
      dim_info[i] = dim_info_v3d[i];
    }
    dim_info[ndims_v2d-1] = dim_info_v3d[ndims_v3d-1];
    dim_info[ndims_v2d-2] = dim_info_v3d[ndims_v3d-2];
  }

/* 
 * Allocate space for coercing input arrays.  If the input v3d or z
 * are already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 *
 * The output type defaults to float, unless any of the two input arrays
 * are double.
 */
  type_v2d     = NCL_float;
  type_obj_v2d = nclTypefloatClass;
  if(type_v3d != NCL_double) {
    tmp_v3d = (double *)calloc(nxyz,sizeof(double));
    if(tmp_v3d == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_3d_z: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_v2d     = NCL_double;
    type_obj_v2d = nclTypedoubleClass;
  }
  if(type_z != NCL_double) {
    tmp_z = (double *)calloc(nxyz,sizeof(double));
    if(tmp_z == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_3d_z: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_v2d     = NCL_double;
    type_obj_v2d = nclTypedoubleClass;
  }
/*
 * Coerce loc (tmp_loc) to double if ncessary.
 */
  tmp_loc = coerce_input_double(loc,type_loc,1,0,NULL,NULL);

/*
 * Allocate space for output array.
 */ 
  if(type_v2d == NCL_double) {
    v2d = (double *)calloc(size_v2d,sizeof(double));
    if(v2d == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_3d_z: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
    missing_v2d.doubleval = ((NclTypeClass)nclTypedoubleClass)->type_class.default_mis.doubleval;
    vmsg = missing_v2d.doubleval;
  }
  else {
    v2d     = (float *)calloc(size_v2d,sizeof(float));
    tmp_v2d = (double *)calloc(nxy,sizeof(double));
    if(tmp_v2d == NULL || v2d == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_3d_z: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
    missing_v2d.floatval = ((NclTypeClass)nclTypefloatClass)->type_class.default_mis.floatval;
    vmsg = (double)missing_v2d.floatval;

  }
/*
 * Loop across leftmost dimensions and call the Fortran routine
 * for each three-dimensional subsection.
 */
  index_v2d = index_v3d = 0;
  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of v3d (tmp_v3d) to double if necessary.
 */
    if(type_v3d != NCL_double) {
      coerce_subset_input_double(v3d,tmp_v3d,index_v3d,type_v3d,nxyz,
                                 0,NULL,NULL);
    }
    else {
      tmp_v3d = &((double*)v3d)[index_v3d];
    }
/*
 * Coerce subsection of z (tmp_z) to double if necessary.
 */
    if(type_z != NCL_double) {
      coerce_subset_input_double(z,tmp_z,index_v3d,type_z,nxyz,0,NULL,NULL);
    }
    else {
      tmp_z = &((double*)z)[index_v3d];
    }

/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_v2d == NCL_double) tmp_v2d = &((double*)v2d)[index_v2d];
/*
 * Call Fortran routine.
 */
    NGCALLF(dinterp3dz,DINTERP3DZ)(tmp_v3d,tmp_v2d,tmp_z,tmp_loc,
                                   &inx,&iny,&inz,&vmsg);
/*
 * Coerce output back to float if necessary.
 */
    if(type_v2d == NCL_float) {
      coerce_output_float_only(v2d,tmp_v2d,nxy,index_v2d);
    }

    index_v3d += nxyz;
    index_v2d += nxy;
  }
/*
 * Free up memory.
 */
  if(type_v3d != NCL_double) NclFree(tmp_v3d);
  if(type_z   != NCL_double) NclFree(tmp_z);
  if(type_loc != NCL_double) NclFree(tmp_loc);
  if(type_v2d != NCL_double) NclFree(tmp_v2d);

/*
 * If v3d had a "description" or units attribute, return them with
 * the output variable as an attribute.  Otherwise, return a
 * blank string for description, and nothing for units.
 */
  if(!found_desc) {
    cdesc = (char *)calloc(2,sizeof(char));
    strcpy(cdesc," ");
  }
/*
 * I don't think we can return "description" or "units" here, because 
 * they are attached to an NCL input parameter. It could screw things up
 * if we try to return it as an attribute with the output variable.
 * Instead, create a new description and units "quark" variable.
 */
  qdesc  = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *qdesc = NrmStringToQuark(cdesc);
  if (!found_desc)
          free(cdesc);

  if(found_units) {
    qunits  = (NclQuark*)NclMalloc(sizeof(NclQuark));
    *qunits = NrmStringToQuark(cunits);
  }

/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)v2d,
                            &missing_v2d,
                            ndims_v2d,
                            dsizes_v2d,
                            TEMPORARY,
                            NULL,
                            type_obj_v2d
                            );
/*
 * Set up attributes to return.
 */
  att_id = _NclAttCreate(NULL,NULL,Ncl_Att,0,NULL);

  dsizes[0] = 1;
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)qdesc,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "description",
             att_md,
             NULL
             );
    
  if(found_units) {
    att_md = _NclCreateVal(
                           NULL,
                           NULL,
                           Ncl_MultiDValData,
                           0,
                           (void*)qunits,
                           NULL,
                           1,
                           dsizes,
                           TEMPORARY,
                           NULL,
                           (NclObjClass)nclTypestringClass
                           );
    _NclAddAtt(
               att_id,
             "units",
               att_md,
               NULL
               );
  }
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          att_id,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );
  NclFree(dsizes_v2d);
  if(dim_info != NULL) NclFree(dim_info);
  NclFree(dim_info_v3d);

/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);
}


NhlErrorTypes wrf_interp_2d_xy_W( void )
{
/*
 * Input array variables
 */
  void *v3d, *xy;
  double *tmp_v3d = NULL; 
  double *tmp_xy = NULL; 
  int ndims_v3d, ndims_xy;
  ng_size_t dsizes_v3d[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_xy[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_v3d, type_xy;

/*
 * Variables for retrieving attributes from "v3d".
 */
  NclAttList  *attr_list;
  NclAtt  attr_obj;
  NclStackEntry   stack_entry;
  NrmQuark *description, *units;
  char *cdesc = NULL;
  char *cunits = NULL;
  logical found_desc = False, found_units = False;
/*
 * Output variable.
 */
  void *v2d;
  double *tmp_v2d = NULL;
  int ndims_v2d;
  ng_size_t *dsizes_v2d, size_v2d;
  NclBasicDataTypes type_v2d;
  NclObjClass type_obj_v2d;
  NclScalar missing_v2d;

/*
 * Variables for returning the output array with attributes attached.
 */
  int att_id;
  ng_size_t dsizes[1];
  NclMultiDValData att_md, return_md;
  NclVar tmp_var;
  NclStackEntry return_data;
  NclQuark *qdesc, *qunits;

/*
 * Various
 */
  ng_size_t i, nx, ny, nz, nxnynz, nxy, nxy_nz , nxy_2, size_leftmost;
  ng_size_t index_v3d, index_v2d, index_xy;
  int inx, iny, inz, inxy;
/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
  v3d = (void*)NclGetArgValue(
           0,
           2,
           &ndims_v3d,
           dsizes_v3d,
           NULL,
           NULL,
           &type_v3d,
           DONT_CARE);

  xy = (void*)NclGetArgValue(
           1,
           2,
           &ndims_xy,
           dsizes_xy,
           NULL,
           NULL,
           &type_xy,
           DONT_CARE);

/*
 * Error checking.
 */
  if(ndims_v3d < 3) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_2d_xy: The v3d array must be at least 3-dimensional");
    return(NhlFATAL);
  }
  if(ndims_v3d != (ndims_xy+1)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_2d_xy: The v3d array must have one more dimension than the xy array");
    return(NhlFATAL);
  }
  if(dsizes_xy[ndims_xy-1] != 2) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_2d_xy: The rightmost dimension of xy must be 2");
    return(NhlFATAL);
  }
  nz  = dsizes_v3d[ndims_v3d-3];
  ny  = dsizes_v3d[ndims_v3d-2];
  nx  = dsizes_v3d[ndims_v3d-1];
  nxy = dsizes_xy[ndims_xy-2];
  nxnynz   = nx * ny * nz;
  nxy_nz   = nxy * nz;
  nxy_2    = nxy * 2;

/*
 * Test dimension sizes.
 */
  if((nxy > INT_MAX) || (nx > INT_MAX) || (ny > INT_MAX) || (nz > INT_MAX)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_2d_xy: one or more dimension sizes is greater than INT_MAX", nxy);
    return(NhlFATAL);
  }
  inx = (int) nx;
  iny = (int) ny;
  inz = (int) nz;
  inxy = (int) nxy;

/*
 * Check leftmost dimensions, if any, and calculate their size.
 * Also set dimension sizes for output array.
 */
  ndims_v2d = ndims_xy;     /* leftmost dims x nz x nxy */
  dsizes_v2d = (ng_size_t*)calloc(ndims_v2d,sizeof(ng_size_t));  
  if( dsizes_v2d == NULL ) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_2d_xy: Unable to allocate memory for holding dimension sizes");
    return(NhlFATAL);
  }
  size_leftmost = 1;
  for(i = 0; i < ndims_v3d-3; i++) {
    if(dsizes_v3d[i] != dsizes_xy[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_2d_xy: The leftmost dimensions of v3d and xy must be the same");
      return(NhlFATAL);
    }
    dsizes_v2d[i] = dsizes_v3d[i];
    size_leftmost *= dsizes_v3d[i];
  }
  dsizes_v2d[ndims_v2d-2] = nz;
  dsizes_v2d[ndims_v2d-1] = nxy;

  size_v2d = size_leftmost * nxy_nz;

/*
 * Check if v3d has any attributes, namely "description" or "units".
 * These attributes will be attached to the return variable v2d.
 */
  stack_entry = _NclGetArg(0, 2, DONT_CARE);
  switch (stack_entry.kind) {
  case NclStk_VAR:
    if (stack_entry.u.data_var->var.att_id != -1) {
      attr_obj = (NclAtt) _NclGetObj(stack_entry.u.data_var->var.att_id);
      if (attr_obj == NULL) {
        break;
      }
    }
    else {
/*
 * att_id == -1 ==> no optional args given.
 */
      break;
    }
/* 
 * Get optional arguments. If none are specified, then return
 * missing values.
 */
    if (attr_obj->att.n_atts == 0) {
      break;
    }
    else {
/*
 * Get list of attributes.
 */
      attr_list = attr_obj->att.att_list;
/*
 * Loop through attributes and check them.
 */
      while (attr_list != NULL) {
        if ((strcmp(attr_list->attname, "description")) == 0) {
          description = (NrmQuark *) attr_list->attvalue->multidval.val;
          cdesc       = NrmQuarkToString(*description);
          found_desc  = True;
        }
        if ((strcmp(attr_list->attname, "units")) == 0) {
          units  = (NrmQuark *) attr_list->attvalue->multidval.val;
          cunits = NrmQuarkToString(*units);
          found_units  = True;
        }
        attr_list = attr_list->next;
      }
    }
  default:
    break;
  }

/* 
 * Allocate space for coercing input arrays.  If the input v3d or xy
 * are already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 *
 * The output type defaults to float, unless any of the two input arrays
 * are double.
 */
  type_v2d     = NCL_float;
  type_obj_v2d = nclTypefloatClass;
  if(type_v3d != NCL_double) {
    tmp_v3d = (double *)calloc(nxnynz,sizeof(double));
    if(tmp_v3d == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_2d_xy: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_v2d     = NCL_double;
    type_obj_v2d = nclTypedoubleClass;
  }
  if(type_xy != NCL_double) {
    tmp_xy = (double *)calloc(nxy_2,sizeof(double));
    if(tmp_xy == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_2d_xy: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_v2d     = NCL_double;
    type_obj_v2d = nclTypedoubleClass;
  }

/*
 * Allocate space for output array.
 */ 
  if(type_v2d == NCL_double) {
    v2d = (double *)calloc(size_v2d,sizeof(double));
    if(v2d == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_2d_xy: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
    missing_v2d.doubleval = ((NclTypeClass)nclTypedoubleClass)->type_class.default_mis.doubleval;
  }
  else {
    v2d     = (float *)calloc(size_v2d,sizeof(float));
    tmp_v2d = (double *)calloc(nxy_nz,sizeof(double));
    if(tmp_v2d == NULL || v2d == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_2d_xy: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
    missing_v2d.floatval = ((NclTypeClass)nclTypefloatClass)->type_class.default_mis.floatval;
  }

/*
 * Loop across leftmost dimensions and call the Fortran routine
 * for reach three-dimensional subsection.
 */
  index_v3d = index_v2d = index_xy = 0;
  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of v3d (tmp_v3d) to double if necessary.
 */
    if(type_v3d != NCL_double) {
      coerce_subset_input_double(v3d,tmp_v3d,index_v3d,type_v3d,nxnynz,
                                 0,NULL,NULL);
    }
    else {
      tmp_v3d = &((double*)v3d)[index_v3d];
    }
/*
 * Coerce subsection of xy (tmp_xy) to double if necessary.
 */
    if(type_xy != NCL_double) {
      coerce_subset_input_double(xy,tmp_xy,index_xy,type_xy,nxy_2,0,NULL,NULL);
    }
    else {
      tmp_xy = &((double*)xy)[index_xy];
    }

/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_v2d == NCL_double) tmp_v2d = &((double*)v2d)[index_v2d];
/*
 * Call Fortran routine.
 */
    NGCALLF(dinterp2dxy,DINTERP2DXY)(tmp_v3d,tmp_v2d,tmp_xy,&inx,&iny,&inz,&inxy);
/*
 * Coerce output back to float if necessary.
 */
    if(type_v2d == NCL_float) {
      coerce_output_float_only(v2d,tmp_v2d,nxy_nz,index_v2d);
    }

    index_v3d += nxnynz;    /* Increment indices */
    index_v2d += nxy_nz;
    index_xy  += nxy_2;
  }
/*
 * Free up memory.
 */
  if(type_v3d != NCL_double) NclFree(tmp_v3d);
  if(type_xy  != NCL_double) NclFree(tmp_xy);
  if(type_v2d != NCL_double) NclFree(tmp_v2d);

/*
 * If v3d had a "description" or units attribute, return them with
 * the output variable as an attribute.  Otherwise, return a
 * blank string for description, and nothing for units.
 */
  if(!found_desc) {
    cdesc = (char *)calloc(2,sizeof(char));
    strcpy(cdesc," ");
  }
/*
 * I don't think we can return "description" or "units" here, because 
 * they are attached to an NCL input parameter. It could screw things up
 * if we try to return it as an attribute with the output variable.
 * Instead, create a new description and units "quark" variable.
 */
  qdesc  = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *qdesc = NrmStringToQuark(cdesc);
  if (!found_desc)
          free(cdesc);
  if(found_units) {
    qunits  = (NclQuark*)NclMalloc(sizeof(NclQuark));
    *qunits = NrmStringToQuark(cunits);
  }
/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)v2d,
                            &missing_v2d,
                            ndims_v2d,
                            dsizes_v2d,
                            TEMPORARY,
                            NULL,
                            type_obj_v2d
                            );
/*
 * Set up attributes to return.
 */
  att_id = _NclAttCreate(NULL,NULL,Ncl_Att,0,NULL);

  dsizes[0] = 1;
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)qdesc,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "description",
             att_md,
             NULL
             );
    
  if(found_units) {
    att_md = _NclCreateVal(
                           NULL,
                           NULL,
                           Ncl_MultiDValData,
                           0,
                           (void*)qunits,
                           NULL,
                           1,
                           dsizes,
                           TEMPORARY,
                           NULL,
                           (NclObjClass)nclTypestringClass
                           );
    _NclAddAtt(
               att_id,
             "units",
               att_md,
               NULL
               );
  }
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          NULL,
                          att_id,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );
  NclFree(dsizes_v2d);
/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);
}


NhlErrorTypes wrf_interp_1d_W( void )
{
/*
 * Input array variables
 */
  void *v_in, *z_in, *z_out;
  double *tmp_v_in = NULL; 
  double *tmp_z_in = NULL; 
  double *tmp_z_out = NULL; 
  int ndims_v_in, ndims_z_in, ndims_z_out;
  ng_size_t dsizes_v_in[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_z_in[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_z_out[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_v_in, type_z_in, type_z_out;
/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info;

/*
 * Variables for retrieving attributes from "v3d".
 */
  NclAttList  *attr_list;
  NclAtt  attr_obj;
  NclStackEntry   stack_entry;
  NrmQuark *description, *units;
  char *cdesc = NULL;
  char *cunits = NULL;
  logical found_desc = False, found_units = False;

/*
 * Output variable.
 */
  void *v_out;
  double *tmp_v_out = NULL;
  double v_out_msg;
  ng_size_t *dsizes_v_out, size_v_out;
  NclBasicDataTypes type_v_out;
  NclObjClass type_obj_v_out;
  NclScalar missing_v_out;

/*
 * Variables for returning the output array with attributes attached.
 */
  int att_id;
  ng_size_t dsizes[1];
  NclMultiDValData att_md, return_md;
  NclVar tmp_var;
  NclStackEntry return_data;
  NclQuark *qdesc, *qunits;

/*
 * Various
 */
  ng_size_t i, nz_in, nz_out, size_leftmost, index_v_in, index_v_out;
  int inz_in, inz_out;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
  v_in = (void*)NclGetArgValue(
           0,
           3,
           &ndims_v_in,
           dsizes_v_in,
           NULL,
           NULL,
           &type_v_in,
           DONT_CARE);

  z_in = (void*)NclGetArgValue(
           1,
           3,
           &ndims_z_in,
           dsizes_z_in,
           NULL,
           NULL,
           &type_z_in,
           DONT_CARE);

  z_out = (void*)NclGetArgValue(
           2,
           3,
           &ndims_z_out,
           dsizes_z_out,
           NULL,
           NULL,
           &type_z_out,
           DONT_CARE);

/*
 * Error checking.
 */
  if(ndims_v_in != ndims_z_in || ndims_v_in != ndims_z_out) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_1d: The v_in, z_in, and z_out arrays must be the same number of dimensions");
    return(NhlFATAL);
  }
  nz_in  = dsizes_z_in[ndims_z_in-1];
  nz_out = dsizes_z_out[ndims_z_out-1];
  if(dsizes_v_in[ndims_v_in-1] != nz_in) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_1d: The rightmost dimension of v_in and z_in must be the same");
    return(NhlFATAL);
  }

/*
 * Test dimension sizes.
 */
  if((nz_in > INT_MAX) || (nz_out > INT_MAX)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_1d: nz_in and/or nz_out is greater than INT_MAX");
    return(NhlFATAL);
  }
  inz_in = (int) nz_in;
  inz_out = (int) nz_out;

/*
 * Retrieve dimension names from the "v3d" variable, if any.
 * These dimension names will later be attached to the output variable.
 */
  dim_info = get_wrf_dim_info(0,3,ndims_v_in,dsizes_v_in);

/*
 * Check if v_in has any attributes, namely "description" or "units".
 * These attributes will be attached to the return variable v_out.
 */
  stack_entry = _NclGetArg(0, 3, DONT_CARE);
  switch (stack_entry.kind) {
  case NclStk_VAR:
    if (stack_entry.u.data_var->var.att_id != -1) {
      attr_obj = (NclAtt) _NclGetObj(stack_entry.u.data_var->var.att_id);
      if (attr_obj == NULL) {
        break;
      }
    }
    else {
/*
 * att_id == -1 ==> no optional args given.
 */
      break;
    }
/* 
 * Get optional arguments. If none are specified, then return
 * missing values.
 */
    if (attr_obj->att.n_atts == 0) {
      break;
    }
    else {
/*
 * Get list of attributes.
 */
      attr_list = attr_obj->att.att_list;
/*
 * Loop through attributes and check them.
 */
      while (attr_list != NULL) {
        if ((strcmp(attr_list->attname, "description")) == 0) {
          description = (NrmQuark *) attr_list->attvalue->multidval.val;
          cdesc       = NrmQuarkToString(*description);
          found_desc  = True;
        }
        if ((strcmp(attr_list->attname, "units")) == 0) {
          units  = (NrmQuark *) attr_list->attvalue->multidval.val;
          cunits = NrmQuarkToString(*units);
          found_units  = True;
        }
        attr_list = attr_list->next;
      }
    }
  default:
    break;
  }

/*
 * Calculate leftmost dimensions, if any, and check their sizes.
 * Also set dimension sizes for output array.
 */
  dsizes_v_out = (ng_size_t*)calloc(ndims_z_out,sizeof(ng_size_t));  
  if( dsizes_v_out == NULL ) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_1d: Unable to allocate memory for holding dimension sizes");
    return(NhlFATAL);
  }

  size_leftmost = 1;
  for(i = 0; i < ndims_v_in-1; i++ ) {
    if(dsizes_v_in[i] != dsizes_z_in[i] || 
       dsizes_v_in[i] != dsizes_z_out[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_1d: The input arrays must be the same dimensionality");
      return(NhlFATAL);
    }
    dsizes_v_out[i] = dsizes_v_in[i];
    size_leftmost *= dsizes_v_in[i];
  }
  dsizes_v_out[ndims_v_in-1] = nz_out;
  size_v_out = size_leftmost * nz_out;

/* 
 * Allocate space for coercing input arrays.  If the input arrays
 * are already double, then we don't need to allocate space for the
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 *
 * The output type defaults to float, unless any of the two input arrays
 * are double.
 */
  type_v_out     = NCL_float;
  type_obj_v_out = nclTypefloatClass;
  if(type_v_in != NCL_double) {
    tmp_v_in = (double *)calloc(nz_in,sizeof(double));
    if(tmp_v_in == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_1d: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_v_out     = NCL_double;
    type_obj_v_out = nclTypedoubleClass;
  }
  if(type_z_in != NCL_double) {
    tmp_z_in = (double *)calloc(nz_in,sizeof(double));
    if(tmp_z_in == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_1d: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_v_out     = NCL_double;
    type_obj_v_out = nclTypedoubleClass;
  }

  if(type_z_out != NCL_double) {
    tmp_z_out = (double *)calloc(nz_out,sizeof(double));
    if(tmp_z_out == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_1d: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_v_out     = NCL_double;
    type_obj_v_out = nclTypedoubleClass;
  }

/*
 * Allocate space for output array.
 */ 
  if(type_v_out == NCL_double) {
    v_out = (double *)calloc(size_v_out,sizeof(double));
    if(v_out == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_1d: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
    v_out_msg = missing_v_out.doubleval = ((NclTypeClass)nclTypedoubleClass)->type_class.default_mis.doubleval;
  }
  else {
    v_out     = (float *)calloc(size_v_out,sizeof(float));
    tmp_v_out = (double *)calloc(nz_out,sizeof(double));
    if(tmp_v_out == NULL || v_out == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_interp_1d: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
    v_out_msg = (double)((NclTypeClass)nclTypefloatClass)->type_class.default_mis.floatval;
    missing_v_out.floatval = ((NclTypeClass)nclTypefloatClass)->type_class.default_mis.floatval;
  }

/*
 * Loop across leftmost dimensions and call the Fortran routine
 * for reach one-dimensional subsection.
 */
  index_v_out = index_v_in = 0;
  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of v_in (tmp_v_in) to double if necessary.
 */
    if(type_v_in != NCL_double) {
      coerce_subset_input_double(v_in,tmp_v_in,index_v_in,type_v_in,nz_in,
                                 0,NULL,NULL);
    }
    else {
      tmp_v_in = &((double*)v_in)[index_v_in];
    }
/*
 * Coerce subsection of z_in (tmp_z_in) to double if necessary.
 */
    if(type_z_in != NCL_double) {
      coerce_subset_input_double(z_in,tmp_z_in,index_v_in,type_z_in,nz_in,
                                 0,NULL,NULL);
    }
    else {
      tmp_z_in = &((double*)z_in)[index_v_in];
    }

/*
 * Coerce subsection of z_out (tmp_z_out) to double if necessary.
 */
    if(type_z_out != NCL_double) {
      coerce_subset_input_double(z_out,tmp_z_out,index_v_out,type_z_out,
                                 nz_out,0,NULL,NULL);
    }
    else {
      tmp_z_out = &((double*)z_out)[index_v_out];
    }

/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_v_out == NCL_double) tmp_v_out = &((double*)v_out)[index_v_out];
/*
 * Call Fortran routine.
 */
    NGCALLF(dinterp1d,DINTERP1D)(tmp_v_in,tmp_v_out,tmp_z_in,tmp_z_out,&inz_in,
                                 &inz_out,&v_out_msg);
/*
 * Coerce output back to float if necessary.
 */
    if(type_v_out == NCL_float) {
      coerce_output_float_only(v_out,tmp_v_out,nz_out,index_v_out);
    }

    index_v_in  += nz_in;
    index_v_out += nz_out;
  }
/*
 * Free up memory.
 */
  if(type_v_in  != NCL_double) NclFree(tmp_v_in);
  if(type_z_in  != NCL_double) NclFree(tmp_z_in);
  if(type_z_out != NCL_double) NclFree(tmp_z_out);
  if(type_v_out != NCL_double) NclFree(tmp_v_out);

/*
 * If v3d had a "description" or units attribute, return them with
 * the output  variable as an attribute.  Otherwise, return a
 * blank string for description, and nothing for units.
 */
  if(!found_desc) {
    cdesc = (char *)calloc(2,sizeof(char));
    strcpy(cdesc," ");
  }
/*
 * I don't think we can return "description" or "units" here, because 
 * they are attached to an NCL input parameter. It could screw things up
 * if we try to return it as an attribute with the output variable.
 * Instead, create a new description and units "quark" variable.
 */
  qdesc  = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *qdesc = NrmStringToQuark(cdesc);
  if (!found_desc)
          free(cdesc);
  if(found_units) {
    qunits  = (NclQuark*)NclMalloc(sizeof(NclQuark));
    *qunits = NrmStringToQuark(cunits);
  }
/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)v_out,
                            &missing_v_out,
                            ndims_z_out,
                            dsizes_v_out,
                            TEMPORARY,
                            NULL,
                            type_obj_v_out
                            );
/*
 * Set up attributes to return.
 */
  att_id = _NclAttCreate(NULL,NULL,Ncl_Att,0,NULL);

  dsizes[0] = 1;
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)qdesc,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "description",
             att_md,
             NULL
            );
    
  if(found_units) {
    att_md = _NclCreateVal(
                           NULL,
                           NULL,
                           Ncl_MultiDValData,
                           0,
                           (void*)qunits,
                           NULL,
                           1,
                           dsizes,
                           TEMPORARY,
                           NULL,
                           (NclObjClass)nclTypestringClass
                           );
    _NclAddAtt(
               att_id,
             "units",
               att_md,
               NULL
               );
  }
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          att_id,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );
  NclFree(dsizes_v_out);
  NclFree(dim_info);
/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);

}


NhlErrorTypes wrf_smooth_2d_W( void )
{
/*
 * Input variables
 *
 */
  void *a;
  int has_missing_a, ndims_a;
  ng_size_t dsizes_a[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_a;
  NclScalar missing_a;
  int *it;

/*
 * Various
 */
  double *db = NULL;
  float *fb = NULL;
  ng_size_t i, index_a, size_leftmost;
  int ny, nx, nynx;
  double d_missing;
  float  f_missing;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
/*
 * Get argument # 0
 */
  a = (void*)NclGetArgValue(
           0,
           2,
           &ndims_a,
           dsizes_a,
           &missing_a,
           &has_missing_a,
           &type_a,
           1);

/*
 * Check dimension sizes and input type.
 */
  if(ndims_a < 2) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_smooth_2d: The 'a' array must have at least 2 dimensions");
    return(NhlFATAL);
  }
  if(type_a != NCL_double && type_a != NCL_float) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_smooth_2d: The 'a' array must be float or double");
    return(NhlFATAL);
  }
  ny = dsizes_a[ndims_a-2];
  nx = dsizes_a[ndims_a-1];
  nynx = ny * nx;

/*
 * Get argument # 1
 */
  it = (int*)NclGetArgValue(
           1,
           2,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);

/*
 * Calculate size of leftmost dimensions.
 */
  size_leftmost  = 1;
  for(i = 0; i < ndims_a-2; i++) size_leftmost *= dsizes_a[i];

/*
 * Allocate space for "b", which "a" will be copied to inside
 * Fortran routine.
 */
  if(type_a == NCL_double) {
    db = (double *)malloc(nynx*sizeof(double));
    if(db == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_smooth_2d: Unable to allocate memory for temporary array");
      return(NhlFATAL);
    }
    d_missing = has_missing_a ? missing_a.doubleval : ((NclTypeClass)nclTypedoubleClass)->type_class.default_mis.doubleval;
  }
  else {
    fb = (float *)malloc(nynx*sizeof(float));
    if(fb == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_smooth_2d: Unable to allocate memory for temporary array");
      return(NhlFATAL);
    } 
    f_missing = has_missing_a ? missing_a.floatval : ((NclTypeClass)nclTypedoubleClass)->type_class.default_mis.floatval;
 }

/*
 * Loop across leftmost dimensions and call the Fortran routine for each
 * two-dimensional subsection.
 */
  index_a = 0;

  for(i = 0; i < size_leftmost; i++) {
    if(type_a == NCL_double) {
      NGCALLF(dfilter2d,DFILTER2D)(&((double*)a)[index_a], db, &nx, &ny, it,
                                   &d_missing);
    }
    else {
      NGCALLF(filter2d,FILTER2D)(&((float*)a)[index_a], fb, &nx, &ny, it,
                                   &f_missing);
    }
    index_a += nynx;
  }

  if(type_a == NCL_double) {
    NclFree(db);
  }
  else {
    NclFree(fb);
  }
/*
 * This is a procedure, so no values are returned.
 */
  return(NhlNOERROR);
}

NhlErrorTypes wrf_latlon_to_ij_W( void )
{

/*
 * Input variables
 */
  void *lat_array, *lon_array, *lat_loc, *lon_loc;
  double *tmp_lat_array = NULL; 
  double *tmp_lon_array = NULL; 
  double *tmp_lat_loc = NULL; 
  double *tmp_lon_loc = NULL; 
  int ndims_lat_array;
  ng_size_t dsizes_lat_array[NCL_MAX_DIMENSIONS];
  int ndims_lon_array;
  ng_size_t dsizes_lon_array[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_lat_loc[1];
  ng_size_t dsizes_lon_loc[1];
  NclBasicDataTypes type_lat_array, type_lon_array;
  NclBasicDataTypes type_lat_loc, type_lon_loc;
  int is_scalar_latlon_loc;

/*
 * Return variable
 */
  int iret, *ret;
  int ndims_ret;
  ng_size_t *dsizes_ret;
  NclScalar missing_ret;

/*
 * Various
 */
  ng_size_t ny, nx, nynx, nretlocs;
  ng_size_t index_array, index_ret;
  ng_size_t i, j, ndims_leftmost, size_leftmost, size_output;
  int inx, iny;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
/*
 * Get argument # 0
 */
  lat_array = (void*)NclGetArgValue(
           0,
           4,
           &ndims_lat_array,
           dsizes_lat_array,
           NULL,
           NULL,
           &type_lat_array,
           DONT_CARE);

/*
 * Get argument # 1
 */
  lon_array = (void*)NclGetArgValue(
           1,
           4,
           &ndims_lon_array,
           dsizes_lon_array,
           NULL,
           NULL,
           &type_lon_array,
           DONT_CARE);

/*
 * Check dimension sizes of lat,lon arrays and calculate size of
 * leftmost dimensions.
 */
  if(ndims_lat_array < 2 || ndims_lon_array < 2 || 
     ndims_lon_array != ndims_lat_array) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_latlon_to_ij: The lat,lon arrays must have at least two dimensions and the same number of dimensions as each other");
    return(NhlFATAL);
  }

  ny = dsizes_lat_array[ndims_lat_array-2];
  nx = dsizes_lat_array[ndims_lat_array-1];
  nynx = ny * nx;

/*
 * Test dimension sizes.
 */
  if((nx > INT_MAX) || (ny > INT_MAX)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_latlon_to_ij: nx and/or ny is greater than INT_MAX");
    return(NhlFATAL);
  }
  inx = (int) nx;
  iny = (int) ny;

  size_leftmost  = 1;
  ndims_leftmost = ndims_lat_array-2;
  for(i = 0; i < ndims_lon_array; i++) {
    if(dsizes_lon_array[i] != dsizes_lat_array[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_latlon_to_ij: The dimension sizes of the lat,lon arrays must be the same");
      return(NhlFATAL);
    }
    if(i < ndims_leftmost) size_leftmost *= dsizes_lat_array[i];
  }

/*
 * Get argument # 2
 */
  lat_loc = (void*)NclGetArgValue(
           2,
           4,
           NULL,
           dsizes_lat_loc,
           NULL,
           NULL,
           &type_lat_loc,
           DONT_CARE);

/*
 * Get argument # 3
 */
  lon_loc = (void*)NclGetArgValue(
           3,
           4,
           NULL,
           dsizes_lon_loc,
           NULL,
           NULL,
           &type_lon_loc,
           DONT_CARE);

/*
 * Check dimension sizes of lat,lon locations.
 */
  if(dsizes_lon_loc[0] != dsizes_lat_loc[0]) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_latlon_to_ij: The lat,lon locations must be the same length");
    return(NhlFATAL);
  }
  if(dsizes_lon_loc[0] == 1) {
    is_scalar_latlon_loc = 1;
  }
  else {
    is_scalar_latlon_loc = 0;
  }

/* 
 * Allocate space for coercing input arrays.  If any of the input
 * is already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 *
 * Allocate space for tmp_lat_array.
 */
  if(type_lat_array != NCL_double) {
    tmp_lat_array = (double *)calloc(nynx,sizeof(double));
    if(tmp_lat_array == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_latlon_to_ij: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
/*
 * Allocate space for tmp_lon_array.
 */
  if(type_lon_array != NCL_double) {
    tmp_lon_array = (double *)calloc(nynx,sizeof(double));
    if(tmp_lon_array == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_latlon_to_ij: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
/*
 * Allocate space for tmp_lat_loc.
 */
  if(type_lat_loc != NCL_double) {
    tmp_lat_loc = (double *)calloc(1,sizeof(double));
    if(tmp_lat_loc == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_latlon_to_ij: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
/*
 * Allocate space for tmp_lon_loc.
 */
  if(type_lon_loc != NCL_double) {
    tmp_lon_loc = (double *)calloc(1,sizeof(double));
    if(tmp_lon_loc == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_latlon_to_ij: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
/*
 * Calculate size of output array.  The output array will have dimension
 * sizes equal to the leftmost dimensions of the lat,lon arrays (minus
 * the last two dimensions), the length of the lat,lon locations
 * (if not a scalar), and the last dimension will be 2, which holds the
 * i,j location on the grid.
 */
  nretlocs    = size_leftmost * dsizes_lat_loc[0];
  size_output = 2 * nretlocs;

/* 
 * Allocate space for output array.
 */
  ret = (int*)calloc(size_output, sizeof(int));
  if(ret == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_latlon_to_ij: Unable to allocate memory for output array");
    return(NhlFATAL);
  }

/* 
 * Allocate space for output dimension sizes and set them. The last dimension
 * will always be 2, in order to hold the i,j locations. Some examples:
 *
 *  Lat,lon array are 90 x 180, lat,lon locations are scalars:
 *    Output will be array with 2 elements.
 *
 *   Lat,lon array are 5 x 90 x 180, lat,lon locations are scalars:
 *     Output will be array of length 5 x 2.
 *
 *   Lat,lon array are 5 x 90 x 180, lat,lon locations are length 10:
 *     Output will be array of length 5 x 10 x 2. 
 *
 *   Lat,lon array are 3 x 5 x 90 x 180, lat,lon locations are length 4:
 *     Output will be array of length 3 x 5 x 4 x 2.
 */
  if(is_scalar_latlon_loc) {
    ndims_ret = ndims_leftmost + 1;
  }
  else {
    ndims_ret = ndims_leftmost + 2;
  }
  dsizes_ret = (ng_size_t*)calloc(ndims_ret,sizeof(ng_size_t));  
  if( dsizes_ret == NULL ) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_latlon_to_ij: Unable to allocate memory for holding dimension sizes");
    return(NhlFATAL);
  }
/*
 * Fill in dimension sizes for output array.  See above examples.
 */
  for(i = 0; i < ndims_leftmost; i++) {
    dsizes_ret[i] = dsizes_lat_array[i];
  }
  if(!is_scalar_latlon_loc) {
    dsizes_ret[ndims_leftmost] = dsizes_lat_loc[0];
  }
  dsizes_ret[ndims_ret-1] = 2;

/*
 * Loop across leftmost dimensions of lat,lon array, the lat,lon locations,
 * and call the Fortran routine for each subsection of the input arrays.
 */
  index_array = index_ret = 0;
  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of lat_array (tmp_lat_array) to double if necessary.
 */
    if(type_lat_array != NCL_double) {
      coerce_subset_input_double(lat_array,tmp_lat_array,index_array,
                                 type_lat_array,nynx,0,NULL,NULL);
    }
    else {
      tmp_lat_array = &((double*)lat_array)[index_array];
    }
    
/*
 * Coerce subsection of lon_array (tmp_lon_array) to double if necessary.
 */
    if(type_lon_array != NCL_double) {
      coerce_subset_input_double(lon_array,tmp_lon_array,index_array,
                                 type_lon_array,nynx,0,NULL,NULL);
    }
    else {
      tmp_lon_array = &((double*)lon_array)[index_array];
    }

/*
 * Get default integer missing value.
 */
  missing_ret.intval =  ((NclTypeClass)nclTypeintClass)->type_class.default_mis.intval;

/*
 * Loop across lat,lon locations.
 */
    for(j = 0; j < dsizes_lat_loc[0]; j++) {

/*
 * Coerce subsection of lat_loc (tmp_lat_loc) to double if necessary.
 */
      if (type_lat_loc != NCL_double) {
        coerce_subset_input_double(lat_loc,tmp_lat_loc,j,type_lat_loc,1,0,
                                   NULL,NULL);
      }
      else {
        tmp_lat_loc = &((double*)lat_loc)[j];
      }

/*
 * Coerce subsection of lon_loc (tmp_lon_loc) to double if necessary.
 */
      if(type_lon_loc != NCL_double) {
        coerce_subset_input_double(lon_loc,tmp_lon_loc,j,type_lon_loc,1,0,
                                   NULL,NULL);
      }
      else {
        tmp_lon_loc = &((double*)lon_loc)[j];
      }

/*
 * Call the Fortran routine. Make sure you return the i,j index
 * swapped, since we are going from Fortran to C.
 */
      NGCALLF(dgetijlatlong,DGETIJLATLONG)(tmp_lat_array, tmp_lon_array, 
                                           tmp_lat_loc, tmp_lon_loc,
                                           &ret[index_ret+1], 
                                           &ret[index_ret], &inx, &iny,
                                           &missing_ret.intval);
      index_ret+=2;
    }
    index_array += nynx;
  }

/*
 * Free unneeded memory.
 */
  if(type_lat_array != NCL_double) NclFree(tmp_lat_array);
  if(type_lon_array != NCL_double) NclFree(tmp_lon_array);
  if(type_lat_loc   != NCL_double) NclFree(tmp_lat_loc);
  if(type_lon_loc   != NCL_double) NclFree(tmp_lon_loc);

  iret = NclReturnValue(ret,ndims_ret,dsizes_ret,&missing_ret,NCL_int,0);
  NclFree(dsizes_ret);
  return(iret);
}

NhlErrorTypes wrf_uvmet_W( void )
{

/*
 * Input variables
 */
/*
 * Argument # 0
 */
  void *u;
  double *tmp_u = NULL;
  int ndims_u;
  ng_size_t dsizes_u[NCL_MAX_DIMENSIONS];
  int has_missing_u;
  NclBasicDataTypes type_u;
  NclScalar missing_u, missing_du;

/*
 * Argument # 1
 */
  void *v;
  double *tmp_v = NULL;
  int ndims_v;
  ng_size_t dsizes_v[NCL_MAX_DIMENSIONS];
  int has_missing_v;
  NclBasicDataTypes type_v;
  NclScalar missing_v, missing_dv;

/*
 * Argument # 2
 */
  void *lat;
  double *tmp_lat = NULL;
  int ndims_lat;
  ng_size_t dsizes_lat[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_lat;

/*
 * Argument # 3
 */
  void *lon;
  double *tmp_lon = NULL;
  int ndims_lon;
  ng_size_t dsizes_lon[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_lon;

/*
 * Argument # 4
 */
  void *cenlon;
  double *tmp_cenlon;
  NclBasicDataTypes type_cenlon;

/*
 * Argument # 5
 */
  void *cone;
  double *tmp_cone;
  NclBasicDataTypes type_cone;

/*
 * Return variable and attributes.
 */
  void *uvmet;
  double *tmp_uvmet, tmp_uvmet_msg;
  int ndims_uvmet;
  ng_size_t *dsizes_uvmet;
  int has_missing;
  NclScalar missing_uvmet;
  NclBasicDataTypes type_uvmet;
  NclObjClass type_obj_uvmet;
  NclQuark *description, *units;
  char *cdescription, *cunits;
  NclDimRec *dim_info = NULL; 
  NclDimRec *dim_info_u, *dim_info_v;

/*
 * Various
 */
  ng_size_t nx, ny, nz, nxp1, nynxp1, nyp1, nyp1nx, nynx, twonynx;
  ng_size_t index_u, index_v, index_latlon, index_uvmet_u, index_uvmet_v;
  ng_size_t i, j;
  ng_size_t size_leftmost, size_leftmost_uvmet, size_uvmet, size_output;
  double rpd, *longca, *longcb;
  int istag, ndims_leftmost;
  int inx, iny, inxp1, inyp1;

/*
 * Variables for returning the output array with attributes attached.
 */
  int att_id;
  ng_size_t dsizes[1];
  NclMultiDValData att_md, return_md;
  NclVar tmp_var;
  NclStackEntry return_data;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
/*
 * Get argument # 0
 */
  u = (void*)NclGetArgValue(
           0,
           6,
           &ndims_u,
           dsizes_u,
           &missing_u,
           &has_missing_u,
           &type_u,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_u < 2) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_uvmet: The u array must have at least 2 dimensions");
    return(NhlFATAL);
  }
  ny     = dsizes_u[ndims_u-2];
  nxp1   = dsizes_u[ndims_u-1];
  nynxp1 = ny * nxp1;

/*
 * Get argument # 1
 */
  v = (void*)NclGetArgValue(
           1,
           6,
           &ndims_v,
           dsizes_v,
           &missing_v,
           &has_missing_v,
           &type_v,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_v != ndims_u) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_uvmet: The u and v arrays must have the same number of dimensions");
    return(NhlFATAL);
  }

  nyp1   = dsizes_v[ndims_v-2];
  nx     = dsizes_v[ndims_v-1];
  nyp1nx = nyp1 * nx;

/*
 * Test dimension sizes.
 */
  if((nxp1 > INT_MAX) || (nyp1 > INT_MAX) || (nx > INT_MAX) || (ny > INT_MAX)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_uvmet: one or more dimension sizes is greater than INT_MAX");
    return(NhlFATAL);
  }
  inx = (int) nx;
  iny = (int) ny;
  inxp1 = (int) nxp1;
  inyp1 = (int) nyp1;

/*
 * Coerce the missing values.
 */
  coerce_missing(type_u,has_missing_u,&missing_u,&missing_du,NULL);
  coerce_missing(type_v,has_missing_v,&missing_v,&missing_dv,NULL);
  if(has_missing_u || has_missing_v) {
    has_missing = True;
    /*fprintf(stderr, "\n\nfile: %s, line: %d\n", __FILE__, __LINE__);*/
    /* fprintf(stderr, "\tu or v has missing.\n");*/
  }
  else {
    has_missing = False;
  }
/*
 * Check whether we have staggered or unstaggered grids. 
 *
 * If unstaggered:
 *  - The rightmost two dimensions of u and v must be the same.
 *
 * If staggered:
 *  - The rightmost dimension of u must be one more than the
 *    rightmost dimension of v.
 *  - The second rightmost dimension of v must be one more 
 *    than the second rightmost dimension of u.
 */
  if(nxp1 == nx && nyp1 == ny) istag = 0;
  else                         istag = 1;

/*
 * Get argument # 2
 */
  lat = (void*)NclGetArgValue(
           2,
           6,
           &ndims_lat,
           dsizes_lat,
           NULL,
           NULL,
           &type_lat,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_lat != 2 && ndims_lat != ndims_u && 
     (ndims_u > 2 && ndims_lat != (ndims_u-1))) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_uvmet: The lat array must either be 2D, the same dimensions as u,v, or one fewer dimensions than u,v");
    return(NhlFATAL);
  }
  if(dsizes_lat[ndims_lat-2] != ny || dsizes_lat[ndims_lat-1] != nx) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_uvmet: The rightmost 2 dimensions of lat must be ny x nx");
    return(NhlFATAL);
  }
  nynx = ny * nx;

/*
 * Check dimension sizes for lat. It can be:
 *    - 2D (ny x nx)
 *    - Same dimensionality as U,V (but with rightmost dimemsions ny x nx)
 *    - One fewer dimension than U,V, with all leftmost up to the third
 *      rightmost dimensions the same as U,V.
 */
  if(ndims_lat > 2) {
    if(ndims_lat == ndims_u) {
      for(i = 0; i < ndims_u-2; i++) {
        if(dsizes_lat[i] != dsizes_u[i]) {
          NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_uvmet: if u and lat have the same number of dimensions, then all but the rightmost 2 dimensions must be the same");
          return(NhlFATAL);
        }
      }
    }
    else {
      for(i = 0; i < ndims_u-3; i++) {
        if(dsizes_lat[i] != dsizes_u[i]) {
          NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_uvmet: if lat has one fewer dimensions than u, then all but the rightmost 3 dimensions must be the same");
          return(NhlFATAL);
        }
      }
    }
  }

/*
 * Get argument # 3
 */
  lon = (void*)NclGetArgValue(
           3,
           6,
           &ndims_lon,
           dsizes_lon,
           NULL,
           NULL,
           &type_lon,
           DONT_CARE);

/*
 * Check dimension sizes for lon. This should be easier than lat,
 * since we've done all the work for lat, and the lat,lon have to be
 * exactly the same dimensions.
 */
  if(ndims_lon != ndims_lat) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_uvmet: The lat,lon arrays must have the same number of dimensions");
    return(NhlFATAL);
  }
  for(i = 0; i < ndims_lat; i++) {
    if(dsizes_lat[i] != dsizes_lon[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_uvmet: The lat,lon arrays must have the same dimension sizes");
      return(NhlFATAL);
    }
  }

/*
 * Get argument # 4
 */
  cenlon = (void*)NclGetArgValue(
           4,
           6,
           NULL,
           NULL,
           NULL,
           NULL,
           &type_cenlon,
           DONT_CARE);
/*
 * Get argument # 5
 */
  cone = (void*)NclGetArgValue(
           5,
           6,
           NULL,
           NULL,
           NULL,
           NULL,
           &type_cone,
           DONT_CARE);

/*
 * Calculate size of leftmost dimensions. Note that u, v can have an
 * extra leftmost dimension over lat, lon, so we need to separate these
 * out. The third-from-the-rightmost dimension will be called "nz".
 */
  size_leftmost  = 1;
  if(ndims_lat > 2 && ndims_lat == (ndims_u-1)) {
    nz = dsizes_u[ndims_u-3];
    ndims_leftmost = ndims_u-3;
    for(i = 0; i < ndims_leftmost; i++) size_leftmost *= dsizes_u[i];
  }
  else {
    nz = 1;
    ndims_leftmost = ndims_u-2;
    for(i = 0; i < ndims_leftmost; i++) size_leftmost *= dsizes_u[i];
  }

/*
 * The output type defaults to float, unless this input array is double.
 */
  type_uvmet     = NCL_float;
  type_obj_uvmet = nclTypefloatClass;
/* 
 * Allocate space for coercing input arrays.  If any of the input
 * is already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 *
 * Allocate space for tmp_u.
 */
  if(type_u != NCL_double) {
    tmp_u = (double *)calloc(nynxp1,sizeof(double));
    if(tmp_u == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_uvmet: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_uvmet     = NCL_double;
    type_obj_uvmet = nclTypedoubleClass;
  }
/*
 * Allocate space for tmp_v.
 */
  if(type_v != NCL_double) {
    tmp_v = (double *)calloc(nyp1nx,sizeof(double));
    if(tmp_v == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_uvmet: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_uvmet     = NCL_double;
    type_obj_uvmet = nclTypedoubleClass;
  }
/*
 * Allocate space for tmp_lat and tmp_lon, depending on whether
 * they are 2D or not.
 */
  if(ndims_lat == 2) {
    tmp_lat = coerce_input_double(lat,type_lat,nynx,0,NULL,NULL);
    tmp_lon = coerce_input_double(lon,type_lon,nynx,0,NULL,NULL);
  }
  else {
/*
 * Allocate space for tmp_lat
 */
    if(type_lat != NCL_double) {
      tmp_lat = (double *)calloc(nynx,sizeof(double));
      if(tmp_lat == NULL) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_uvmet: Unable to allocate memory for coercing input array to double");
        return(NhlFATAL);
      }
    }
/*
 * Allocate space for tmp_lon.
 */
    if(type_lon != NCL_double) {
      tmp_lon = (double *)calloc(nynx,sizeof(double));
      if(tmp_lon == NULL) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_uvmet: Unable to allocate memory for coercing input array to double");
        return(NhlFATAL);
      }
    }
  }

/*
 * Allocate space for tmp_cenlon.
 */
  tmp_cenlon = coerce_input_double(cenlon,type_cenlon,1,0,NULL,NULL);
  if(tmp_cenlon == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_uvmet: Unable to allocate memory for coercing input array to double");
    return(NhlFATAL);
  }
/*
 * Allocate space for tmp_cone.
 */
  tmp_cone = coerce_input_double(cone,type_cone,1,0,NULL,NULL);
  if(tmp_cone == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_uvmet: Unable to allocate memory for coercing input array to double");
    return(NhlFATAL);
  }

/*
 * Calculate size of output array.
 */
  twonynx             = 2 * ny * nx;
  size_leftmost_uvmet = size_leftmost * nz;
  size_uvmet          = size_leftmost_uvmet * nynx;
  size_output         = 2 * size_uvmet; 

/* 
 * Allocate space for output array.
 */
  tmp_uvmet = (double *)calloc(twonynx,sizeof(double));
  if(type_uvmet != NCL_double) {
    uvmet = (void *)calloc(size_output, sizeof(float));
    missing_uvmet.floatval = ((NclTypeClass)nclTypefloatClass)->type_class.default_mis.floatval;
    tmp_uvmet_msg = (double)missing_uvmet.floatval;
  }
  else {
    uvmet = (void *)calloc(size_output, sizeof(double));
    missing_uvmet.doubleval = ((NclTypeClass)nclTypedoubleClass)->type_class.default_mis.doubleval;
    tmp_uvmet_msg = missing_uvmet.doubleval;
  }
  if(uvmet == NULL || tmp_uvmet == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_uvmet: Unable to allocate memory for output array");
    return(NhlFATAL);
  }

/*
 * Allocate space for some dummy arrays.
 */
  longca = (double*)calloc(nynx,sizeof(double));
  longcb = (double*)calloc(nynx,sizeof(double));
  if( longca == NULL || longcb == NULL ) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_uvmet: Unable to allocate memory for output arrays");
    return(NhlFATAL);
  }

/* 
 * Allocate space for output dimension sizes and set them.
 */
  ndims_uvmet = ndims_u + 1;
  dsizes_uvmet = (ng_size_t*)calloc(ndims_uvmet,sizeof(ng_size_t));  
  if( dsizes_uvmet == NULL ) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_uvmet: Unable to allocate memory for holding dimension sizes");
    return(NhlFATAL);
  }
  dsizes_uvmet[0] = 2;
  for(i = 1; i < ndims_uvmet-2; i++) dsizes_uvmet[i] = dsizes_u[i-1];
  dsizes_uvmet[ndims_uvmet-2] = ny;
  dsizes_uvmet[ndims_uvmet-1] = nx;

/*
 * Loop across leftmost dimensions *and* nz, and call the Fortran
 * routine for each subsection of the input arrays.
 *
 * The input lat, lon arrays don't have the "nz" dimension,
 * so they have to be taken care of outside the nz loop.
 */
  index_u = index_v = index_latlon = index_uvmet_u = 0;
  index_uvmet_v = size_uvmet;

  rpd = 3.14159265/180.;

  for(i = 0; i < size_leftmost; i++) {
    if(ndims_lat > 2) {
/*
 * Coerce subsection of lat (tmp_lat) to double if necessary.
 */
      if(type_lat != NCL_double) {
        coerce_subset_input_double(lat,tmp_lat,index_latlon,type_lat,nynx,
                                   0,NULL,NULL);
      }
      else {
        tmp_lat = &((double*)lat)[index_latlon];
      }
/*
 * Coerce subsection of lon (tmp_lon) to double if necessary.
 */
      if(type_lon != NCL_double) {
        coerce_subset_input_double(lon,tmp_lon,index_latlon,type_lon,nynx,
                                   0,NULL,NULL);
      }
      else {
        tmp_lon = &((double*)lon)[index_latlon];
      }
    }
    for(j = 0; j < nz; j++) {
/*
 * Coerce subsection of u (tmp_u) to double if necessary.
 */
      if(type_u != NCL_double) {
        coerce_subset_input_double(u,tmp_u,index_u,type_u,nynxp1,0,NULL,NULL);
      }
      else {
        tmp_u = &((double*)u)[index_u];
      }
/*
 * Coerce subsection of v (tmp_v) to double if necessary.
 */
      if(type_v != NCL_double) {
        coerce_subset_input_double(v,tmp_v,index_v,type_v,nyp1nx,0,NULL,NULL);
      }
      else {
        tmp_v = &((double*)v)[index_v];
      }
/*
 * Call the Fortran routine.
 */
      NGCALLF(dcomputeuvmet,DCOMPUTEUVMET)(tmp_u, tmp_v, tmp_uvmet, longca, 
                                           longcb, tmp_lon, tmp_lat, 
                                           tmp_cenlon, tmp_cone, &rpd, 
                                           &inx, &iny, &inxp1, &inyp1, &istag,
                                           &has_missing,&missing_du.doubleval,
                                           &missing_du.doubleval,
                                           &tmp_uvmet_msg);

/*
 * Coerce output back to float if necessary.
 */
      coerce_output_float_or_double(uvmet,&tmp_uvmet[0],type_uvmet,nynx,
                                    index_uvmet_u);
      coerce_output_float_or_double(uvmet,&tmp_uvmet[nynx],type_uvmet,nynx,
                                    index_uvmet_v);

      index_u       += nynxp1;
      index_v       += nyp1nx;
      index_uvmet_u += nynx;
      index_uvmet_v += nynx;
    }
    if(ndims_lat > 2) {
      index_latlon += nynx;
    }
  }

/*
 * Free unneeded memory.
 */
  if(type_u      != NCL_double) NclFree(tmp_u);
  if(type_v      != NCL_double) NclFree(tmp_v);
  if(type_lat    != NCL_double) NclFree(tmp_lat);
  if(type_lon    != NCL_double) NclFree(tmp_lon);
  if(type_cenlon != NCL_double) NclFree(tmp_cenlon);
  if(type_cone   != NCL_double) NclFree(tmp_cone);
  NclFree(tmp_uvmet);
  NclFree(longca);
  NclFree(longcb);

/*
 * Set up some attributes ("description" and "units") to return.
 * Note that if the input arrays are anything but 2D, the units
 * will be "Temperature", and "2m Temperature" otherwise.
 */
  cdescription = (char *)calloc(17,sizeof(char));
  strcpy(cdescription,"u,v met velocity");
  cunits       = (char *)calloc(4,sizeof(char));
  strcpy(cunits,"m/s");
  description = (NclQuark*)NclMalloc(sizeof(NclQuark));
  units       = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *description = NrmStringToQuark(cdescription);
  *units       = NrmStringToQuark(cunits);
  free(cdescription);
  free(cunits);

/*
 * Get dimension info of U and V to see if we have named dimensions.
 * This will be used for return variable.  The return value's
 * dimension names will 
 */
  dim_info_u = get_wrf_dim_info(0,6,ndims_u,dsizes_u);
  dim_info_v = get_wrf_dim_info(1,6,ndims_v,dsizes_v);
  if(dim_info_u != NULL && dim_info_v != NULL) {
    dim_info = malloc(sizeof(NclDimRec)*ndims_uvmet);
    if(dim_info == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_uvmet: Unable to allocate memory for holding dimension information");
      return(NhlFATAL);
    }
    for(i = 0; i < ndims_uvmet; i++ ) {
      dim_info[i].dim_num  = i;
      dim_info[i].dim_size = dsizes_uvmet[i];
      if(i != 0) dim_info[i].dim_quark = dim_info_u[i-1].dim_quark;
      else       dim_info[0].dim_quark = NrmStringToQuark("u_v");
    }
/* 
 * Just the rightmost dimension is different from u's named dimensions. 
 */
    dim_info[ndims_uvmet-1].dim_quark = dim_info_v[ndims_v-1].dim_quark;
  }

/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)uvmet,
                            &missing_uvmet,
                            ndims_uvmet,
                            dsizes_uvmet,
                            TEMPORARY,
                            NULL,
                            type_obj_uvmet
                            );
/*
 * Set up attributes to return.
 */
  att_id = _NclAttCreate(NULL,NULL,Ncl_Att,0,NULL);

  dsizes[0] = 1;
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)description,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "description",
             att_md,
             NULL
             );
    
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)units,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "units",
             att_md,
             NULL
             );
    
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          att_id,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );
  NclFree(dsizes_uvmet);
  NclFree(dim_info);
  NclFree(dim_info_u);
  NclFree(dim_info_v);
/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);
}

NhlErrorTypes wrf_dbz_W( void )
{
/*
 * Input variables
 */
/*
 * Argument # 0
 */
  void *prs;
  double *tmp_prs = NULL;
  int ndims_prs;
  ng_size_t dsizes_prs[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_prs;

/*
 * Argument # 1
 */
  void *tmk;
  double *tmp_tmk = NULL;
  int ndims_tmk;
  ng_size_t dsizes_tmk[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_tmk;

/*
 * Argument # 2
 */
  void *qvp;
  double *tmp_qvp = NULL;
  int ndims_qvp;
  ng_size_t dsizes_qvp[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_qvp;

/*
 * Argument # 3
 */
  void *qra;
  double *tmp_qra = NULL;
  int ndims_qra;
  ng_size_t dsizes_qra[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_qra;

/*
 * Argument # 4
 */
  void *qsn;
  double *tmp_qsn;
  double *tmp1_qsn = NULL;
  int is_scalar_qsn, ndims_qsn;
  ng_size_t dsizes_qsn[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_qsn;

/*
 * Argument # 5
 */
  void *qgr;
  double *tmp_qgr = NULL;
  double *tmp1_qgr = NULL;
  int is_scalar_qgr, ndims_qgr;
  ng_size_t dsizes_qgr[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_qgr;

/*
 * Argument # 6
 */
  int *ivarint;
/*
 * Argument # 7
 */
  int *iliqskin;
/*
 * Return variable
 */
  void *dbz;
  double *tmp_dbz = NULL;
  NclBasicDataTypes type_dbz;
  NclObjClass type_obj_dbz;
  NclQuark *description, *units;
  char *cdescription, *cunits;
/*
 * Variables for returning the output array with dimension names attached.
 */
  int att_id;
  ng_size_t dsizes[1];
  NclMultiDValData return_md, att_md;
  NclVar tmp_var;
  NclStackEntry return_data;

/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info;

/*
 * Various
 */
  ng_size_t btdim, sndim, wedim, nbtsnwe, index_dbz;
  ng_size_t i, j, size_leftmost, size_output;
  int sn0 = 0, iwedim, isndim, ibtdim;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
/*
 * Get argument # 0
 */
  prs = (void*)NclGetArgValue(
           0,
           8,
           &ndims_prs,
           dsizes_prs,
           NULL,
           NULL,
           &type_prs,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_prs < 3) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_dbz: The prs array must have at least 3 dimensions");
    return(NhlFATAL);
  }
  btdim = dsizes_prs[ndims_prs-3];
  sndim = dsizes_prs[ndims_prs-2];
  wedim = dsizes_prs[ndims_prs-1];
  nbtsnwe = btdim * sndim * wedim;
  
/*
 * Test dimension sizes.
 */
  if((wedim > INT_MAX) || (sndim > INT_MAX) || (btdim > INT_MAX)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_dbz: one or more dimension sizes is greater than INT_MAX");    
    return(NhlFATAL);
  }
  iwedim = (int) wedim;
  isndim = (int) sndim;
  ibtdim = (int) btdim;

/*
 * Get argument # 1
 */
  tmk = (void*)NclGetArgValue(
           1,
           8,
           &ndims_tmk,
           dsizes_tmk,
           NULL,
           NULL,
           &type_tmk,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_tmk != ndims_prs) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_dbz: The tmk array must have the same number of dimensions as the prs array");
    return(NhlFATAL);
  }

/*
 * Get argument # 2
 */
  qvp = (void*)NclGetArgValue(
           2,
           8,
           &ndims_qvp,
           dsizes_qvp,
           NULL,
           NULL,
           &type_qvp,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_qvp != ndims_prs) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_dbz: The qv array must have the same number of dimensions as the prs array");
    return(NhlFATAL);
  }


/*
 * Get argument # 3
 */
  qra = (void*)NclGetArgValue(
           3,
           8,
           &ndims_qra,
           dsizes_qra,
           NULL,
           NULL,
           &type_qra,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_qra != ndims_prs) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_dbz: The qr array must have the same number of dimensions as the prs array");
    return(NhlFATAL);
  }


/*
 * Get argument # 4
 */
  qsn = (void*)NclGetArgValue(
           4,
           8,
           &ndims_qsn,
           dsizes_qsn,
           NULL,
           NULL,
           &type_qsn,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  is_scalar_qsn = is_scalar(ndims_qsn,dsizes_qsn);
  if(!is_scalar_qsn && ndims_qsn != ndims_prs) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_dbz: qs must either be a scalar or have the same number of dimensions as the prs array");
    return(NhlFATAL);
  }

/*
 * Get argument # 5
 */
  qgr = (void*)NclGetArgValue(
           5,
           8,
           &ndims_qgr,
           dsizes_qgr,
           NULL,
           NULL,
           &type_qgr,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  is_scalar_qgr = is_scalar(ndims_qgr,dsizes_qgr);
  if(!is_scalar_qgr && ndims_qgr != ndims_prs) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_dbz: qg must either be a scalar or have the same number of dimensions as the prs array");
    return(NhlFATAL);
  }

/*
 * Check that the first 6 input arrays all have the same dimensionality.
 */
  for(i = 0; i < ndims_prs; i++) {
    if(dsizes_tmk[i] != dsizes_prs[i] || dsizes_qvp[i] != dsizes_prs[i] || 
       dsizes_qra[i] != dsizes_prs[i] || 
      (!is_scalar_qsn && dsizes_qsn[i] != dsizes_prs[i]) || 
      (!is_scalar_qgr && dsizes_qgr[i] != dsizes_prs[i])) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_dbz: The prs, tmk, qv, qr, qs, and qg arrays must have the same dimensions (qs and qg can be scalars)");
      return(NhlFATAL);
    }
  }

/*
 * Get argument # 6
 */
  ivarint = (int*)NclGetArgValue(
           6,
           8,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);
/*
 * Get argument # 7
 */
  iliqskin = (int*)NclGetArgValue(
           7,
           8,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);

/*
 * Calculate size of leftmost dimensions.
 */
  size_leftmost  = 1;
  for(i = 0; i < ndims_prs-3; i++) size_leftmost *= dsizes_prs[i];

/*
 * The output type defaults to float, unless this input array is double.
 */
  type_dbz     = NCL_float;
  type_obj_dbz = nclTypefloatClass;

/* 
 * Allocate space for coercing input arrays.  If any of the input
 * is already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 */
/*
 * Allocate space for tmp_prs.
 */
  if(type_prs != NCL_double) {
    tmp_prs = (double *)calloc(nbtsnwe,sizeof(double));
    if(tmp_prs == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_dbz: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_dbz     = NCL_double;
    type_obj_dbz = nclTypedoubleClass;
  }
/*
 * Allocate space for tmp_tmk.
 */
  if(type_tmk != NCL_double) {
    tmp_tmk = (double *)calloc(nbtsnwe,sizeof(double));
    if(tmp_tmk == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_dbz: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_dbz     = NCL_double;
    type_obj_dbz = nclTypedoubleClass;
  }
/*
 * Allocate space for tmp_qvp.
 */
  if(type_qvp != NCL_double) {
    tmp_qvp = (double *)calloc(nbtsnwe,sizeof(double));
    if(tmp_qvp == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_dbz: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_dbz     = NCL_double;
    type_obj_dbz = nclTypedoubleClass;
  }
/*
 * Allocate space for tmp_qra no matter what, because qra might be
 * changed by the Fortran routine, and we don't want those changes
 * to propagate back here.
 */
  tmp_qra = (double *)calloc(nbtsnwe,sizeof(double));
  if(tmp_qra == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_dbz: Unable to allocate memory for coercing input array to double");
    return(NhlFATAL);
  }
  if(type_qra == NCL_double) {
    type_dbz     = NCL_double;
    type_obj_dbz = nclTypedoubleClass;
  }
/*
 * Allocate space for tmp_qsn no matter what, because qsn might be
 * changed by the Fortran routine, and we don't want those changes
 * to propagate back here.
 *
 * qsn could be a scalar. If so, we'll need to propagate it to a full
 * array. We'll do this later inside the do loop where the Fortran
 * routine is called.
 */
  tmp_qsn = (double *)calloc(nbtsnwe,sizeof(double));
  if(tmp_qsn == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_dbz: Unable to allocate memory for coercing input array to double");
    return(NhlFATAL);
  }
  if(is_scalar_qsn) {
    if(type_qsn != NCL_double) {
      tmp1_qsn = coerce_input_double(qsn,type_qsn,1,0,NULL,NULL);
    }
    else {
      tmp1_qsn  = (double*)malloc(sizeof(double));
      *tmp1_qsn = ((double*)qsn)[0];
    }
  }
  if(type_qsn == NCL_double) {
    type_dbz     = NCL_double;
    type_obj_dbz = nclTypedoubleClass;
  }
/*
 * Allocate space for tmp_qgr.
 *
 * If it is a scalar, then propagate the scalar to an array.
 */
  if(is_scalar_qgr || type_qgr != NCL_double) { 
    tmp_qgr = (double *)calloc(nbtsnwe,sizeof(double));
    if(tmp_qgr == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_dbz: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  if(is_scalar_qgr) {
    tmp1_qgr = coerce_input_double(qgr,type_qgr,1,0,NULL,NULL);
    for(i = 0; i < nbtsnwe; i++) tmp_qgr[i] = *tmp1_qgr;
  }
  if(type_qgr == NCL_double) {
    type_dbz     = NCL_double;
    type_obj_dbz = nclTypedoubleClass;
  }

/* 
 * Allocate space for output array.
 */
  size_output = size_leftmost * nbtsnwe;

  if(type_dbz != NCL_double) {
    dbz = (void *)calloc(size_output, sizeof(float));
    tmp_dbz = (double *)calloc(nbtsnwe,sizeof(double));
    if(tmp_dbz == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_dbz: Unable to allocate memory for temporary output array");
      return(NhlFATAL);
    }
  }
  else {
    dbz = (void *)calloc(size_output, sizeof(double));
  }
  if(dbz == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_dbz: Unable to allocate memory for output array");
    return(NhlFATAL);
  }

/*
 * Loop across leftmost dimensions and call the Fortran routine for each
 * subsection of the input arrays.
 */
  index_dbz = 0;

  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of prs (tmp_prs) to double if necessary.
 */
    if(type_prs != NCL_double) {
      coerce_subset_input_double(prs,tmp_prs,index_dbz,type_prs,nbtsnwe,
                                 0,NULL,NULL);
    }
    else {
      tmp_prs = &((double*)prs)[index_dbz];
    }

/*
 * Coerce subsection of tmk (tmp_tmk) to double if necessary.
 */
    if(type_tmk != NCL_double) {
      coerce_subset_input_double(tmk,tmp_tmk,index_dbz,type_tmk,nbtsnwe,
                                 0,NULL,NULL);
    }
    else {
      tmp_tmk = &((double*)tmk)[index_dbz];
    }

/*
 * Coerce subsection of qvp (tmp_qvp) to double if necessary.
 */
    if(type_qvp != NCL_double) {
      coerce_subset_input_double(qvp,tmp_qvp,index_dbz,type_qvp,nbtsnwe,
                                 0,NULL,NULL);
    }
    else {
      tmp_qvp = &((double*)qvp)[index_dbz];
    }

/*
 * If qsn is a scalar, then propagate it to a full array.
 */
    if(is_scalar_qsn) {
      for(j = 0; j < nbtsnwe; j++) tmp_qsn[j] = *tmp1_qsn;
      if(*tmp1_qsn == 0.) {
        sn0 = 0;
      }
    }
    else {
/*
 * Force the coercion of qsn to tmp_qsn, because the original arrays may
 * get changed by the Fortran routine, and we don't want those changes to
 * propagate back here.
 */
      coerce_subset_input_double(qsn,tmp_qsn,index_dbz,type_qsn,nbtsnwe,
                                 0,NULL,NULL);
/*
 * Check values for qsn array. If all zero, then set sn0 to 0. Otherwise
 * set sn0 to 1.
 */
      j   = 0;
      sn0 = 0;
      while( (j < nbtsnwe) && !sn0) {
        if(tmp_qsn[j] != 0.) sn0 = 1;
        j++;
      }
    }

/*
 * Force the coercion of qra to tmp_qra, because the original arrays may
 * get changed by the Fortran routine, and we don't want those changes to
 * propagate back here.
 */
    coerce_subset_input_double(qra,tmp_qra,index_dbz,type_qra,nbtsnwe,
                               0,NULL,NULL);
/*
 * Coerce subsection of qgr (tmp_qgr) to double if necessary.
 */
    if(!is_scalar_qgr) {
      double *tmp_qgr_save = tmp_qgr;
      if(type_qgr != NCL_double) {
        coerce_subset_input_double(qgr,tmp_qgr,index_dbz,type_qgr,nbtsnwe,
                                   0,NULL,NULL);
      }
      else {
        tmp_qgr = &((double*)qgr)[index_dbz];
      }
      if (tmp_qgr_save != NULL && tmp_qgr_save != tmp_qgr)
              NclFree(tmp_qgr);
    }

/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_dbz == NCL_double) tmp_dbz = &((double*)dbz)[index_dbz];

/*
 * Call the Fortran routine.
 */
    NGCALLF(calcdbz,CALCDBZ)(tmp_dbz, tmp_prs, tmp_tmk, tmp_qvp, tmp_qra, 
                             tmp_qsn, tmp_qgr, &iwedim, &isndim, &ibtdim, 
                             &sn0, ivarint, iliqskin);
/*
 * Coerce output back to float if necessary.
 */
    if(type_dbz == NCL_float) {
      coerce_output_float_only(dbz,tmp_dbz,nbtsnwe,index_dbz);
    }
    index_dbz += nbtsnwe;
  }

/*
 * Free unneeded memory.
 */
  if(type_prs != NCL_double) NclFree(tmp_prs);
  if(type_tmk != NCL_double) NclFree(tmp_tmk);
  if(type_qvp != NCL_double) NclFree(tmp_qvp);
  NclFree(tmp_qra);
  NclFree(tmp_qsn);
  if(type_qgr != NCL_double)  NclFree(tmp_qgr);
  if(type_dbz != NCL_double) NclFree(tmp_dbz);
  if(is_scalar_qsn) NclFree(tmp1_qsn);
  if(is_scalar_qgr && type_qgr != NCL_double) NclFree(tmp1_qgr);

/*
 * Retrieve dimension names from the "tmk" variable, if any.
 * These dimension names will later be attached to the output variable.
 */
  dim_info = get_wrf_dim_info(1,8,ndims_tmk,dsizes_tmk);

/*
 * Set up return value.
 */
/*
 * Return value back to NCL script.
 */

  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)dbz,
                            NULL,
                            ndims_tmk,
                            dsizes_tmk,
                            TEMPORARY,
                            NULL,
                            type_obj_dbz
                            );
/*
 * Set up some attributes ("description" and "units") to return.
 */
  cdescription = (char *)calloc(13,sizeof(char));
  strcpy(cdescription,"Reflectivity");
  description  = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *description = NrmStringToQuark(cdescription);
  free(cdescription);

  cunits       = (char *)calloc(4,sizeof(char));
  strcpy(cunits,"dBZ");
  units        = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *units       = NrmStringToQuark(cunits);
  free(cunits);


/*
 * Set up attributes to return.
 */
  att_id = _NclAttCreate(NULL,NULL,Ncl_Att,0,NULL);

  dsizes[0] = 1;
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)description,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "description",
             att_md,
             NULL
             );
    
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)units,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "units",
             att_md,
             NULL
             );
    
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          att_id,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );
  NclFree(dim_info);

/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);

}

NhlErrorTypes wrf_pvo_W( void )
{

/*
 * Input variables
 *
 * Argument # 0
 */
  void *u;
  double *tmp_u = NULL;
  int ndims_u;
  ng_size_t dsizes_u[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_u;

/*
 * Argument # 1
 */
  void *v;
  double *tmp_v = NULL;
  int ndims_v;
  ng_size_t dsizes_v[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_v;

/*
 * Argument # 2
 */
  void *th;
  double *tmp_th = NULL;
  int ndims_th;
  ng_size_t dsizes_th[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_th;

/*
 * Argument # 3
 */
  void *p;
  double *tmp_p = NULL;
  int ndims_p;
  ng_size_t dsizes_p[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_p;

/*
 * Argument # 4
 */
  void *msfu;
  double *tmp_msfu = NULL;
  int ndims_msfu;
  ng_size_t dsizes_msfu[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_msfu;

/*
 * Argument # 5
 */
  void *msfv;
  double *tmp_msfv = NULL;
  int ndims_msfv;
  ng_size_t dsizes_msfv[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_msfv;

/*
 * Argument # 6
 */
  void *msft;
  double *tmp_msft = NULL;
  int ndims_msft;
  ng_size_t dsizes_msft[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_msft;

/*
 * Argument # 7
 */
  void *cor;
  double *tmp_cor = NULL;
  int ndims_cor;
  ng_size_t dsizes_cor[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_cor;

/*
 * Argument # 8
 */
  void *dx;
  double *tmp_dx = NULL;
  NclBasicDataTypes type_dx;

/*
 * Argument # 9
 */
  void *dy;
  double *tmp_dy = NULL;
  NclBasicDataTypes type_dy;

/*
 * Argument # 10
 */
  int *opt;

/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info;

/*
 * Return variable
 */
  void *pv;
  double *tmp_pv = NULL;
  int att_id;
  NclBasicDataTypes type_pv;
  NclObjClass type_obj_pv;
  NclQuark *description, *units;
  char *cdescription, *cunits;

/*
 * Various
 */
  ng_size_t nx, ny, nz, nxp1, nyp1;
  ng_size_t nznynxp1, nznyp1nx, nznynx, nynxp1, nyp1nx, nynx;
  ng_size_t i, size_pv, size_leftmost;
  ng_size_t index_u, index_v, index_th, index_msfu, index_msfv, index_msft;
  int inx, iny, inz, inxp1, inyp1;
/*
 * Variables for returning the output array with dimension names attached.
 */
  ng_size_t dsizes[1];
  NclMultiDValData att_md, return_md;
  NclVar tmp_var;
  NclStackEntry return_data;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
/*
 * Get argument # 0
 */
  u = (void*)NclGetArgValue(
           0,
           11,
           &ndims_u,
           dsizes_u,
           NULL,
           NULL,
           &type_u,
           DONT_CARE);

/*
 * Error checking on dimensions.
 */
  if(ndims_u < 3) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: u must have at least 3 dimensions");
    return(NhlFATAL);
  }
  nz   = dsizes_u[ndims_u-3];
  ny   = dsizes_u[ndims_u-2];
  nxp1 = dsizes_u[ndims_u-1];

/*
 * Get argument # 1
 */
  v = (void*)NclGetArgValue(
           1,
           11,
           &ndims_v,
           dsizes_v,
           NULL,
           NULL,
           &type_v,
           DONT_CARE);

/*
 * Error checking on dimensions.
 */
  if(ndims_v != ndims_u) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: u, v, th, and p must have the same number of dimensions");
    return(NhlFATAL);
  }
  if(dsizes_v[ndims_v-3] != nz) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: The third-from-the-right dimension of v must be the same as the third-from-the-right dimension of u");
    return(NhlFATAL);
  }
/*
 * Error checking on leftmost dimension sizes.
 */
  for(i = 0; i < ndims_u-3; i++) {
    if(dsizes_u[i] != dsizes_v[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: The leftmost dimensions of u and v must be the same");
      return(NhlFATAL);
    }
  }

  nyp1 = dsizes_v[ndims_v-2];
  nx   = dsizes_v[ndims_v-1];

/*
 * Get argument # 2
 */
  th = (void*)NclGetArgValue(
           2,
           11,
           &ndims_th,
           dsizes_th,
           NULL,
           NULL,
           &type_th,
           DONT_CARE);

/*
 * Error checking on dimensions.
 */
  if(ndims_th != ndims_u) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: u, v, th, and p must have the same number of dimensions");
    return(NhlFATAL);
  }

  if(dsizes_th[ndims_th-3] != nz || dsizes_th[ndims_th-2] != ny ||
     dsizes_th[ndims_th-1] != nx) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: The rightmost dimensions of th must be a combination of the dimensions of u and v (see documentation)");
    return(NhlFATAL);
  }

/*
 * Error checking on leftmost dimension sizes.
 */
  for(i = 0; i < ndims_u-3; i++) {
    if(dsizes_th[i] != dsizes_u[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: The leftmost dimensions of th and u must be the same");
      return(NhlFATAL);
    }
  }

/*
 * Get argument # 3
 */
  p = (void*)NclGetArgValue(
           3,
           11,
           &ndims_p,
           dsizes_p,
           NULL,
           NULL,
           &type_p,
           DONT_CARE);

/*
 * Error checking on dimensions.
 */
  if(ndims_p != ndims_u) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: u, v, th, and p must have the same number of dimensions");
    return(NhlFATAL);
  }

/*
 * Error checking on dimension sizes.
 */
  for(i = 0; i < ndims_th; i++) {
    if(dsizes_p[i] != dsizes_th[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: The dimensions of p and th must be the same");
      return(NhlFATAL);
    }
  }

/*
 * Get argument # 4
 */
  msfu = (void*)NclGetArgValue(
           4,
           11,
           &ndims_msfu,
           dsizes_msfu,
           NULL,
           NULL,
           &type_msfu,
           DONT_CARE);

/*
 * Error checking on dimensions.
 */
  if(ndims_msfu < 2) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: msfu must have at least 2 dimensions");
    return(NhlFATAL);
  }
  if(ndims_msfu !=2 && ndims_msfu != (ndims_u-1)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: msfu must be 2D or have one fewer dimensions than u");
    return(NhlFATAL);
  }
  if(dsizes_msfu[ndims_msfu-2] != ny || dsizes_msfu[ndims_msfu-1] != nxp1) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: The rightmost 2 dimensions of msfu must be the same as the rightmost 2 dimensions of u");
    return(NhlFATAL);
  }

/*
 * Error checking on leftmost dimension sizes. msfu, msfv, msft, and 
 * cor can be 2D or nD.  If they are nD, they must have same leftmost
 * dimensions as other input arrays.
 */
  if(ndims_msfu > 2) {
    for(i = 0; i < ndims_u-3; i++) {
      if(dsizes_msfu[i] != dsizes_u[i]) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: If msfu is not 2-dimensional, then the leftmost dimensions of msfu and u must be the same");
        return(NhlFATAL);
      }
    }
  }

/*
 * Get argument # 5
 */
  msfv = (void*)NclGetArgValue(
           5,
           11,
           &ndims_msfv,
           dsizes_msfv,
           NULL,
           NULL,
           &type_msfv,
           DONT_CARE);

/*
 * Error checking on dimensions.
 */
  if(ndims_msfv != ndims_msfu) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: msfu, msfv, msft, and cor must have the same number of dimensions");
    return(NhlFATAL);
  }
  if(dsizes_msfv[ndims_msfv-2] != nyp1 || dsizes_msfv[ndims_msfv-1] != nx) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: The rightmost 2 dimensions of msfv must be the same as the rightmost 2 dimensions of v");
    return(NhlFATAL);
  }

/*
 * Error checking on leftmost dimension sizes.
 */
  for(i = 0; i < ndims_msfu-2; i++) {
    if(dsizes_msfv[i] != dsizes_msfu[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: The leftmost dimensions of msfv and msfu must be the same");
      return(NhlFATAL);
    }
  }

/*
 * Get argument # 6
 */
  msft = (void*)NclGetArgValue(
           6,
           11,
           &ndims_msft,
           dsizes_msft,
           NULL,
           NULL,
           &type_msft,
           DONT_CARE);

/*
 * Error checking on dimensions.
 */
  if(ndims_msft != ndims_msfu) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: msfu, msfv, msft, and cor must have the same number of dimensions");
    return(NhlFATAL);
  }
  if(dsizes_msft[ndims_msft-2] != ny || dsizes_msft[ndims_msft-1] != nx) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: The rightmost 2 dimensions of msft must be the same as the rightmost 2 dimensions of th");
    return(NhlFATAL);
  }

/*
 * Error checking on leftmost dimension sizes.
 */
  for(i = 0; i < ndims_msfu-2; i++) {
    if(dsizes_msft[i] != dsizes_msfu[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: The leftmost dimensions of msft and msfu must be the same");
      return(NhlFATAL);
    }
  }

/*
 * Get argument # 7
 */
  cor = (void*)NclGetArgValue(
           7,
           11,
           &ndims_cor,
           dsizes_cor,
           NULL,
           NULL,
           &type_cor,
           DONT_CARE);

/*
 * Error checking on dimensions.
 */
  if(ndims_cor != ndims_msft) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: msfu, msfv, msft, and cor must have the same number of dimensions");
    return(NhlFATAL);
  }

/*
 * Error checking on dimension sizes.
 */
  for(i = 0; i < ndims_msft; i++) {
    if(dsizes_cor[i] != dsizes_msft[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: The dimensions of cor and msft must be the same");
      return(NhlFATAL);
    }
  }

/*
 * Get argument # 8
 */
  dx = (void*)NclGetArgValue(
           8,
           11,
           NULL,
           NULL,
           NULL,
           NULL,
           &type_dx,
           DONT_CARE);
  tmp_dx = coerce_input_double(dx,type_dx,1,0,NULL,NULL);

/*
 * Get argument # 9
 */
  dy = (void*)NclGetArgValue(
           9,
           11,
           NULL,
           NULL,
           NULL,
           NULL,
           &type_dy,
           DONT_CARE);
  tmp_dy = coerce_input_double(dy,type_dy,1,0,NULL,NULL);

/*
 * Get argument # 10
 */
  opt = (int*)NclGetArgValue(
           10,
           11,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);

  nynx     = ny * nx;
  nznynx   = nz * nynx;
  nynxp1   = ny * nxp1;
  nyp1nx   = nyp1 * nx;
  nznynxp1 = nz * nynxp1;
  nznyp1nx = nz * nyp1nx;

/*
 * Test dimension sizes.
 */
    if((nxp1 > INT_MAX) || (nyp1 > INT_MAX) || (nz > INT_MAX) || 
       (nx > INT_MAX) ||(ny > INT_MAX)) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: one or more dimension sizes is greater than INT_MAX");
      return(NhlFATAL);
    }
    inx = (int) nx;
    iny = (int) ny;
    inz = (int) nz;
    inxp1 = (int) nxp1;
    inyp1 = (int) nyp1;

/*
 * Calculate size of leftmost dimensions.
 */
  size_leftmost = 1;
  for(i = 0; i < ndims_u-3; i++) size_leftmost *= dsizes_u[i];
  size_pv = size_leftmost * nznynx;

/*
 * Retrieve dimension names from the "th" variable, if any.
 * These dimension names will later be attached to the output variable.
 */
  dim_info = get_wrf_dim_info(2,11,ndims_th,dsizes_th);

/* 
 * Allocate space for coercing input arrays.  If any of the input
 * is already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 */

/*
 * Allocate space for tmp_u.
 */
  if(type_u != NCL_double) {
    tmp_u = (double *)calloc(nznynxp1,sizeof(double));
    if(tmp_u == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: Unable to allocate memory for coercing u to double");
      return(NhlFATAL);
    }
  }

/*
 * Allocate space for tmp_v.
 */
  if(type_v != NCL_double) {
    tmp_v = (double *)calloc(nznyp1nx,sizeof(double));
    if(tmp_v == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: Unable to allocate memory for coercing v to double");
      return(NhlFATAL);
    }
  }
/*
 * Allocate space for tmp_th.
 */
  if(type_th != NCL_double) {
    tmp_th = (double *)calloc(nznynx,sizeof(double));
    if(tmp_th == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: Unable to allocate memory for coercing th to double");
      return(NhlFATAL);
    }
  }
/*
 * Allocate space for tmp_p.
 */
  if(type_p != NCL_double) {
    tmp_p = (double *)calloc(nznynx,sizeof(double));
    if(tmp_p == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: Unable to allocate memory for coercing p to double");
      return(NhlFATAL);
    }
  }
/*
 * Allocate space for coercing msfu, msfv, and cor to double precision.
 * These arrays can be 2D or nD, so take this into account. If one of
 * them is 2D, then all three of them have to be 2D.
 */
  if(ndims_msfu == 2) {
    tmp_msfu = coerce_input_double(msfu,type_msfu,nynxp1,0,NULL,NULL);
    tmp_msfv = coerce_input_double(msfv,type_msfv,nyp1nx,0,NULL,NULL);
    tmp_msft = coerce_input_double(msft,type_msft,nynx,0,NULL,NULL);
    tmp_cor  = coerce_input_double(cor,type_cor,nynx,0,NULL,NULL);
  }
  else {
/*
 * Allocate space for tmp_msfu.
 */
    if(type_msfu != NCL_double) {
      tmp_msfu = (double*)calloc(nynxp1,sizeof(double));
      if(tmp_msfu == NULL) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: Unable to allocate memory for coercing msfu to double");
        return(NhlFATAL);
      }
    }
/*
 * Allocate space for tmp_msfv.
 */
    if(type_msfv != NCL_double) {
      tmp_msfv = (double*)calloc(nyp1nx,sizeof(double));
      if(tmp_msfv == NULL) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: Unable to allocate memory for coercing msfv to double");
        return(NhlFATAL);
      }
    }
/*
 * Allocate space for tmp_msft.
 */
    if(type_msft != NCL_double) {
      tmp_msft = (double*)calloc(nynx,sizeof(double));
      if(tmp_msft == NULL) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: Unable to allocate memory for coercing msft to double");
        return(NhlFATAL);
      }
    }
/*
 * Allocate space for tmp_cor.
 */
    if(type_cor != NCL_double) {
      tmp_cor = (double *)calloc(nynx,sizeof(double));
      if(tmp_cor == NULL) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: Unable to allocate memory for coercing cor to double");
        return(NhlFATAL);
      }
    }
  }

/*
 * The output type defaults to float, unless any input arrays are double.
 */
  if(type_u    == NCL_double || type_v    == NCL_double || 
     type_th   == NCL_double || type_p    == NCL_double || 
     type_msfu == NCL_double || type_msfv == NCL_double || 
     type_msft == NCL_double || type_cor  == NCL_double) {
    type_pv     = NCL_double;
    type_obj_pv = nclTypedoubleClass;
  }
  else { 
    type_pv     = NCL_float;
    type_obj_pv = nclTypefloatClass;
  }

/* 
 * Allocate space for output array.
 */
  if(type_pv != NCL_double) {
    pv     = (void *)calloc(size_pv, sizeof(float));
    tmp_pv = (double *)calloc(nznynx,sizeof(double));
    if(pv == NULL || tmp_pv == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }
  else {
    pv = (void *)calloc(size_pv, sizeof(double));
    if(pv == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_pvo: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }

/*
 * Call the Fortran routine.
 */
  index_u = index_v = index_th = index_msfu = index_msfv = index_msft = 0;
  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of u (tmp_u) to double if necessary.
 */
    if(type_u != NCL_double) {
      coerce_subset_input_double(u,tmp_u,index_u,type_u,nznynxp1,0,NULL,NULL);
    }
    else {
      tmp_u = &((double*)u)[index_u];
    }

/*
 * Coerce subsection of v (tmp_v) to double if necessary.
 */
    if(type_v != NCL_double) {
      coerce_subset_input_double(v,tmp_v,index_v,type_v,nznyp1nx,0,NULL,NULL);
    }
    else {
      tmp_v = &((double*)v)[index_v];
    }

/*
 * Coerce subsection of th (tmp_th) to double if necessary.
 */
    if(type_th != NCL_double) {
      coerce_subset_input_double(th,tmp_th,index_th,type_th,nznynx,0,NULL,NULL);
    }
    else {
      tmp_th = &((double*)th)[index_th];
    }

/*
 * Coerce subsection of p (tmp_p) to double if necessary.
 */
    if(type_p != NCL_double) {
      coerce_subset_input_double(p,tmp_p,index_th,type_p,nznynx,0,NULL,NULL);
    }
    else {
      tmp_p = &((double*)p)[index_th];
    }

/*
 * msfu, msfv, msft, and cor can be 2D or nD, so account
 * for that here. If they are 2D, they've already been coerced
 * before the loop.
 */
    if(ndims_msfu > 2) {
/*
 * Coerce subsection of msfu (tmp_msfu) to double if necessary.
 */
      if(type_msfu != NCL_double) {
        coerce_subset_input_double(msfu,tmp_msfu,index_msfu,type_msfu,nynxp1,0,NULL,NULL);
      }
      else {
        tmp_msfu = &((double*)msfu)[index_msfu];
      }
/*
 * Coerce subsection of msfv (tmp_msfv) to double if necessary.
 */
      if(type_msfv != NCL_double) {
        coerce_subset_input_double(msfv,tmp_msfv,index_msfv,type_msfv,nyp1nx,0,NULL,NULL);
      }
      else {
        tmp_msfv = &((double*)msfv)[index_msfv];
      }
/*
 * Coerce subsection of msft (tmp_msft) to double if necessary.
 */
      if(type_msft != NCL_double) {
        coerce_subset_input_double(msft,tmp_msft,index_msft,type_msft,nynx,0,NULL,NULL);
      }
      else {
        tmp_msft = &((double*)msft)[index_msft];
      }

/*
 * Coerce subsection of cor (tmp_cor) to double if necessary.
 */
      if(type_cor != NCL_double) {
        coerce_subset_input_double(cor,tmp_cor,index_msft,type_cor,nynx,0,NULL,NULL);
      }
      else {
        tmp_cor = &((double*)cor)[index_msft];
      }
    }

/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_pv == NCL_double) tmp_pv = &((double*)pv)[index_th];

      NGCALLF(dcomputepv,DCOMPUTEPV)(tmp_pv, tmp_u, tmp_v, tmp_th, tmp_p, 
                                     tmp_msfu, tmp_msfv, tmp_msft, tmp_cor, 
                                     tmp_dx, tmp_dy, &inx, &iny, &inz, &inxp1, &inyp1);
    if(type_pv != NCL_double) {
      coerce_output_float_only(pv,tmp_pv,nznynx,index_th);
    }
    index_u    += nznynxp1;
    index_v    += nznyp1nx;
    index_th   += nznynx;
    if(ndims_msfu > 2) {
      index_msfu += nynxp1;
      index_msfv += nyp1nx;
      index_msft += nynx;
    }
  }

/*
 * Free unneeded memory.
 */
  if(type_u    != NCL_double) NclFree(tmp_u);
  if(type_v    != NCL_double) NclFree(tmp_v);
  if(type_th   != NCL_double) NclFree(tmp_th);
  if(type_p    != NCL_double) NclFree(tmp_p);
  if(type_msfu != NCL_double) NclFree(tmp_msfu);
  if(type_msfv != NCL_double) NclFree(tmp_msfv);
  if(type_msft != NCL_double) NclFree(tmp_msft);
  if(type_cor  != NCL_double) NclFree(tmp_cor);
  if(type_pv   != NCL_double) NclFree(tmp_pv);

/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)pv,
                            NULL,
                            ndims_th,
                            dsizes_th,
                            TEMPORARY,
                            NULL,
                            type_obj_pv
                            );

/*
 * Set up some attributes ("description" and "units") to return.
 */
  cdescription = (char *)calloc(20,sizeof(char));
  strcpy(cdescription,"Potential Vorticity");
  description  = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *description = NrmStringToQuark(cdescription);
  free(cdescription);

  cunits       = (char *)calloc(4,sizeof(char));
  strcpy(cunits,"PVU");
  units        = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *units       = NrmStringToQuark(cunits);
  free(cunits);

/*
 * Set up attributes to return.
 */
  att_id = _NclAttCreate(NULL,NULL,Ncl_Att,0,NULL);

  dsizes[0] = 1;
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)description,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "description",
             att_md,
             NULL
             );
    
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)units,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "units",
             att_md,
             NULL
             );
    
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          att_id,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );

  NclFree(dim_info);
/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);
}

NhlErrorTypes wrf_avo_W( void )
{

/*
 * Input variables
 *
 * Argument # 0
 */
  void *u;
  double *tmp_u = NULL;
  int ndims_u;
  ng_size_t dsizes_u[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_u;

/*
 * Argument # 1
 */
  void *v;
  double *tmp_v = NULL;
  int ndims_v;
  ng_size_t dsizes_v[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_v;

/*
 * Argument # 2
 */
  void *msfu;
  double *tmp_msfu = NULL;
  int ndims_msfu;
  ng_size_t dsizes_msfu[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_msfu;

/*
 * Argument # 3
 */
  void *msfv;
  double *tmp_msfv = NULL;
  int ndims_msfv;
  ng_size_t dsizes_msfv[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_msfv;

/*
 * Argument # 4
 */
  void *msft;
  double *tmp_msft = NULL;
  int ndims_msft;
  ng_size_t dsizes_msft[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_msft;

/*
 * Argument # 5
 */
  void *cor;
  double *tmp_cor = NULL;
  int ndims_cor;
  ng_size_t dsizes_cor[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_cor;

/*
 * Argument # 6
 */
  void *dx;
  double *tmp_dx = NULL;
  NclBasicDataTypes type_dx;

/*
 * Argument # 7
 */
  void *dy;
  double *tmp_dy = NULL;
  NclBasicDataTypes type_dy;

/*
 * Argument # 8
 */
  int *opt;

/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info, *dim_info_v;

/*
 * Return variable
 */
  void *av;
  double *tmp_av = NULL;
  int att_id;
  ng_size_t *dsizes_av;
  NclBasicDataTypes type_av;
  NclObjClass type_obj_av;
  NclQuark *description, *units;
  char *cdescription, *cunits;

/*
 * Various
 */
  ng_size_t nx, ny, nz, nxp1, nyp1;
  ng_size_t nznynxp1, nznyp1nx, nznynx, nynxp1, nyp1nx, nynx;
  ng_size_t i, size_av, size_leftmost;
  ng_size_t index_u, index_v, index_msfu, index_msfv, index_msft, index_av;
  int inx, iny, inz, inxp1, inyp1;

/*
 * Variables for returning the output array with dimension names attached.
 */
  ng_size_t dsizes[1];
  NclMultiDValData att_md, return_md;
  NclVar tmp_var;
  NclStackEntry return_data;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
/*
 * Get argument # 0
 */
  u = (void*)NclGetArgValue(
           0,
           9,
           &ndims_u,
           dsizes_u,
           NULL,
           NULL,
           &type_u,
           DONT_CARE);

/*
 * Error checking on dimensions.
 */
  if(ndims_u < 3) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: u must have at least 3 dimensions");
    return(NhlFATAL);
  }
  nz   = dsizes_u[ndims_u-3];
  ny   = dsizes_u[ndims_u-2];
  nxp1 = dsizes_u[ndims_u-1];

/*
 * Get argument # 1
 */
  v = (void*)NclGetArgValue(
           1,
           9,
           &ndims_v,
           dsizes_v,
           NULL,
           NULL,
           &type_v,
           DONT_CARE);

/*
 * Error checking on dimensions.
 */
  if(ndims_v != ndims_u) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: u and v must have the same number of dimensions");
    return(NhlFATAL);
  }
  if(dsizes_v[ndims_v-3] != nz) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: The third-from-the-right dimension of v must be the same as the third-from-the-right dimension of u");
    return(NhlFATAL);
  }
/*
 * Error checking on leftmost dimension sizes.
 */
  for(i = 0; i < ndims_u-3; i++) {
    if(dsizes_u[i] != dsizes_v[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: The leftmost dimensions of u and v must be the same");
      return(NhlFATAL);
    }
  }

  nyp1 = dsizes_v[ndims_v-2];
  nx   = dsizes_v[ndims_v-1];

/*
 * Get argument # 2
 */
  msfu = (void*)NclGetArgValue(
           2,
           9,
           &ndims_msfu,
           dsizes_msfu,
           NULL,
           NULL,
           &type_msfu,
           DONT_CARE);

/*
 * Error checking on dimensions.
 */
  if(ndims_msfu < 2) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: msfu must have at least 2 dimensions");
    return(NhlFATAL);
  }
  if(ndims_msfu !=2 && ndims_msfu != (ndims_u-1)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: msfu must be 2D or have one fewer dimensions than u");
    return(NhlFATAL);
  }
  if(dsizes_msfu[ndims_msfu-2] != ny || dsizes_msfu[ndims_msfu-1] != nxp1) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: The rightmost 2 dimensions of msfu must be the same as the rightmost 2 dimensions of u");
    return(NhlFATAL);
  }

/*
 * Error checking on leftmost dimension sizes. msfu, msfv, msft, and 
 * cor can be 2D or nD.  If they are nD, they must have same leftmost
 * dimensions as other input arrays.
 */
  if(ndims_msfu > 2) {
    for(i = 0; i < ndims_u-3; i++) {
      if(dsizes_msfu[i] != dsizes_u[i]) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: If msfu is not 2-dimensional, then the leftmost dimensions of msfu and u must be the same");
        return(NhlFATAL);
      }
    }
  }

/*
 * Get argument # 3
 */
  msfv = (void*)NclGetArgValue(
           3,
           9,
           &ndims_msfv,
           dsizes_msfv,
           NULL,
           NULL,
           &type_msfv,
           DONT_CARE);

/*
 * Error checking on dimensions.
 */
  if(ndims_msfv != ndims_msfu) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: msfu, msfv, msft, and cor must have the same number of dimensions");
    return(NhlFATAL);
  }
  if(dsizes_msfv[ndims_msfv-2] != nyp1 || dsizes_msfv[ndims_msfv-1] != nx) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: The rightmost 2 dimensions of msfv must be the same as the rightmost 2 dimensions of v");
    return(NhlFATAL);
  }

/*
 * Error checking on leftmost dimension sizes.
 */
  for(i = 0; i < ndims_msfu-2; i++) {
    if(dsizes_msfv[i] != dsizes_msfu[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: The leftmost dimensions of msfv and msfu must be the same");
      return(NhlFATAL);
    }
  }

/*
 * Get argument # 4
 */
  msft = (void*)NclGetArgValue(
           4,
           9,
           &ndims_msft,
           dsizes_msft,
           NULL,
           NULL,
           &type_msft,
           DONT_CARE);

/*
 * Error checking on dimensions.
 */
  if(ndims_msft != ndims_msfu) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: msfu, msfv, msft, and cor must have the same number of dimensions");
    return(NhlFATAL);
  }
  if(dsizes_msft[ndims_msft-2] != ny || dsizes_msft[ndims_msft-1] != nx) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: The rightmost 2 dimensions of msft must be the same as the rightmost 2 dimensions of th");
    return(NhlFATAL);
  }

/*
 * Error checking on leftmost dimension sizes.
 */
  for(i = 0; i < ndims_msfu-2; i++) {
    if(dsizes_msft[i] != dsizes_msfu[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: The leftmost dimensions of msft and msfu must be the same");
      return(NhlFATAL);
    }
  }

/*
 * Get argument # 5
 */
  cor = (void*)NclGetArgValue(
           5,
           9,
           &ndims_cor,
           dsizes_cor,
           NULL,
           NULL,
           &type_cor,
           DONT_CARE);

/*
 * Error checking on dimensions.
 */
  if(ndims_cor != ndims_msft) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: msfu, msfv, msft, and cor must have the same number of dimensions");
    return(NhlFATAL);
  }

/*
 * Error checking on dimension sizes.
 */
  for(i = 0; i < ndims_msft; i++) {
    if(dsizes_cor[i] != dsizes_msft[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: The dimensions of cor and msft must be the same");
      return(NhlFATAL);
    }
  }

/*
 * Get argument # 6
 */
  dx = (void*)NclGetArgValue(
           6,
           9,
           NULL,
           NULL,
           NULL,
           NULL,
           &type_dx,
           DONT_CARE);
  tmp_dx = coerce_input_double(dx,type_dx,1,0,NULL,NULL);

/*
 * Get argument # 7
 */
  dy = (void*)NclGetArgValue(
           7,
           9,
           NULL,
           NULL,
           NULL,
           NULL,
           &type_dy,
           DONT_CARE);
  tmp_dy = coerce_input_double(dy,type_dy,1,0,NULL,NULL);

/*
 * Get argument # 8
 */
  opt = (int*)NclGetArgValue(
           8,
           9,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);

  nynx     = ny * nx;
  nznynx   = nz * nynx;
  nynxp1   = ny * nxp1;
  nyp1nx   = nyp1 * nx;
  nznynxp1 = nz * nynxp1;
  nznyp1nx = nz * nyp1nx;

/*
 * Test dimension sizes.
 */
    if((nxp1 > INT_MAX) || (nyp1 > INT_MAX) || (nz > INT_MAX) || 
       (nx > INT_MAX) ||(ny > INT_MAX)) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: one or more dimension sizes is greater than INT_MAX");
      return(NhlFATAL);
    }
    inx = (int) nx;
    iny = (int) ny;
    inz = (int) nz;
    inxp1 = (int) nxp1;
    inyp1 = (int) nyp1;

/*
 * Calculate size of leftmost dimensions, and set
 * dimension sizes for output array.
 */
  dsizes_av = (ng_size_t*)calloc(ndims_u,sizeof(ng_size_t));  
  if( dsizes_av == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: Unable to allocate memory for holding dimension sizes");
    return(NhlFATAL);
  }

  size_leftmost = 1;
  for(i = 0; i < ndims_u-3; i++) {
    size_leftmost *= dsizes_u[i];
    dsizes_av[i] = dsizes_u[i];
  }
  size_av = size_leftmost * nznynx;
  dsizes_av[ndims_u-1] = nx;
  dsizes_av[ndims_u-2] = ny;
  dsizes_av[ndims_u-3] = nz;

/*
 * Retrieve dimension names from the "u" and "v" variables, if any.
 *
 * U's dimension names will be used for the output, except for the
 * rightmost dimension which will be replaced by V's rightmost dimension
 * name.
 */
  dim_info   = get_wrf_dim_info(0,9,ndims_u,dsizes_u);
  dim_info_v = get_wrf_dim_info(1,9,ndims_v,dsizes_v);

  dim_info[ndims_u-1].dim_size = nx;
  dim_info[ndims_u-2].dim_size = ny;
  dim_info[ndims_u-3].dim_size = nz;
  dim_info[ndims_u-1].dim_quark = dim_info_v[ndims_v-1].dim_quark;

/* 
 * Allocate space for coercing input arrays.  If any of the input
 * is already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 */

/*
 * Allocate space for tmp_u.
 */
  if(type_u != NCL_double) {
    tmp_u = (double *)calloc(nznynxp1,sizeof(double));
    if(tmp_u == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: Unable to allocate memory for coercing u to double");
      return(NhlFATAL);
    }
  }

/*
 * Allocate space for tmp_v.
 */
  if(type_v != NCL_double) {
    tmp_v = (double *)calloc(nznyp1nx,sizeof(double));
    if(tmp_v == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: Unable to allocate memory for coercing v to double");
      return(NhlFATAL);
    }
  }
/*
 * Allocate space for coercing msfu, msfv, and cor to double precision.
 * These arrays can be 2D or nD, so take this into account. If one of
 * them is 2D, then all three of them have to be 2D.
 */
  if(ndims_msfu == 2) {
    tmp_msfu = coerce_input_double(msfu,type_msfu,nynxp1,0,NULL,NULL);
    tmp_msfv = coerce_input_double(msfv,type_msfv,nyp1nx,0,NULL,NULL);
    tmp_msft = coerce_input_double(msft,type_msft,nynx,0,NULL,NULL);
    tmp_cor  = coerce_input_double(cor,type_cor,nynx,0,NULL,NULL);
  }
  else {
/*
 * Allocate space for tmp_msfu.
 */
    if(type_msfu != NCL_double) {
      tmp_msfu = (double*)calloc(nynxp1,sizeof(double));
      if(tmp_msfu == NULL) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: Unable to allocate memory for coercing msfu to double");
        return(NhlFATAL);
      }
    }
/*
 * Allocate space for tmp_msfv.
 */
    if(type_msfv != NCL_double) {
      tmp_msfv = (double*)calloc(nyp1nx,sizeof(double));
      if(tmp_msfv == NULL) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: Unable to allocate memory for coercing msfv to double");
        return(NhlFATAL);
      }
    }
/*
 * Allocate space for tmp_msft.
 */
    if(type_msft != NCL_double) {
      tmp_msft = (double*)calloc(nynx,sizeof(double));
      if(tmp_msft == NULL) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: Unable to allocate memory for coercing msft to double");
        return(NhlFATAL);
      }
    }
/*
 * Allocate space for tmp_cor.
 */
    if(type_cor != NCL_double) {
      tmp_cor = (double *)calloc(nynx,sizeof(double));
      if(tmp_cor == NULL) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: Unable to allocate memory for coercing cor to double");
        return(NhlFATAL);
      }
    }
  }
/*
 * The output type defaults to float, unless any input arrays are double.
 */
  if(type_u    == NCL_double || type_v    == NCL_double || 
     type_msfu == NCL_double || type_msfv == NCL_double || 
     type_msft == NCL_double || type_cor  == NCL_double) {
    type_av     = NCL_double;
    type_obj_av = nclTypedoubleClass;
  }
  else { 
    type_av     = NCL_float;
    type_obj_av = nclTypefloatClass;
  }

/* 
 * Allocate space for output array.
 */
  if(type_av != NCL_double) {
    av     = (void *)calloc(size_av, sizeof(float));
    tmp_av = (double *)calloc(nznynx,sizeof(double));
    if(av == NULL || tmp_av == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }
  else {
    av = (void *)calloc(size_av, sizeof(double));
    if(av == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_avo: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }

/*
 * Call the Fortran routine.
 */
  index_u = index_v = index_msfu = index_msfv = index_msft = index_av = 0;
  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of u (tmp_u) to double if necessary.
 */
    if(type_u != NCL_double) {
      coerce_subset_input_double(u,tmp_u,index_u,type_u,nznynxp1,0,NULL,NULL);
    }
    else {
      tmp_u = &((double*)u)[index_u];
    }

/*
 * Coerce subsection of v (tmp_v) to double if necessary.
 */
    if(type_v != NCL_double) {
      coerce_subset_input_double(v,tmp_v,index_v,type_v,nznyp1nx,0,NULL,NULL);
    }
    else {
      tmp_v = &((double*)v)[index_v];
    }

/*
 * msfu, msfv, msft, and cor can be 2D or nD, so account
 * for that here. If they are 2D, they've already been coerced
 * before the loop.
 */
    if(ndims_msfu > 2) {
/*
 * Coerce subsection of msfu (tmp_msfu) to double if necessary.
 */
      if(type_msfu != NCL_double) {
        coerce_subset_input_double(msfu,tmp_msfu,index_msfu,type_msfu,nynxp1,0,NULL,NULL);
      }
      else {
        tmp_msfu = &((double*)msfu)[index_msfu];
      }
/*
 * Coerce subsection of msfv (tmp_msfv) to double if necessary.
 */
      if(type_msfv != NCL_double) {
        coerce_subset_input_double(msfv,tmp_msfv,index_msfv,type_msfv,nyp1nx,0,NULL,NULL);
      }
      else {
        tmp_msfv = &((double*)msfv)[index_msfv];
      }
/*
 * Coerce subsection of msft (tmp_msft) to double if necessary.
 */
      if(type_msft != NCL_double) {
        coerce_subset_input_double(msft,tmp_msft,index_msft,type_msft,nynx,0,NULL,NULL);
      }
      else {
        tmp_msft = &((double*)msft)[index_msft];
      }

/*
 * Coerce subsection of cor (tmp_cor) to double if necessary.
 */
      if(type_cor != NCL_double) {
        coerce_subset_input_double(cor,tmp_cor,index_msft,type_cor,nynx,0,NULL,NULL);
      }
      else {
        tmp_cor = &((double*)cor)[index_msft];
      }
    }

/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_av == NCL_double) tmp_av = &((double*)av)[index_av];

    NGCALLF(dcomputeabsvort,DCOMPUTEABSVORT)(tmp_av, tmp_u, tmp_v, tmp_msfu,
                                             tmp_msfv, tmp_msft, tmp_cor,
                                             tmp_dx, tmp_dy, &inx, &iny, &inz,
                                             &inxp1, &inyp1);
    if(type_av != NCL_double) {
      coerce_output_float_only(av,tmp_av,nznynx,index_av);
    }
    index_u    += nznynxp1;
    index_v    += nznyp1nx;
    index_av   += nznynx;
    if(ndims_msfu > 2) {
      index_msfu += nynxp1;
      index_msfv += nyp1nx;
      index_msft += nynx;
    }
  }

/*
 * Free unneeded memory.
 */
  if(type_u    != NCL_double) NclFree(tmp_u);
  if(type_v    != NCL_double) NclFree(tmp_v);
  if(type_msfu != NCL_double) NclFree(tmp_msfu);
  if(type_msfv != NCL_double) NclFree(tmp_msfv);
  if(type_msft != NCL_double) NclFree(tmp_msft);
  if(type_cor  != NCL_double) NclFree(tmp_cor);
  if(type_av   != NCL_double) NclFree(tmp_av);

/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)av,
                            NULL,
                            ndims_u,
                            dsizes_av,
                            TEMPORARY,
                            NULL,
                            type_obj_av
                            );

  NclFree(dsizes_av);
/*
 * Set up some attributes ("description" and "units") to return.
 */
  cdescription = (char *)calloc(19,sizeof(char));
  strcpy(cdescription,"Absolute Vorticity");
  description  = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *description = NrmStringToQuark(cdescription);
  free(cdescription);

  cunits       = (char *)calloc(9,sizeof(char));
  strcpy(cunits,"10-5 s-1");
  units        = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *units       = NrmStringToQuark(cunits);
  free(cunits);

/*
 * Set up attributes to return.
 */
  att_id = _NclAttCreate(NULL,NULL,Ncl_Att,0,NULL);

  dsizes[0] = 1;
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)description,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "description",
             att_md,
             NULL
             );
    
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)units,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "units",
             att_md,
             NULL
             );
    
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          att_id,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );

/*
 * Return output grid and attributes to NCL.
 */
  NclFree(dim_info);
  NclFree(dim_info_v);
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);

}

NhlErrorTypes wrf_helicity_W( void )
{

/*
 * Input variables
 *
 * Argument # 0
 */
  void *u;
  double *tmp_u = NULL;
  int ndims_u;
  ng_size_t dsizes_u[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_u;

/*
 * Argument # 1
 */
  void *v;
  double *tmp_v = NULL;
  int ndims_v;
  ng_size_t dsizes_v[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_v;

/*
 * Argument # 2
 */
  void *z;
  double *tmp_z = NULL;
  int ndims_z;
  ng_size_t dsizes_z[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_z;

/*
 * Argument # 3
 */
  void *ter;
  double *tmp_ter = NULL;
  int ndims_ter;
  ng_size_t dsizes_ter[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_ter;

/*
 * Argument # 4
 */
  void *top;
  double *tmp_top = NULL;
  NclBasicDataTypes type_top;

/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info;

/*
 * Return variable
 */
  void *sreh;
  double *tmp_sreh = NULL;
  int att_id;
  NclBasicDataTypes type_sreh;
  NclObjClass type_obj_sreh;
  NclQuark *description, *units;
  char *cdescription, *cunits;

/*
 * Various
 */
  ng_size_t i, miy, mjx, mkzh, mxy, mxyz;
  ng_size_t size_sreh, size_leftmost, index_u, index_ter;
  int imiy, imjx, imkzh;
/*
 * Variables for returning the output array with dimension names attached.
 */
  ng_size_t dsizes[1];
  NclMultiDValData att_md, return_md;
  NclVar tmp_var;
  NclStackEntry return_data;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
/*
 * Get argument # 0
 */
  u = (void*)NclGetArgValue(
           0,
           5,
           &ndims_u,
           dsizes_u,
           NULL,
           NULL,
           &type_u,
           DONT_CARE);

/*
 * Error checking on dimensions.
 */
  if(ndims_u < 3) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_helicity: u must have at least 3 dimensions");
    return(NhlFATAL);
  }
  mkzh = dsizes_u[ndims_u-3];
  mjx  = dsizes_u[ndims_u-2];
  miy  = dsizes_u[ndims_u-1];

/*
 * Get argument # 1
 */
  v = (void*)NclGetArgValue(
           1,
           5,
           &ndims_v,
           dsizes_v,
           NULL,
           NULL,
           &type_v,
           DONT_CARE);

/*
 * Get argument # 2
 */
  z = (void*)NclGetArgValue(
           2,
           5,
           &ndims_z,
           dsizes_z,
           NULL,
           NULL,
           &type_z,
           DONT_CARE);

/*
 * Error checking on dimension sizes.
 */
  if(ndims_u != ndims_v || ndims_u != ndims_z) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_helicity: u, v, and z must have the same dimensions");
    return(NhlFATAL);
  }

  for(i = 0; i < ndims_u; i++) {
    if(dsizes_u[i] != dsizes_v[i] || dsizes_u[i] != dsizes_z[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_helicity: u, v, and z must have the same dimensions");
      return(NhlFATAL);
    }
  }

/*
 * Get argument # 3
 */
  ter = (void*)NclGetArgValue(
           3,
           5,
           &ndims_ter,
           dsizes_ter,
           NULL,
           NULL,
           &type_ter,
           DONT_CARE);

/*
 * Error checking on dimensions.
 */
  if(ndims_ter != (ndims_u-1)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_helicity: ter must have one fewer dimension sizes than u, v, z");
    return(NhlFATAL);
  }

  if(dsizes_ter[ndims_ter-2] != mjx || dsizes_ter[ndims_ter-1] != miy) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_helicity: The rightmost two dimensions of ter must be the same as the rightmost two dimensions of u, v, z");
    return(NhlFATAL);
  }

/*
 * Error checking on leftmost dimension sizes.
 */
  for(i = 0; i < ndims_ter-2; i++) {
    if(dsizes_ter[i] != dsizes_u[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_helicity: The leftmost dimensions of ter and u, v, z must be the same");
      return(NhlFATAL);
    }
  }

/*
 * Get argument # 4
 */
  top = (void*)NclGetArgValue(
           4,
           5,
           NULL,
           NULL,
           NULL,
           NULL,
           &type_top,
           DONT_CARE);
  tmp_top = coerce_input_double(top,type_top,1,0,NULL,NULL);

  mxy  = mjx * miy;
  mxyz = mxy * mkzh;

  if((miy > INT_MAX) || (mjx > INT_MAX) || (mkzh > INT_MAX)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_helicity: one or more dimension sizes is greater than INT_MAX");
    return(NhlFATAL);
  }
  imiy = (int) miy;
  imjx = (int) mjx;
  imkzh = (int) mkzh;

/*
 * Calculate size of leftmost dimensions.
 */
  size_leftmost = 1;
  for(i = 0; i < ndims_ter-2; i++) size_leftmost *= dsizes_ter[i];

  size_sreh = size_leftmost * mxy;

/*
 * Retrieve dimension names from the "ter", if any.
 *
 * ter's dimension names will be used for the output.
 */
  dim_info   = get_wrf_dim_info(3,5,ndims_ter,dsizes_ter);

/* 
 * Allocate space for coercing input arrays.  If any of the input
 * is already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 */

/*
 * Allocate space for tmp_u.
 */
  if(type_u != NCL_double) {
    tmp_u = (double *)calloc(mxyz,sizeof(double));
    if(tmp_u == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_helicity: Unable to allocate memory for coercing u to double");
      return(NhlFATAL);
    }
  }

/*
 * Allocate space for tmp_v.
 */
  if(type_v != NCL_double) {
    tmp_v = (double *)calloc(mxyz,sizeof(double));
    if(tmp_v == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_helicity: Unable to allocate memory for coercing v to double");
      return(NhlFATAL);
    }
  }

/*
 * Allocate space for tmp_z.
 */
  if(type_z != NCL_double) {
    tmp_z = (double *)calloc(mxyz,sizeof(double));
    if(tmp_z == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_helicity: Unable to allocate memory for coercing z to double");
      return(NhlFATAL);
    }
  }

/*
 * Allocate space for tmp_ter.
 */
  if(type_ter != NCL_double) {
    tmp_ter = (double *)calloc(mxy,sizeof(double));
    if(tmp_ter == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_helicity: Unable to allocate memory for coercing ter to double");
      return(NhlFATAL);
    }
  }

/*
 * The output type defaults to float, unless any input arrays are double.
 */
  if(type_u == NCL_double || type_v   == NCL_double || 
     type_z == NCL_double || type_ter == NCL_double) {
    type_sreh     = NCL_double;
    type_obj_sreh = nclTypedoubleClass;
  }
  else { 
    type_sreh     = NCL_float;
    type_obj_sreh = nclTypefloatClass;
  }

/* 
 * Allocate space for output array.
 */
  if(type_sreh != NCL_double) {
    sreh     = (void *)calloc(size_sreh, sizeof(float));
    tmp_sreh = (double *)calloc(mxy,sizeof(double));
    if(sreh == NULL || tmp_sreh == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_helicity: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }
  else {
    sreh = (void *)calloc(size_sreh, sizeof(double));
    if(sreh == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_helicity: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }

/*
 * Call the Fortran routine.
 */
  index_u = index_ter = 0;
  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of u (tmp_u) to double if necessary.
 */
    if(type_u != NCL_double) {
      coerce_subset_input_double(u,tmp_u,index_u,type_u,mxyz,0,NULL,NULL);
    }
    else {
      tmp_u = &((double*)u)[index_u];
    }

/*
 * Coerce subsection of v (tmp_v) to double if necessary.
 */
    if(type_v != NCL_double) {
      coerce_subset_input_double(v,tmp_v,index_u,type_v,mxyz,0,NULL,NULL);
    }
    else {
      tmp_v = &((double*)v)[index_u];
    }

/*
 * Coerce subsection of z (tmp_z) to double if necessary.
 */
    if(type_z != NCL_double) {
      coerce_subset_input_double(z,tmp_z,index_u,type_z,mxyz,0,NULL,NULL);
    }
    else {
      tmp_z = &((double*)z)[index_u];
    }

/*
 * Coerce subsection of ter (tmp_ter) to double if necessary.
 */
    if(type_ter != NCL_double) {
      coerce_subset_input_double(ter,tmp_ter,index_ter,type_ter,mxy,0,NULL,NULL);
    }
    else {
      tmp_ter = &((double*)ter)[index_ter];
    }


/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_sreh == NCL_double) tmp_sreh = &((double*)sreh)[index_ter];

    NGCALLF(dcalrelhl,DCALRELHL)(tmp_u, tmp_v, tmp_z, tmp_ter, tmp_top,
                                 tmp_sreh, &imiy, &imjx, &imkzh);
    if(type_sreh != NCL_double) {
      coerce_output_float_only(sreh,tmp_sreh,mxy,index_ter);
    }
    index_u   += mxyz;
    index_ter += mxy;
  }

/*
 * Free unneeded memory.
 */
  if(type_u    != NCL_double) NclFree(tmp_u);
  if(type_v    != NCL_double) NclFree(tmp_v);
  if(type_z    != NCL_double) NclFree(tmp_z);
  if(type_ter  != NCL_double) NclFree(tmp_ter);
  if(type_sreh != NCL_double) NclFree(tmp_sreh);
  if(type_top != NCL_double) NclFree(tmp_top);

/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)sreh,
                            NULL,
                            ndims_ter,
                            dsizes_ter,
                            TEMPORARY,
                            NULL,
                            type_obj_sreh
                            );

/*
 * Set up some attributes ("description" and "units") to return.
 */
  cdescription = (char *)calloc(24,sizeof(char));
  strcpy(cdescription,"Storm Relative Helicity");
  description  = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *description = NrmStringToQuark(cdescription);
  free(cdescription);

  cunits       = (char *)calloc(8,sizeof(char));
  strcpy(cunits,"m-2/s-2");
  units        = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *units       = NrmStringToQuark(cunits);
  free(cunits);

/*
 * Set up attributes to return.
 */
  att_id = _NclAttCreate(NULL,NULL,Ncl_Att,0,NULL);

  dsizes[0] = 1;
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)description,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "description",
             att_md,
             NULL
             );
    
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)units,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "units",
             att_md,
             NULL
             );
    
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          att_id,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );
  NclFree(dim_info);

/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);

}


NhlErrorTypes wrf_updraft_helicity_W( void )
{

/*
 * Input variables
 */
/*
 * Argument # 0
 */
  void *zp;
  double *tmp_zp = NULL;
  int ndims_zp;
  ng_size_t dsizes_zp[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_zp;

/*
 * Argument # 1
 */
  void *mapfct;
  double *tmp_mapfct = NULL;
  int ndims_mapfct;
  ng_size_t dsizes_mapfct[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_mapfct;

/*
 * Argument # 2
 */
  void *dx;
  double *tmp_dx = NULL;
  NclBasicDataTypes type_dx;

/*
 * Argument # 3
 */
  void *dy;
  double *tmp_dy = NULL;
  NclBasicDataTypes type_dy;

/*
 * Argument # 4
 */
  void *us;
  double *tmp_us = NULL;
  int ndims_us;
  ng_size_t dsizes_us[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_us;

/*
 * Argument # 5
 */
  void *vs;
  double *tmp_vs = NULL;
  int ndims_vs;
  ng_size_t dsizes_vs[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_vs;

/*
 * Argument # 6
 */
  void *w;
  double *tmp_w = NULL;
  int ndims_w;
  ng_size_t dsizes_w[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_w;

/*
 * Argument # 7
 */
  logical *opt;

/*
 * Possible attributes.
 */
  void *uhmnhgt = NULL;
  void *uhmxhgt = NULL;
  double *tmp_uhmnhgt, *tmp_uhmxhgt;
  logical set_uhmnhgt, set_uhmxhgt;
  NclBasicDataTypes type_uhmnhgt = NCL_none;
  NclBasicDataTypes type_uhmxhgt = NCL_none;

/*
 * Variables for retrieving attributes from "opt".
 */
  NclAttList  *attr_list;
  NclAtt  attr_obj;
  NclStackEntry stack_entry;

/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info = NULL;
  NclDimRec *dim_info_us;

/*
 * Return variable
 */
  void *uh;
  double *tmp_uh = NULL;
  int ndims_uh;
  ng_size_t *dsizes_uh;
  NclBasicDataTypes type_uh;
  NclObjClass type_obj_uh;
  int att_id;
  NclQuark *description, *units;
  char *cdescription, *cunits;

/*
 * Various
 */
  ng_size_t nx, ny, nzp1, nzp1nynx, nynx, nz, nznynx;
  ng_size_t index_zp, index_uh, index_us;
  double *tem1, *tem2;
  ng_size_t i, ndims_leftmost, size_leftmost, size_output;
  int inx, iny, inz, inzp1;

/*
 * Variables for returning the output array with dimension names attached.
 */
  ng_size_t dsizes[1];
  NclMultiDValData att_md, return_md;
  NclVar tmp_var;
  NclStackEntry return_data;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
/*
 * Get argument # 0
 */
  zp = (void*)NclGetArgValue(
           0,
           8,
           &ndims_zp,
           dsizes_zp,
           NULL,
           NULL,
           &type_zp,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_zp < 3) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: The zp array must have at least 3 dimensions");
    return(NhlFATAL);
  }
  nzp1 = dsizes_zp[ndims_zp-3];
  ny   = dsizes_zp[ndims_zp-2];
  nx   = dsizes_zp[ndims_zp-1];
  nynx = ny * nx;
  nzp1nynx = nynx * nzp1;

/*
 * Get argument # 1
 */
  mapfct = (void*)NclGetArgValue(
           1,
           8,
           &ndims_mapfct,
           dsizes_mapfct,
           NULL,
           NULL,
           &type_mapfct,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_mapfct != 2 && ndims_mapfct != (ndims_zp-1)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: The mapfct array must either be two-dimensional or have one fewer dimensions than zp");
    return(NhlFATAL);
  }
  if(dsizes_mapfct[ndims_mapfct-2] != ny || 
     dsizes_mapfct[ndims_mapfct-1] != nx) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: The rightmost dimensions of mapfct must be ny x nx");
    return(NhlFATAL);
  }

/*
 * Get argument # 2
 */
  us = (void*)NclGetArgValue(
           2,
           8,
           &ndims_us,
           dsizes_us,
           NULL,
           NULL,
           &type_us,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_us != ndims_zp) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: The us array must have the same number of dimensions as zp");
    return(NhlFATAL);
  }
  if(dsizes_us[ndims_us-2] != ny || dsizes_us[ndims_us-1] != nx) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: The rightmost dimensions of us must be nz x ny x ny");
    return(NhlFATAL);
  }
  nz = dsizes_us[ndims_us-3];
  nznynx = nz * nynx;

/*
 * Test dimension sizes.
 */
  if((nx > INT_MAX) || (ny > INT_MAX) || (nzp1 > INT_MAX) || (nz > INT_MAX)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: one or more dimension sizes is greater than INT_MAX");
    return(NhlFATAL);
  }
  inx = (int) nx;
  iny = (int) ny;
  inz = (int) nz;
  inzp1 = (int) nzp1;

/*
 * Get argument # 3
 */
  vs = (void*)NclGetArgValue(
           3,
           8,
           &ndims_vs,
           dsizes_vs,
           NULL,
           NULL,
           &type_vs,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_vs != ndims_us) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: The vs array must have the same number of dimensions as us");
    return(NhlFATAL);
  }
  if(dsizes_vs[ndims_vs-3] != nz || dsizes_vs[ndims_vs-2] != ny || 
     dsizes_vs[ndims_vs-1] != nx) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: The rightmost dimensions of vs must be nz x ny x ny");
    return(NhlFATAL);
  }
/*
 * Get argument # 4
 */
  w = (void*)NclGetArgValue(
           4,
           8,
           &ndims_w,
           dsizes_w,
           NULL,
           NULL,
           &type_w,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_w != ndims_zp || (dsizes_w[ndims_w-3] != nzp1 || 
                             dsizes_w[ndims_w-2] != ny ||
                             dsizes_w[ndims_w-1] != nx)) {
  NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: The w array must be the same dimensionality as zp");
  return(NhlFATAL);
 }

/*
 * Get argument # 5
 */
  dx = (void*)NclGetArgValue(
           5,
           8,
           NULL,
           NULL,
           NULL,
           NULL,
           &type_dx,
           DONT_CARE);
/*
 * Get argument # 6
 */
  dy = (void*)NclGetArgValue(
           6,
           8,
           NULL,
           NULL,
           NULL,
           NULL,
           &type_dy,
           DONT_CARE);

/*
 * Get argument # 7
 */
  opt = (logical*)NclGetArgValue(
           7,
           8,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);

/*
 * Start checking for attributes attached to "opt"
 */
  set_uhmnhgt = set_uhmxhgt = False;

  stack_entry = _NclGetArg(7, 8, DONT_CARE);
  switch (stack_entry.kind) {
  case NclStk_VAR:
    if (stack_entry.u.data_var->var.att_id != -1) {
      attr_obj = (NclAtt) _NclGetObj(stack_entry.u.data_var->var.att_id);
      if (attr_obj == NULL) {
        break;
      }
    }
    else {
/*
 * att_id == -1 ==> no optional args given.
 */
      break;
    }
/* 
 * Get optional arguments.
 */
    if (attr_obj->att.n_atts > 0) {
/*
 * Get list of attributes.
 */
      attr_list = attr_obj->att.att_list;
/*
 * Loop through attributes and check them. We are looking for:
 *
 *  uhmnhgt or uhmxhgt
 */
      while (attr_list != NULL) {
        if(!strcasecmp(attr_list->attname, "uhmnhgt")) {
          uhmnhgt      = attr_list->attvalue->multidval.val;
          type_uhmnhgt = attr_list->attvalue->multidval.data_type;
          set_uhmnhgt  = True;
        }
        else if(!strcasecmp(attr_list->attname, "uhmxhgt")) {
          uhmxhgt      = attr_list->attvalue->multidval.val;
          type_uhmxhgt = attr_list->attvalue->multidval.data_type;
          set_uhmxhgt  = True;
        }
        attr_list = attr_list->next;
      }
    default:
      break;
    }
  }
  if(set_uhmnhgt) {
    tmp_uhmnhgt = coerce_input_double(uhmnhgt,type_uhmnhgt,1,0,NULL,NULL);
  }
  else {
    type_uhmnhgt = NCL_double;
    tmp_uhmnhgt  = (double *)calloc(1,sizeof(double));
    *tmp_uhmnhgt = 2000.;
  }
  if(set_uhmxhgt) {
    tmp_uhmxhgt = coerce_input_double(uhmxhgt,type_uhmxhgt,1,0,NULL,NULL);
  }
  else {
    type_uhmxhgt = NCL_double;
    tmp_uhmxhgt  = (double *)calloc(1,sizeof(double));
    *tmp_uhmxhgt = 5000.;
  }

/*
 * Calculate size of leftmost dimensions.
 */
  size_leftmost  = 1;
  ndims_leftmost = ndims_zp-3;
  for(i = 0; i < ndims_leftmost; i++) {
    if(dsizes_us[i] != dsizes_zp[i] || dsizes_vs[i] != dsizes_zp[i] ||
       dsizes_w[i]  != dsizes_zp[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: The leftmost dimensions of zp, us, vs, and w must be the same");
      return(NhlFATAL);
    }
    if(ndims_mapfct > 2 && dsizes_mapfct[i] != dsizes_zp[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: mapfct must either be two-dimensional or have the same leftmost dimensions as zp");
      return(NhlFATAL);
    }
    size_leftmost *= dsizes_zp[i];
  }


/* 
 * Allocate space for coercing input arrays.  If any of the input
 * is already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 *
 * Allocate space for tmp_zp.
 */
  if(type_zp != NCL_double) {
    type_uh     = NCL_float;
    type_obj_uh = nclTypefloatClass;

    tmp_zp = (double *)calloc(nzp1nynx,sizeof(double));
    if(tmp_zp == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_uh     = NCL_double;
    type_obj_uh = nclTypedoubleClass;
  }
/*
 * Allocate space for tmp_mapfct. This array can be 2D or nD. If it
 * is 2D, go ahead and coerce the values now. Otherwise, create a temp
 * 2D array and we'll coerce the values in the loop below along with
 * everybody else.
 */
  if(ndims_mapfct == 2) {
    tmp_mapfct = coerce_input_double(mapfct,type_mapfct,nynx,0,NULL,NULL);
  }
  else {
    if(type_mapfct != NCL_double) {
      tmp_mapfct = (double *)calloc(nynx,sizeof(double));
      if(tmp_mapfct == NULL) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: Unable to allocate memory for coercing input array to double");
        return(NhlFATAL);
      }
    }
  }
/*
 * Allocate space for tmp_dx.
 */
  tmp_dx = coerce_input_double(dx,type_dx,1,0,NULL,NULL);
  if(tmp_dx == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: Unable to allocate memory for coercing input array to double");
    return(NhlFATAL);
  }
/*
 * Allocate space for tmp_dy.
 */
  tmp_dy = coerce_input_double(dy,type_dy,1,0,NULL,NULL);
  if(tmp_dy == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: Unable to allocate memory for coercing input array to double");
    return(NhlFATAL);
  }

/*
 * Allocate space for tmp_us.
 */
  if(type_us != NCL_double) {
    tmp_us = (double *)calloc(nznynx,sizeof(double));
    if(tmp_us == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
/*
 * Allocate space for tmp_vs.
 */
  if(type_vs != NCL_double) {
    tmp_vs = (double *)calloc(nznynx,sizeof(double));
    if(tmp_vs == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
/*
 * Allocate space for tmp_w.
 */
  if(type_w != NCL_double) {
    tmp_w = (double *)calloc(nzp1nynx,sizeof(double));
    if(tmp_w == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }

/*
 * Calculate size of output array.
 */
  size_output = size_leftmost * nynx;

/* 
 * Allocate space for output array.
 */
  if(type_uh != NCL_double) {
    uh = (void *)calloc(size_output, sizeof(float));
    tmp_uh = (double *)calloc(nynx,sizeof(double));
    if(tmp_uh == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: Unable to allocate memory for temporary output array");
      return(NhlFATAL);
    }
  }
  else {
    uh = (void *)calloc(size_output, sizeof(double));
  }
  if(uh == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: Unable to allocate memory for output array");
    return(NhlFATAL);
  }

/*
 * Allocate space for work arrays.
 */
  tem1 = (void *)calloc(nznynx, sizeof(double));
  tem2 = (void *)calloc(nznynx, sizeof(double));
  if( tem1 == NULL || tem2 == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: Unable to allocate memory for work arrays");
    return(NhlFATAL);
  }

/* 
 * Allocate space for output dimension sizes and set them.
 */
  ndims_uh = ndims_leftmost + 2;
  dsizes_uh = (ng_size_t*)calloc(ndims_uh,sizeof(ng_size_t));  
  if( dsizes_uh == NULL ) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: Unable to allocate memory for holding dimension sizes");
    return(NhlFATAL);
  }
  for(i = 0; i < ndims_uh-2; i++) dsizes_uh[i] = dsizes_zp[i];
  dsizes_uh[ndims_uh-2] = ny;
  dsizes_uh[ndims_uh-1] = nx;

/*
 * Retrieve dimension names from "u", if any.
 *
 * u's dimension names will be used for the output.
 */
  dim_info_us = get_wrf_dim_info(2,8,ndims_us,dsizes_us);

  if(dim_info_us != NULL) {
    dim_info = malloc(sizeof(NclDimRec)*ndims_uh);
    if(dim_info == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_updraft_helicity: Unable to allocate memory for holding dimension information");
      return(NhlFATAL);
    }
    for(i = 0; i < ndims_uh-2; i++) dim_info[i] = dim_info_us[i];
    dim_info[ndims_uh-1] = dim_info_us[ndims_us-1];
    dim_info[ndims_uh-2] = dim_info_us[ndims_us-2];
  }

/*
 * Loop across leftmost dimensions and call the Fortran routine for each
 * subsection of the input arrays.
 */
  index_zp = index_uh = index_us = 0;

  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of zp (tmp_zp) to double if necessary.
 */
    if(type_zp != NCL_double) {
      coerce_subset_input_double(zp,tmp_zp,index_zp,type_zp,nzp1nynx,
                                 0,NULL,NULL);
    }
    else {
      tmp_zp = &((double*)zp)[index_zp];
    }

/*
 * Coerce subsection of mapfct (tmp_mapfct) to double if necessary.
 */
    if(ndims_mapfct > 2) {
      if(type_mapfct != NCL_double) {
        coerce_subset_input_double(mapfct,tmp_mapfct,index_uh,
                                   type_mapfct,nynx,0,NULL,NULL);
      }
      else {
        tmp_mapfct = &((double*)mapfct)[index_uh];
      }
    }

/*
 * Coerce subsection of us (tmp_us) to double if necessary.
 */
    if(type_us != NCL_double) {
      coerce_subset_input_double(us,tmp_us,index_us,type_us,nznynx,
                                 0,NULL,NULL);
    }
    else {
      tmp_us = &((double*)us)[index_us];
    }

/*
 * Coerce subsection of vs (tmp_vs) to double if necessary.
 */
    if(type_vs != NCL_double) {
      coerce_subset_input_double(vs,tmp_vs,index_us,type_vs,nznynx,
                                 0,NULL,NULL);
    }
    else {
      tmp_vs = &((double*)vs)[index_us];
    }

/*
 * Coerce subsection of w (tmp_w) to double if necessary.
 */
    if(type_w != NCL_double) {
      coerce_subset_input_double(w,tmp_w,index_zp,type_w,nzp1nynx,0,NULL,NULL);
    }
    else {
      tmp_w = &((double*)w)[index_zp];
    }


/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_uh == NCL_double) tmp_uh = &((double*)uh)[index_uh];


/*
 * Call the Fortran routine.
 */
    NGCALLF(dcalcuh,DCALCUH)(&inx, &iny, &inz, &inzp1, tmp_zp, tmp_mapfct, 
                             tmp_dx, tmp_dy, tmp_uhmnhgt, tmp_uhmxhgt, 
                             tmp_us, tmp_vs, tmp_w, tmp_uh, tem1, tem2);
/*
 * Coerce output back to float if necessary.
 */
    if(type_uh == NCL_float) {
      coerce_output_float_only(uh,tmp_uh,nynx,index_uh);
    }
    index_zp += nzp1nynx;
    index_uh += nynx;
    index_us += nznynx;
  }

/*
 * Free unneeded memory.
 */
  if(type_zp != NCL_double)      NclFree(tmp_zp);
  if(type_mapfct != NCL_double)  NclFree(tmp_mapfct);
  if(type_dx != NCL_double)      NclFree(tmp_dx);
  if(type_dy != NCL_double)      NclFree(tmp_dy);
  if(type_us != NCL_double)      NclFree(tmp_us);
  if(type_vs != NCL_double)      NclFree(tmp_vs);
  if(type_w != NCL_double)       NclFree(tmp_w);
  if(type_uh != NCL_double)      NclFree(tmp_uh);
  if(type_uhmnhgt != NCL_double || ! set_uhmnhgt) NclFree(tmp_uhmnhgt);
  if(type_uhmxhgt != NCL_double || ! set_uhmxhgt) NclFree(tmp_uhmxhgt);
  NclFree(tem1);
  NclFree(tem2);

/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)uh,
                            NULL,
                            ndims_uh,
                            dsizes_uh,
                            TEMPORARY,
                            NULL,
                            type_obj_uh
                            );

  NclFree(dsizes_uh);
/*
 * Set up some attributes ("description" and "units") to return.
 */
  cdescription = (char *)calloc(17,sizeof(char));
  strcpy(cdescription,"Updraft Helicity");
  description  = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *description = NrmStringToQuark(cdescription);
  free(cdescription);

  cunits       = (char *)calloc(8,sizeof(char));
  strcpy(cunits,"m-2/s-2");
  units        = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *units       = NrmStringToQuark(cunits);
  free(cunits);

/*
 * Set up attributes to return.
 */
  att_id = _NclAttCreate(NULL,NULL,Ncl_Att,0,NULL);

  dsizes[0] = 1;
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)description,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "description",
             att_md,
             NULL
             );
    
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)units,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "units",
             att_md,
             NULL
             );
    
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          att_id,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );
  NclFree(dim_info);
  NclFree(dim_info_us);

/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);
}

NhlErrorTypes wrf_ll_to_ij_W( void )
{

/*
 * Input variables
 */
/*
 * Argument # 0
 */
  void *lon;
  double *tmp_lon = NULL;
  int ndims_lon;
  ng_size_t dsizes_lon[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_lon;
/*
 * Argument # 1
 */
  void *lat;
  double *tmp_lat = NULL;
  int ndims_lat;
  ng_size_t dsizes_lat[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_lat;

/*
 * Argument # 2
 */
  logical *opt;

/*
 * Variables for retrieving attributes from "opt".
 */
  NclAttList  *attr_list;
  NclAtt  attr_obj;
  NclStackEntry stack_entry;

/*
 * Variables that can be set via attributes.
 */
  int map_proj;
  void *truelat1 = NULL;
  void *truelat2 = NULL;
  void *stand_lon = NULL;
  void *ref_lat = NULL;
  void *ref_lon = NULL;
  void *pole_lat = NULL;
  void *pole_lon = NULL;
  void *knowni = NULL;
  void *knownj = NULL;
  void *dx = NULL;
  void *dy = NULL;
  void *latinc = NULL;
  void *loninc = NULL;

  double *tmp_truelat1, *tmp_truelat2, *tmp_stand_lon;
  double *tmp_ref_lat, *tmp_ref_lon, *tmp_pole_lat, *tmp_pole_lon;
  double *tmp_knowni, *tmp_knownj, *tmp_dx, *tmp_dy, *tmp_latinc, *tmp_loninc;

  NclBasicDataTypes type_truelat1 = NCL_none;
  NclBasicDataTypes type_truelat2 = NCL_none;
  NclBasicDataTypes type_stand_lon = NCL_none;
  NclBasicDataTypes type_ref_lat = NCL_none;
  NclBasicDataTypes type_ref_lon = NCL_none;
  NclBasicDataTypes type_pole_lat = NCL_none;
  NclBasicDataTypes type_pole_lon = NCL_none;
  NclBasicDataTypes type_knowni = NCL_none;
  NclBasicDataTypes type_knownj = NCL_none;
  NclBasicDataTypes type_dx = NCL_none;
  NclBasicDataTypes type_dy = NCL_none;
  NclBasicDataTypes type_latinc = NCL_none;
  NclBasicDataTypes type_loninc = NCL_none;

  logical set_map_proj, set_truelat1, set_truelat2, set_stand_lon, set_ref_lat;
  logical set_ref_lon, set_pole_lat, set_pole_lon,  set_knowni, set_knownj;
  logical set_dx, set_dy, set_latinc, set_loninc;

/*
 * Return variable
 */
  void *loc;
  double *tmp_loc;
  int ndims_loc;
  ng_size_t *dsizes_loc;
  NclBasicDataTypes type_loc;
  NclObjClass type_obj_loc;
/*
 * Variables for returning the output array with attributes attached.
 */
  NclMultiDValData return_md;
  NclVar tmp_var;
  NclStackEntry return_data;
/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info;

/*
 * Various
 */
  int npts, i;
  int errstat;
  char* errmsg;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
/*
 * Get argument # 0
 */
  lon = (void*)NclGetArgValue(
           0,
           3,
           &ndims_lon,
           dsizes_lon,
           NULL,
           NULL,
           &type_lon,
           DONT_CARE);

/*
 * Get argument # 1
 */
  lat = (void*)NclGetArgValue(
           1,
           3,
           &ndims_lat,
           dsizes_lat,
           NULL,
           NULL,
           &type_lat,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_lon != ndims_lat) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ll_to_ij: lat and lon must have the same number of dimensions");
    return(NhlFATAL);
  }

  for(i = 0; i < ndims_lat; i++) {
    if(dsizes_lon[i] != dsizes_lat[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ll_to_ij: lat and lon must have the same dimension sizes");
      return(NhlFATAL);
    }
  }

/*
 * Get argument # 2
 */
  opt = (logical*)NclGetArgValue(
           2,
           3,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);

/*
 * Calculate size of lat/lon dimensions.
 */
  npts = 1;
  for(i = 0; i < ndims_lat; i++) npts *= dsizes_lat[i];

/*
 * Start checking for attributes attached to "opt". Some are optional,
 * and some are not.  We'll check them later.
 */

  set_map_proj = set_truelat1 = set_truelat2 = set_stand_lon = False;
  set_ref_lat = set_ref_lon = set_pole_lat = set_pole_lon = False;
  set_knowni = set_knownj = set_dx = set_dy = set_latinc = set_loninc = False;

  stack_entry = _NclGetArg(2, 3, DONT_CARE);
  switch (stack_entry.kind) {
  case NclStk_VAR:
    if (stack_entry.u.data_var->var.att_id != -1) {
      attr_obj = (NclAtt) _NclGetObj(stack_entry.u.data_var->var.att_id);
      if (attr_obj == NULL) {
        break;
      }
    }
    else {
/*
 * att_id == -1 ==> no optional args given.
 */
      break;
    }
/* 
 * Get optional arguments.
 */
    if (attr_obj->att.n_atts > 0) {
/*
 * Get list of attributes.
 */
      attr_list = attr_obj->att.att_list;
/*
 * Loop through attributes and check them. We are looking for:
 *
 *   map_proj
 *   truelat1, truelat2
 *   stand_lon
 *   ref_lat, ref_lon
 *   pole_lat, pole_lon
 *   knowni, knownj
 *   dx, dy
 *   latinc, loninc
 */
      while (attr_list != NULL) {
        if(!strcasecmp(attr_list->attname, "map_proj")) {
          map_proj = *(int *)attr_list->attvalue->multidval.val;
          set_map_proj = True;
        }
        else if(!strcasecmp(attr_list->attname, "truelat1")) {
          truelat1      = attr_list->attvalue->multidval.val;
          type_truelat1 = attr_list->attvalue->multidval.data_type;
          set_truelat1  = True;
        }
        else if(!strcasecmp(attr_list->attname, "truelat2")) {
          truelat2      = attr_list->attvalue->multidval.val;
          type_truelat2 = attr_list->attvalue->multidval.data_type;
          set_truelat2  = True;
        }
        else if(!strcasecmp(attr_list->attname, "stand_lon")) {
          stand_lon      = attr_list->attvalue->multidval.val;
          type_stand_lon = attr_list->attvalue->multidval.data_type;
          set_stand_lon  = True;
        }
        else if(!strcasecmp(attr_list->attname, "ref_lat")) {
          ref_lat      = attr_list->attvalue->multidval.val;
          type_ref_lat = attr_list->attvalue->multidval.data_type;
          set_ref_lat  = True;
        }
        else if(!strcasecmp(attr_list->attname, "ref_lon")) {
          ref_lon      = attr_list->attvalue->multidval.val;
          type_ref_lon = attr_list->attvalue->multidval.data_type;
          set_ref_lon  = True;
        }
        else if(!strcasecmp(attr_list->attname, "pole_lat")) {
          pole_lat      = attr_list->attvalue->multidval.val;
          type_pole_lat = attr_list->attvalue->multidval.data_type;
          set_pole_lat  = True;
        }
        else if(!strcasecmp(attr_list->attname, "pole_lon")) {
          pole_lon      = attr_list->attvalue->multidval.val;
          type_pole_lon = attr_list->attvalue->multidval.data_type;
          set_pole_lon  = True;
        }
        else if(!strcasecmp(attr_list->attname, "knowni")) {
          knowni      = attr_list->attvalue->multidval.val;
          type_knowni = attr_list->attvalue->multidval.data_type;
          set_knowni  = True;
        }
        else if(!strcasecmp(attr_list->attname, "knownj")) {
          knownj      = attr_list->attvalue->multidval.val;
          type_knownj = attr_list->attvalue->multidval.data_type;
          set_knownj  = True;
        }
        else if(!strcasecmp(attr_list->attname, "dx")) {
          dx      = attr_list->attvalue->multidval.val;
          type_dx = attr_list->attvalue->multidval.data_type;
          set_dx  = True;
        }
        else if(!strcasecmp(attr_list->attname, "dy")) {
          dy      = attr_list->attvalue->multidval.val;
          type_dy = attr_list->attvalue->multidval.data_type;
          set_dy  = True;
        }
        else if(!strcasecmp(attr_list->attname, "latinc")) {
          latinc      = attr_list->attvalue->multidval.val;
          type_latinc = attr_list->attvalue->multidval.data_type;
          set_latinc  = True;
        }
        else if(!strcasecmp(attr_list->attname, "loninc")) {
          loninc      = attr_list->attvalue->multidval.val;
          type_loninc = attr_list->attvalue->multidval.data_type;
          set_loninc  = True;
        }
        attr_list = attr_list->next;
      }
    default:
      break;
    }
  }

/*
 * Check for attributes that need to be set, or set to a certain value.
 *
 * Check MAP_PROJ. Must be set.
 */
  if(!set_map_proj) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ll_to_ij: The MAP_PROJ attribute must be set");
    return(NhlFATAL);
  }
  else if(map_proj != 1 && map_proj != 2 && map_proj != 3 && map_proj != 6) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ll_to_ij: The MAP_PROJ attribute must be set to 1, 2, 3, or 6");
    return(NhlFATAL);
  }

/*
 * Check TRUELAT1. Must be set in some cases.
 */
  if( (map_proj == 1 || map_proj == 2 || map_proj == 3) && !set_truelat1) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ll_to_ij: The TRUELAT1 attribute must be set if MAP_PROJ is 1, 2, or 3");
    return(NhlFATAL);
  }
  if(set_truelat1) {
    tmp_truelat1 = coerce_input_double(truelat1,type_truelat1,1,0,NULL,NULL);
  }
  else {
    tmp_truelat1  = (double *)calloc(1,sizeof(double));
    *tmp_truelat1 = 0.;
  }

/*
 * Check TRUELAT2. Must be set in some cases.
 */
  if( map_proj == 1 && !set_truelat2) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ll_to_ij: The TRUELAT2 attribute must be set if MAP_PROJ is 1");
    return(NhlFATAL);
  }
  if(set_truelat2) {
    tmp_truelat2 = coerce_input_double(truelat2,type_truelat2,1,0,NULL,NULL);
  }
  else {
    tmp_truelat2  = (double *)calloc(1,sizeof(double));
    *tmp_truelat2 = 0.;
  }

/*
 * Check STAND_LON. Must be set.
 */
  if(!set_stand_lon) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ll_to_ij: The STAND_LON attribute must be set");
    return(NhlFATAL);
  }
  else {
    tmp_stand_lon = coerce_input_double(stand_lon,type_stand_lon,1,0,NULL,NULL);
  }

/*
 * Check REF_LAT/REF_LON. Must be set.
 */
  if(!set_ref_lat || !set_ref_lon) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ll_to_ij: The REF_LAT/REF_LON attributes must be set");
    return(NhlFATAL);
  }
  else {
    tmp_ref_lat = coerce_input_double(ref_lat,type_ref_lat,1,0,NULL,NULL);
    tmp_ref_lon = coerce_input_double(ref_lon,type_ref_lon,1,0,NULL,NULL);
  }

/*
 * Check POLE_LAT/POLE_LON.
 */
  if(set_pole_lat) {
    tmp_pole_lat = coerce_input_double(pole_lat,type_pole_lat,1,0,NULL,NULL);
  }
  else {
    tmp_pole_lat  = (double *)calloc(1,sizeof(double));
    *tmp_pole_lat = 90.;
  }

  if(set_pole_lon) {
    tmp_pole_lon = coerce_input_double(pole_lon,type_pole_lon,1,0,NULL,NULL);
  }
  else {
    tmp_pole_lon  = (double *)calloc(1,sizeof(double));
    *tmp_pole_lon = 0.;
  }

/*
 * Check KNOWNI/KNOWNJ. Must be set.
 */
  if(!set_knowni || !set_knownj) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ll_to_ij: The KNOWNI/KNOWNJ attributes must be set");
    return(NhlFATAL);
  }
  else {
    tmp_knowni = coerce_input_double(knowni,type_knowni,1,0,NULL,NULL);
    tmp_knownj = coerce_input_double(knownj,type_knownj,1,0,NULL,NULL);
  }

/*
 * Check DX/DY. Must be set in some cases.
 */
  if( (map_proj == 1 || map_proj == 2 || map_proj == 3) &&
      (!set_dx || !set_dy)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ll_to_ij: The DX/DY attributes must be set if MAP_PROJ is 1, 2, or 3");
    return(NhlFATAL);
  }
  if(set_dx) {
    tmp_dx = coerce_input_double(dx,type_dx,1,0,NULL,NULL);
  }
  else {
    tmp_dx  = (double *)calloc(1,sizeof(double));
    *tmp_dx = 0.;
  }

  if(set_dy) {
    tmp_dy = coerce_input_double(dy,type_dy,1,0,NULL,NULL);
  }
  else {
    tmp_dy  = (double *)calloc(1,sizeof(double));
    *tmp_dy = 0.;
  }

/*
 * Check LATINC/LONINC. Must be set in some cases.
 */
  if( map_proj == 6 && (!set_latinc || !set_loninc)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ll_to_ij: The LATINC/LONINC attributes must be set if MAP_PROJ is 6");
    return(NhlFATAL);
  }
  if(set_latinc) {
    tmp_latinc = coerce_input_double(latinc,type_latinc,1,0,NULL,NULL);
  }
  else {
    tmp_latinc  = (double *)calloc(1,sizeof(double));
    *tmp_latinc = 0.;
  }

  if(set_loninc) {
    tmp_loninc = coerce_input_double(loninc,type_loninc,1,0,NULL,NULL);
  }
  else {
    tmp_loninc  = (double *)calloc(1,sizeof(double));
    *tmp_loninc = 0.;
  }

/*
 * The output type defaults to float, unless either of the lat/lon arrays
 * are double.
 */
  type_loc     = NCL_float;
  type_obj_loc = nclTypefloatClass;

/* 
 * Allocate space for coercing input arrays.  If any of the input
 * is already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 */
/*
 * Allocate space for tmp_lat.
 */
  if(type_lat != NCL_double) {
    tmp_lat = (double *)calloc(npts,sizeof(double));
    if(tmp_lat == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ll_to_ij: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_loc     = NCL_double;
    type_obj_loc = nclTypedoubleClass;
  }

/*
 * Allocate space for tmp_lon.
 */
  if(type_lon != NCL_double) {
    tmp_lon = (double *)calloc(npts,sizeof(double));
    if(tmp_lon == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ll_to_ij: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_loc     = NCL_double;
    type_obj_loc = nclTypedoubleClass;
  }

/* 
 * Allocate space for output array.
 */
  tmp_loc = (double *)calloc(2,sizeof(double));
  if(type_loc != NCL_double) loc = (void *)calloc(2*npts, sizeof(float));
  else                       loc = (void *)calloc(2*npts, sizeof(double));
  if(loc == NULL || tmp_loc == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ll_to_ij: Unable to allocate memory for output array");
    return(NhlFATAL);
  }

/* 
 * Allocate space for output dimension sizes and set them.
 */
  if(is_scalar(ndims_lat,dsizes_lat)) {
    ndims_loc = 1;
  }
  else {
    ndims_loc = ndims_lat + 1;
  }
  dsizes_loc = (ng_size_t*)calloc(ndims_loc,sizeof(ng_size_t));  
  if( dsizes_loc == NULL ) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ll_to_ij: Unable to allocate memory for holding dimension sizes");
    return(NhlFATAL);
  }
  for(i = 0; i < ndims_loc-1; i++) dsizes_loc[i+1] = dsizes_lat[i];
  dsizes_loc[0] = 2;

  /* Allocate space for errmsg*/
  errmsg = (char *) calloc(ERRLEN, sizeof(char))

/*
 * Loop across all lat/lon points and call the Fortran routine for each
 * point.
 */
  for(i = 0; i < npts; i++) {
/*
 * Coerce subsection of lat (tmp_lat) to double if necessary.
 */
    if(type_lat != NCL_double) {
      coerce_subset_input_double(lat,tmp_lat,i,type_lat,1,0,NULL,NULL);
    }
    else {
      tmp_lat = &((double*)lat)[i];
    }

/*
 * Coerce subsection of lon (tmp_lon) to double if necessary.
 */
    if(type_lon != NCL_double) {
      coerce_subset_input_double(lon,tmp_lon,i,type_lon,1,0,NULL,NULL);
    }
    else {
      tmp_lon = &((double*)lon)[i];
    }

/*
 * Call the Fortran routine.
 */
    errstat=0;
    errmsg="";
    NGCALLF(dlltoij,DLLTOIJ)(&map_proj, tmp_truelat1, tmp_truelat2, 
                             tmp_stand_lon, tmp_ref_lat, tmp_ref_lon, 
                             tmp_pole_lat, tmp_pole_lon, tmp_knowni,
                             tmp_knownj, tmp_dx, tmp_dy, tmp_latinc, 
                             tmp_loninc, tmp_lat, tmp_lon, tmp_loc,
							 &errstat, errmsg, ERRLEN);

    /* Terminate if there was an error */
	if (errstat != 0) {
		fprintf(stderr, errmsg);
		exit(errstat);
	}

/*
 * Coerce output back to float or double. What's returned is in
 * j,i order, so be sure to return i,j order.
 */
    coerce_output_float_or_double(loc,&tmp_loc[1],type_loc,1,i);
    coerce_output_float_or_double(loc,&tmp_loc[0],type_loc,1,i+npts);
  }

/*
 * Free unneeded memory.
 */
  if(type_lat       != NCL_double) NclFree(tmp_lat);
  if(type_lon       != NCL_double) NclFree(tmp_lon);
  if(type_truelat1  != NCL_double) NclFree(tmp_truelat1);
  if(type_truelat2  != NCL_double) NclFree(tmp_truelat2);
  if(type_stand_lon != NCL_double) NclFree(tmp_stand_lon);
  if(type_ref_lat   != NCL_double) NclFree(tmp_ref_lat);
  if(type_ref_lon   != NCL_double) NclFree(tmp_ref_lon);
  if(type_pole_lat  != NCL_double) NclFree(tmp_pole_lat);
  if(type_pole_lon  != NCL_double) NclFree(tmp_pole_lon);
  if(type_knowni    != NCL_double) NclFree(tmp_knowni);
  if(type_knownj    != NCL_double) NclFree(tmp_knownj);
  if(type_dx        != NCL_double) NclFree(tmp_dx);
  if(type_dy        != NCL_double) NclFree(tmp_dy);
  if(type_latinc    != NCL_double) NclFree(tmp_latinc);
  if(type_loninc    != NCL_double) NclFree(tmp_loninc);
  NclFree(tmp_loc);

  dim_info = malloc(sizeof(NclDimRec)*ndims_loc);
  if(dim_info == NULL) {
    NhlPError(NhlWARNING,NhlEUNKNOWN,"wrf_ll_to_ij: Unable to allocate memory for setting dimension names");
    return(NhlFATAL);
  }
  for(i = 0; i < ndims_loc; i++ ) {
    dim_info[i].dim_num   = i;
    dim_info[i].dim_quark = -1;
    dim_info[i].dim_size  = dsizes_loc[i];
  }
  dim_info[0].dim_quark = NrmStringToQuark("i_j_location");

/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)loc,
                            NULL,
                            ndims_loc,
                            dsizes_loc,
                            TEMPORARY,
                            NULL,
                            type_obj_loc
                            );
    
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          -1,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );

  NclFree(dsizes_loc);
  NclFree(dim_info);

/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);
}


NhlErrorTypes wrf_ij_to_ll_W( void )
{
/*
 * Input variables
 */
/*
 * Argument # 0
 */
  void *iloc;
  double *tmp_iloc = NULL;
  int ndims_iloc;
  ng_size_t dsizes_iloc[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_iloc;

/*
 * Argument # 1
 */
  void *jloc;
  double *tmp_jloc = NULL;
  int ndims_jloc;
  ng_size_t dsizes_jloc[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_jloc;

/*
 * Argument # 2
 */
  logical *opt;

/*
 * Variables for retrieving attributes from "opt".
 */
  NclAttList  *attr_list;
  NclAtt  attr_obj;
  NclStackEntry stack_entry;

/*
 * Variables that can be set via attributes.
 */
  int map_proj;
  void *truelat1 = NULL;
  void *truelat2 = NULL;
  void *stand_lon = NULL;
  void *ref_lat = NULL;
  void *ref_lon = NULL;
  void *pole_lat = NULL;
  void *pole_lon = NULL;
  void *knowni = NULL;
  void *knownj = NULL;
  void *dx = NULL;
  void *dy = NULL;
  void *latinc = NULL;
  void *loninc = NULL;

  double *tmp_truelat1 = NULL;
  double *tmp_truelat2 = NULL;
  double *tmp_stand_lon = NULL;
  double *tmp_ref_lat = NULL;
  double *tmp_ref_lon = NULL;
  double *tmp_pole_lat = NULL;
  double *tmp_pole_lon = NULL;
  double *tmp_knowni = NULL;
  double *tmp_knownj = NULL;
  double *tmp_dx = NULL;
  double *tmp_dy = NULL;
  double *tmp_latinc = NULL;
  double *tmp_loninc = NULL;

  NclBasicDataTypes type_truelat1 = NCL_none;
  NclBasicDataTypes type_truelat2 = NCL_none;
  NclBasicDataTypes type_stand_lon = NCL_none;
  NclBasicDataTypes type_ref_lat = NCL_none;
  NclBasicDataTypes type_ref_lon = NCL_none;
  NclBasicDataTypes type_pole_lat = NCL_none;
  NclBasicDataTypes type_pole_lon = NCL_none;
  NclBasicDataTypes type_knowni = NCL_none;
  NclBasicDataTypes type_knownj = NCL_none;
  NclBasicDataTypes type_dx = NCL_none;
  NclBasicDataTypes type_dy = NCL_none;
  NclBasicDataTypes type_latinc = NCL_none;
  NclBasicDataTypes type_loninc = NCL_none;

  logical set_map_proj, set_truelat1, set_truelat2, set_stand_lon, set_ref_lat;
  logical set_ref_lon, set_pole_lat, set_pole_lon,  set_knowni, set_knownj;
  logical set_dx, set_dy, set_latinc, set_loninc;

/*
 * Return variable
 */
  void *loc;
  double *tmp_loc;
  int ndims_loc;
  ng_size_t *dsizes_loc;
  NclBasicDataTypes type_loc;
  NclObjClass type_obj_loc;
/*
 * Variables for returning the output array with attributes attached.
 */
  NclMultiDValData return_md;
  NclVar tmp_var;
  NclStackEntry return_data;
/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info;

/*
 * Various
 */
  int npts, i;
  int errstat;
  char* errmsg;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
/*
 * Get argument # 0
 */
  iloc = (void*)NclGetArgValue(
           0,
           3,
           &ndims_iloc,
           dsizes_iloc,
           NULL,
           NULL,
           &type_iloc,
           DONT_CARE);

/*
 * Get argument # 1
 */
  jloc = (void*)NclGetArgValue(
           1,
           3,
           &ndims_jloc,
           dsizes_jloc,
           NULL,
           NULL,
           &type_jloc,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_jloc != ndims_iloc) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ij_to_ll: lat and lon must have the same number of dimensions");
    return(NhlFATAL);
  }

  for(i = 0; i < ndims_iloc; i++) {
    if(dsizes_jloc[i] != dsizes_iloc[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ij_to_ll: lat and lon must have the same dimensionality");
      return(NhlFATAL);
    }
  }

/*
 * Get argument # 2
 */
  opt = (logical*)NclGetArgValue(
           2,
           3,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);

/*
 * Calculate size of i/j dimensions.
 */
  npts = 1;
  for(i = 0; i < ndims_iloc; i++) npts *= dsizes_iloc[i];

/*
 * Start checking for attributes attached to "opt". Some are optional,
 * and some are not.  We'll check them later.
 */

  set_map_proj = set_truelat1 = set_truelat2 = set_stand_lon = False;
  set_ref_lat = set_ref_lon = set_pole_lat = set_pole_lon = False;
  set_knowni = set_knownj = set_dx = set_dy = set_latinc = set_loninc = False;

  stack_entry = _NclGetArg(2, 3, DONT_CARE);
  switch (stack_entry.kind) {
  case NclStk_VAR:
    if (stack_entry.u.data_var->var.att_id != -1) {
      attr_obj = (NclAtt) _NclGetObj(stack_entry.u.data_var->var.att_id);
      if (attr_obj == NULL) {
        break;
      }
    }
    else {
/*
 * att_id == -1 ==> no optional args given.
 */
      break;
    }
/* 
 * Get optional arguments.
 */
    if (attr_obj->att.n_atts > 0) {
/*
 * Get list of attributes.
 */
      attr_list = attr_obj->att.att_list;
/*
 * Loop through attributes and check them. We are looking for:
 *
 *   map_proj
 *   truelat1, truelat2
 *   stand_lon
 *   ref_lat, ref_lon
 *   pole_lat, pole_lon
 *   knowni, knownj
 *   dx, dy
 *   latinc, loninc
 */
      while (attr_list != NULL) {
        if(!strcasecmp(attr_list->attname, "map_proj")) {
          map_proj = *(int *)attr_list->attvalue->multidval.val;
          set_map_proj = True;
        }
        else if(!strcasecmp(attr_list->attname, "truelat1")) {
          truelat1      = attr_list->attvalue->multidval.val;
          type_truelat1 = attr_list->attvalue->multidval.data_type;
          set_truelat1  = True;
        }
        else if(!strcasecmp(attr_list->attname, "truelat2")) {
          truelat2      = attr_list->attvalue->multidval.val;
          type_truelat2 = attr_list->attvalue->multidval.data_type;
          set_truelat2  = True;
        }
        else if(!strcasecmp(attr_list->attname, "stand_lon")) {
          stand_lon      = attr_list->attvalue->multidval.val;
          type_stand_lon = attr_list->attvalue->multidval.data_type;
          set_stand_lon  = True;
        }
        else if(!strcasecmp(attr_list->attname, "ref_lat")) {
          ref_lat      = attr_list->attvalue->multidval.val;
          type_ref_lat = attr_list->attvalue->multidval.data_type;
          set_ref_lat  = True;
        }
        else if(!strcasecmp(attr_list->attname, "ref_lon")) {
          ref_lon      = attr_list->attvalue->multidval.val;
          type_ref_lon = attr_list->attvalue->multidval.data_type;
          set_ref_lon  = True;
        }
        else if(!strcasecmp(attr_list->attname, "pole_lat")) {
          pole_lat      = attr_list->attvalue->multidval.val;
          type_pole_lat = attr_list->attvalue->multidval.data_type;
          set_pole_lat  = True;
        }
        else if(!strcasecmp(attr_list->attname, "pole_lon")) {
          pole_lon      = attr_list->attvalue->multidval.val;
          type_pole_lon = attr_list->attvalue->multidval.data_type;
          set_pole_lon  = True;
        }
        else if(!strcasecmp(attr_list->attname, "knowni")) {
          knowni      = attr_list->attvalue->multidval.val;
          type_knowni = attr_list->attvalue->multidval.data_type;
          set_knowni  = True;
        }
        else if(!strcasecmp(attr_list->attname, "knownj")) {
          knownj      = attr_list->attvalue->multidval.val;
          type_knownj = attr_list->attvalue->multidval.data_type;
          set_knownj  = True;
        }
        else if(!strcasecmp(attr_list->attname, "dx")) {
          dx      = attr_list->attvalue->multidval.val;
          type_dx = attr_list->attvalue->multidval.data_type;
          set_dx  = True;
        }
        else if(!strcasecmp(attr_list->attname, "dy")) {
          dy      = attr_list->attvalue->multidval.val;
          type_dy = attr_list->attvalue->multidval.data_type;
          set_dy  = True;
        }
        else if(!strcasecmp(attr_list->attname, "latinc")) {
          latinc      = attr_list->attvalue->multidval.val;
          type_latinc = attr_list->attvalue->multidval.data_type;
          set_latinc  = True;
        }
        else if(!strcasecmp(attr_list->attname, "loninc")) {
          loninc      = attr_list->attvalue->multidval.val;
          type_loninc = attr_list->attvalue->multidval.data_type;
          set_loninc  = True;
        }
        attr_list = attr_list->next;
      }
    default:
      break;
    }
  }

/*
 * Check for attributes that need to be set, or set to a certain value.
 *
 * Check MAP_PROJ. Must be set.
 */
  if(!set_map_proj) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ij_to_ll: The MAP_PROJ attribute must be set");
    return(NhlFATAL);
  }
  else if(map_proj != 1 && map_proj != 2 && map_proj != 3 && map_proj != 6) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ij_to_ll: The MAP_PROJ attribute must be set to 1, 2, 3, or 6");
    return(NhlFATAL);
  }

/*
 * Check TRUELAT1. Must be set in some cases.
 */
  if( (map_proj == 1 || map_proj == 2 || map_proj == 3) && !set_truelat1) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ij_to_ll: The TRUELAT1 attribute must be set if MAP_PROJ is 1, 2, or 3");
    return(NhlFATAL);
  }
  if(set_truelat1) {
    tmp_truelat1 = coerce_input_double(truelat1,type_truelat1,1,0,NULL,NULL);
  }
  else {
    tmp_truelat1  = (double *)calloc(1,sizeof(double));
    *tmp_truelat1 = 0.;
  }

/*
 * Check TRUELAT2. Must be set in some cases.
 */
  if( map_proj == 1 && !set_truelat2) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ij_to_ll: The TRUELAT2 attribute must be set if MAP_PROJ is 1");
    return(NhlFATAL);
  }
  if(set_truelat2) {
    tmp_truelat2 = coerce_input_double(truelat2,type_truelat2,1,0,NULL,NULL);
  }
  else {
    tmp_truelat2  = (double *)calloc(1,sizeof(double));
    *tmp_truelat2 = 0.;
  }

/*
 * Check STAND_LON. Must be set.
 */
  if(!set_stand_lon) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ij_to_ll: The STAND_LON attribute must be set");
    return(NhlFATAL);
  }
  else {
    tmp_stand_lon = coerce_input_double(stand_lon,type_stand_lon,1,0,NULL,NULL);
  }

/*
 * Check REF_LAT/REF_LON. Must be set.
 */
  if(!set_ref_lat || !set_ref_lon) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ij_to_ll: The REF_LAT/REF_LON attributes must be set");
    return(NhlFATAL);
  }
  else {
    tmp_ref_lat = coerce_input_double(ref_lat,type_ref_lat,1,0,NULL,NULL);
    tmp_ref_lon = coerce_input_double(ref_lon,type_ref_lon,1,0,NULL,NULL);
  }

/*
 * Check POLE_LAT/POLE_LON.
 */
  if(set_pole_lat) {
    tmp_pole_lat = coerce_input_double(pole_lat,type_pole_lat,1,0,NULL,NULL);
  }
  else {
    tmp_pole_lat  = (double *)calloc(1,sizeof(double));
    *tmp_pole_lat = 90.;
  }

  if(set_pole_lon) {
    tmp_pole_lon = coerce_input_double(pole_lon,type_pole_lon,1,0,NULL,NULL);
  }
  else {
    tmp_pole_lon  = (double *)calloc(1,sizeof(double));
    *tmp_pole_lon = 0.;
  }

/*
 * Check KNOWNI/KNOWNJ. Must be set.
 */
  if(!set_knowni || !set_knownj) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ij_to_ll: The KNOWNI/KNOWNJ attributes must be set");
    return(NhlFATAL);
  }
  else {
    tmp_knowni = coerce_input_double(knowni,type_knowni,1,0,NULL,NULL);
    tmp_knownj = coerce_input_double(knownj,type_knownj,1,0,NULL,NULL);
  }

/*
 * Check DX/DY. Must be set in some cases.
 */
  if( (map_proj == 1 || map_proj == 2 || map_proj == 3) &&
      (!set_dx || !set_dy)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ij_to_ll: The DX/DY attributes must be set if MAP_PROJ is 1, 2, or 3");
    return(NhlFATAL);
  }
  if(set_dx) {
    tmp_dx = coerce_input_double(dx,type_dx,1,0,NULL,NULL);
  }
  else {
    tmp_dx  = (double *)calloc(1,sizeof(double));
    *tmp_dx = 0.;
  }

  if(set_dy) {
    tmp_dy = coerce_input_double(dy,type_dy,1,0,NULL,NULL);
  }
  else {
    tmp_dy  = (double *)calloc(1,sizeof(double));
    *tmp_dy = 0.;
  }

/*
 * Check LATINC/LONINC. Must be set in some cases.
 */
  if( map_proj == 6 && (!set_latinc || !set_loninc)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ij_to_ll: The LATINC/LONINC attributes must be set if MAP_PROJ is 6");
    return(NhlFATAL);
  }
  if(set_latinc) {
    tmp_latinc = coerce_input_double(latinc,type_latinc,1,0,NULL,NULL);
  }
  else {
    tmp_latinc  = (double *)calloc(1,sizeof(double));
    *tmp_latinc = 0.;
  }

  if(set_loninc) {
    tmp_loninc = coerce_input_double(loninc,type_loninc,1,0,NULL,NULL);
  }
  else {
    tmp_loninc  = (double *)calloc(1,sizeof(double));
    *tmp_loninc = 0.;
  }

/*
 * The output type defaults to float, unless either of the lat/lon arrays
 * are double.
 */
  type_loc     = NCL_float;
  type_obj_loc = nclTypefloatClass;

/* 
 * Allocate space for coercing input arrays.  If any of the input
 * is already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 */
/*
 * Allocate space for tmp_iloc.
 */
  if(type_iloc != NCL_double) {
    tmp_iloc = (double *)calloc(npts,sizeof(double));
    if(tmp_iloc == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ij_to_ll: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_loc     = NCL_double;
    type_obj_loc = nclTypedoubleClass;
  }

/*
 * Allocate space for tmp_jloc.
 */
  if(type_jloc != NCL_double) {
    tmp_jloc = (double *)calloc(npts,sizeof(double));
    if(tmp_jloc == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ij_to_ll: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_loc     = NCL_double;
    type_obj_loc = nclTypedoubleClass;
  }

/* 
 * Allocate space for output array.
 */
  tmp_loc = (double *)calloc(2,sizeof(double));
  if(type_loc != NCL_double) loc = (void *)calloc(2*npts, sizeof(float));
  else                       loc = (void *)calloc(2*npts, sizeof(double));
  if(loc == NULL || tmp_loc == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ij_to_ll: Unable to allocate memory for output array");
    return(NhlFATAL);
  }

/* 
 * Allocate space for output dimension sizes and set them.
 */
  if(is_scalar(ndims_iloc,dsizes_iloc)) {
    ndims_loc = 1;
  }
  else {
    ndims_loc = ndims_iloc + 1;
  }
  dsizes_loc = (ng_size_t*)calloc(ndims_loc,sizeof(ng_size_t));  
  if( dsizes_loc == NULL ) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ij_to_ll: Unable to allocate memory for holding dimension sizes");
    return(NhlFATAL);
  }
  for(i = 0; i < ndims_loc-1; i++) dsizes_loc[i+1] = dsizes_iloc[i];
  dsizes_loc[0] = 2;

  /* Allocate space for errmsg*/
  errmsg = (char *) calloc(ERRLEN, sizeof(char))

/*
 * Loop across all lat/lon points and call the Fortran routine for each
 * point.
 */
  for(i = 0; i < npts; i++) {
/*
 * Coerce subsection of iloc (tmp_iloc) to double if necessary.
 */
    if(type_iloc != NCL_double) {
      coerce_subset_input_double(iloc,tmp_iloc,i,type_iloc,1,0,NULL,NULL);
    }
    else {
      tmp_iloc = &((double*)iloc)[i];
    }

/*
 * Coerce subsection of lon (tmp_jloc) to double if necessary.
 */
    if(type_jloc != NCL_double) {
      coerce_subset_input_double(jloc,tmp_jloc,i,type_jloc,1,0,NULL,NULL);
    }
    else {
      tmp_jloc = &((double*)jloc)[i];
    }

/*
 * Call the Fortran routine.
 */
    errstat = 0;
    errmsg = "";
    NGCALLF(dijtoll,DIJTOLL)(&map_proj, tmp_truelat1, tmp_truelat2, 
                             tmp_stand_lon, tmp_ref_lat, tmp_ref_lon, 
                             tmp_pole_lat, tmp_pole_lon, tmp_knowni,
                             tmp_knownj, tmp_dx, tmp_dy, tmp_latinc, 
                             tmp_loninc, tmp_iloc, tmp_jloc, tmp_loc,
							 &errstat, errmsg, ERRLEN);

    /* Terminate if there was an error */
	if (errstat != 0) {
		fprintf(stderr, errmsg);
		exit(errstat);
	}

/*
 * Coerce output back to float or double. What's returned is in
 * lat,lon (j,i) order, so be sure to return lon,lat (i,j) order.
 */
    coerce_output_float_or_double(loc,&tmp_loc[1],type_loc,1,i);
    coerce_output_float_or_double(loc,&tmp_loc[0],type_loc,1,i+npts);
  }

/*
 * Free unneeded memory.
 */
  if(type_iloc      != NCL_double) NclFree(tmp_iloc);
  if(type_jloc      != NCL_double) NclFree(tmp_jloc);
  if(type_truelat1  != NCL_double) NclFree(tmp_truelat1);
  if(type_truelat2  != NCL_double) NclFree(tmp_truelat2);
  if(type_stand_lon != NCL_double) NclFree(tmp_stand_lon);
  if(type_ref_lat   != NCL_double) NclFree(tmp_ref_lat);
  if(type_ref_lon   != NCL_double) NclFree(tmp_ref_lon);
  if(type_pole_lat  != NCL_double) NclFree(tmp_pole_lat);
  if(type_pole_lon  != NCL_double) NclFree(tmp_pole_lon);
  if(type_knowni    != NCL_double) NclFree(tmp_knowni);
  if(type_knownj    != NCL_double) NclFree(tmp_knownj);
  if(type_dx        != NCL_double) NclFree(tmp_dx);
  if(type_dy        != NCL_double) NclFree(tmp_dy);
  if(type_latinc    != NCL_double) NclFree(tmp_latinc);
  if(type_loninc    != NCL_double) NclFree(tmp_loninc);
  NclFree(tmp_loc);

  dim_info = malloc(sizeof(NclDimRec)*ndims_loc);
  if(dim_info == NULL) {
    NhlPError(NhlWARNING,NhlEUNKNOWN,"wrf_ij_to_ll: Unable to allocate memory for setting dimension names");
    return(NhlFATAL);
  }
  for(i = 0; i < ndims_loc; i++ ) {
    dim_info[i].dim_num   = i;
    dim_info[i].dim_quark = -1;
    dim_info[i].dim_size  = dsizes_loc[i];
  }
  dim_info[0].dim_quark = NrmStringToQuark("lon_lat_location");

/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)loc,
                            NULL,
                            ndims_loc,
                            dsizes_loc,
                            TEMPORARY,
                            NULL,
                            type_obj_loc
                            );
    
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          -1,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );

  NclFree(dsizes_loc);
  NclFree(dim_info);

/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);
}


/*
 * Function for calculating cape (from the RIP code). This function
 * depends on the "psadilookup.dat" file, which by default will be
 * searched for in $NCARG_ROOT/lib/ncarg/data/asc/), unless
 * NCARG_PSADILOOKUP is set to the location of this file.
 */

/*
 * The wrf_cape_3d wrapper is for the case where I3DFLAG is set to
 * 1 in the Fortran rip_cape.f file. This wrapper is similar to
 * rip_cape_3d except it will flip the first four input arrays if
 * the pressure values are not decreasing. It will also multiple the
 * pressure values by 0.01 to convert from hPa to Pa.
 * 
 */
NhlErrorTypes wrf_cape_3d_W( void )
{
/*
 * Input array variables
 */
  void *p, *t, *q, *z, *zsfc, *psfc;
  logical *ter_follow;
  double *tmp_p = NULL;
  double *tmp_t = NULL;
  double *tmp_q = NULL;
  double *tmp_z = NULL;
  double *tmp_zsfc = NULL;
  double *tmp_psfc = NULL;
  double *tmp_p_orig = NULL;
  double *tmp_t_orig = NULL;
  double *tmp_q_orig = NULL;
  double *tmp_z_orig = NULL;
  int ndims_p, ndims_t, ndims_q, ndims_z, ndims_zsfc, ndims_psfc;
  ng_size_t dsizes_p[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_t[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_q[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_z[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_zsfc[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_psfc[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_p, type_t, type_q, type_z, type_zsfc, type_psfc;

/*
 * Output array variables
 */
  void *cape;
  double *tmp_cape_orig, *tmp_cin_orig, cmsg;
  double *tmp_cape, *tmp_cin;
  NclBasicDataTypes type_cape;
  NclObjClass type_obj_cape;
  int ndims_cape;
  NclScalar missing_cape;
  ng_size_t *dsizes_cape;

/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info = NULL;
  NclDimRec *dim_info_t;

/*
 * Variables for returning the output array with attributes and/or
 * dimension names attached.
 */
  NclMultiDValData return_md;
  NclVar tmp_var;
  NclStackEntry return_data;
/*
 * Declare various variables for random purposes.
 */
  ng_size_t i;
  ng_size_t miy = 0;
  ng_size_t mjx = 0;
  ng_size_t mkzh = 0;
  ng_size_t ntime = 0;
  ng_size_t nz = 0;
  ng_size_t size_cape, size_output, size_zsfc;
  int i3dflag=1, scalar_zsfc;
  ng_size_t index_cape, index_zsfc, index_cin;
  int iter;
  logical flip;
  int imiy, imjx, imkzh;
  char *psa_file;
  int errstat;
  char *errmsg;

/*
 * Retrieve parameters
 *
 * Note that any of the pointer parameters can be set to NULL,
 * which implies you don't care about its value.
 *
 */
  p = (void*)NclGetArgValue(
          0,
          7,
          &ndims_p,
          dsizes_p,
          NULL,
          NULL,
          &type_p,
          DONT_CARE);

  t = (void*)NclGetArgValue(
          1,
          7,
          &ndims_t,
          dsizes_t,
          NULL,
          NULL,
          &type_t,
          DONT_CARE);


  q = (void*)NclGetArgValue(
          2,
          7,
          &ndims_q,
          dsizes_q,
          NULL,
          NULL,
          &type_q,
          DONT_CARE);

  z = (void*)NclGetArgValue(
          3,
          7,
          &ndims_z,
          dsizes_z,
          NULL,
          NULL,
          &type_z,
          DONT_CARE);

  zsfc = (void*)NclGetArgValue(
          4,
          7,
          &ndims_zsfc,
          dsizes_zsfc,
          NULL,
          NULL,
          &type_zsfc,
          DONT_CARE);

  psfc = (void*)NclGetArgValue(
          5,
          7,
          &ndims_psfc,
          dsizes_psfc,
          NULL,
          NULL,
          &type_psfc,
          DONT_CARE);

  ter_follow = (logical*)NclGetArgValue(
          6,
          7,
          NULL,
          NULL,
          NULL,
          NULL,
          NULL,
          DONT_CARE);
  
  if(*ter_follow) iter = 1;
  else            iter = 0;

/*
 * Check the input dimension sizes. There are four possible cases
 * for the input dimension sizes:
 *
 *  - p,t,q,z (nz,time,lev,lat,lon) and psfc,zsfc (nz,time,lat,lon)
 *  - p,t,q,z (time,lev,lat,lon) and psfc,zsfc (time,lat,lon)
 *  - p,t,q,z (lev,lat,lon) and psfc,zsfc (lat,lon)
 *  - p,t,q,z (lev) and psfc,zsfc (scalars)
 */
  if(ndims_p != ndims_t || ndims_p != ndims_q || ndims_p != ndims_z) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_3d: The p, t, q, and z arrays must all have the same number of dimensions");
    return(NhlFATAL);
  }
  if(ndims_p != 1 && ndims_p != 3 && ndims_p != 4 && ndims_p != 5) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_3d: The p, t, q, and z arrays must be 1-, 3-, 4-, or 5-dimensional\n");
    return(NhlFATAL);
  }
/*
 * zsfc and psfc can be scalars, if the other input arrays are 1D.
 */
  scalar_zsfc = is_scalar(ndims_zsfc,dsizes_zsfc);

  if((ndims_zsfc != ndims_psfc) || (scalar_zsfc && ndims_p != 1) || 
     (!scalar_zsfc && ndims_zsfc != ndims_p-1)) { 
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_3d: The zsfc and psfc arrays must have the same number of dimensions, and either be scalars or one less dimension than the other input arrays");
    return(NhlFATAL);
  }

/*
 * Now check that the dimension sizes are equal to each other.
 */
  for(i = 0; i < ndims_p; i++) {
    if(dsizes_p[i] != dsizes_t[i] || dsizes_p[i] != dsizes_q[i] || 
       dsizes_p[i] != dsizes_z[i]) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_3d: p, t, q, and z must be the same dimensionality");
    return(NhlFATAL);
    }
  }

  for(i = 0; i < ndims_psfc; i++) {
    if(dsizes_psfc[i] != dsizes_zsfc[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_3d: psfc and zsfc must be the same dimensionality");
      return(NhlFATAL);
    }
  }
/*
 * Get sizes of input arrays.
 */
  if(ndims_p == 5) {
    nz    = dsizes_p[0];       /* nz, serves as leftmost dimension */
    ntime = dsizes_p[1];       /* time, also serves as leftmost dimension */
    mkzh  = dsizes_p[2];       /* lev */
    mjx   = dsizes_p[3];       /* lat */
    miy   = dsizes_p[4];       /* lon */
  }
  else if(ndims_p == 4) {
    nz    = 1;
    ntime = dsizes_p[0];       /* time, serves as a leftmost dimension */
    mkzh  = dsizes_p[1];       /* lev */
    mjx   = dsizes_p[2];       /* lat */
    miy   = dsizes_p[3];       /* lon */
  }
  else if(ndims_p == 3) {
    nz    = 1;
    ntime = 1;
    mkzh  = dsizes_p[0];       /* lev */
    mjx   = dsizes_p[1];       /* lat */
    miy   = dsizes_p[2];       /* lon */
  }
  else if(ndims_p == 1) {
    nz    = 1;
    ntime = 1;
    mkzh  = dsizes_p[0];       /* lev */
    mjx   = 1;                 /* lat */
    miy   = 1;                 /* lon */
  }

/*
 * Test input dimension sizes.
 */
  if((miy > INT_MAX) || (mjx > INT_MAX) || (mkzh > INT_MAX)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_3d: one or more dimension sizes is greater than INT_MAX");
    return(NhlFATAL);
  }
  imiy = (int) miy;
  imjx = (int) mjx;
  imkzh = (int) mkzh;

/*
 * Check some more dimension sizes.
 */
  if(ndims_p == 5) {
    if(dsizes_psfc[0] != nz || dsizes_psfc[1] != ntime || 
       dsizes_psfc[2] != mjx || dsizes_psfc[3] != miy) { 
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_3d: If p,q,t,z are 4-dimensional (time x lev x lat x lon), psfc,zsfc must be 3-dimensional (time x lat x lon)");
      return(NhlFATAL);
    }
  }
  else if(ndims_p == 4) {
    if(dsizes_psfc[0] != ntime || dsizes_psfc[1] != mjx || 
       dsizes_psfc[2] != miy) { 
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_3d: If p,q,t,z are 4-dimensional (time x lev x lat x lon), psfc,zsfc must be 3-dimensional (time x lat x lon)");
      return(NhlFATAL);
    }
  }
  else if(ndims_p == 3) {
    if(dsizes_psfc[0] != mjx || dsizes_psfc[1] != miy) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_3d: If p,q,t,z are 3-dimensional (time x lev x lat x lon), psfc,zsfc must be 2-dimensional (lat x lon)");
      return(NhlFATAL);
    }
  }
/*
 * Calculate size of output array. The output array size depends on
 * the size of p,t,q,z:
 *
 *  - p,t,q,z (nz,time,lev,lat,lon) and psfc,zsfc (nz,time,lat,lon)
 *       output array: (2,nz,time,lev,lat,lon)
 *  - p,t,q,z (time,lev,lat,lon) and psfc,zsfc (time,lat,lon)
 *       output array: (2,time,lev,lat,lon)
 *  - p,t,q,z (lev,lat,lon) and psfc,zsfc (lat,lon)
 *       output array: (2,lev,lat,lon)
 *  - p,t,q,z (lev) and psfc,zsfc (scalars)
 *       output array: (2,lev)
 */
  ndims_cape = ndims_p+1;
  dsizes_cape = (ng_size_t *)calloc(ndims_cape,sizeof(ng_size_t));
  if(dsizes_cape == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_3d: Unable to allocate memory for array dimensionality");
    return(NhlFATAL);
  }

  dsizes_cape[0] = 2;                /* 0 = cape, 1 = cin */
  for(i = 0; i < ndims_p; i++ ) {
    dsizes_cape[i+1] = dsizes_p[i];
  }
  size_zsfc   = mjx * miy;
  size_cape   = mkzh * size_zsfc;       /* Also size of cin array */
  size_output = 2 * size_cape * ntime * nz;

/* 
 * Allocate space for output arrays. We are allocating space for 
 * tmp_cape_orig and tmp_cin_orig even if the output will be double,
 * because we may also need to flip the values before we're done.
 *
 * The addition of missing values was added in V6.1.0.
 */
  if(type_p == NCL_double || type_t == NCL_double || type_q == NCL_double ||
     type_z == NCL_double) {
    type_cape     = NCL_double;
    type_obj_cape = nclTypedoubleClass;
    cape = (double *)calloc(size_output,sizeof(double));
    missing_cape.doubleval = ((NclTypeClass)nclTypedoubleClass)->type_class.default_mis.doubleval;
    cmsg = missing_cape.doubleval;
  }
  else {
    type_cape     = NCL_float;
    type_obj_cape = nclTypefloatClass;
    cape          = (float *)calloc(size_output,sizeof(float));
    missing_cape.floatval = ((NclTypeClass)nclTypefloatClass)->type_class.default_mis.floatval;
    cmsg = (double)missing_cape.floatval;
  }
  tmp_cape_orig = (double *)calloc(size_cape,sizeof(double));
  tmp_cin_orig  = (double *)calloc(size_cape,sizeof(double));
  if(cape == NULL || tmp_cape_orig == NULL || tmp_cin_orig == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_3d: Unable to allocate memory for output arrays");
    return(NhlFATAL);
  }

/*
 * Allocate memory for coercing input arrays to double, if necessary.
 * Force a copy of variable p, because we need to multiply it by 0.01,
 * and we don't want this to propagate back to the NCL script.
 */
  tmp_p_orig = (double *)calloc(size_cape,sizeof(double));
  if(tmp_p_orig == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_3d: Unable to allocate memory for coercing input arrays to double");
    return(NhlFATAL);
  }

  if(type_t != NCL_double) {
    tmp_t_orig = (double *)calloc(size_cape,sizeof(double));
    if(tmp_t_orig == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_3d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  if(type_q != NCL_double) {
    tmp_q_orig = (double *)calloc(size_cape,sizeof(double));
    if(tmp_q_orig == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_3d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  if(type_z != NCL_double) {
    tmp_z_orig = (double *)calloc(size_cape,sizeof(double));
    if(tmp_z_orig == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_3d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  if(type_zsfc != NCL_double) {
    tmp_zsfc = (double *)calloc(size_zsfc,sizeof(double));
    if(tmp_zsfc == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_3d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  if(type_psfc != NCL_double) {
    tmp_psfc = (double *)calloc(size_zsfc,sizeof(double));
    if(tmp_psfc == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_3d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

/*
 * We need to coerce the pressure array once outside the loop to
 * check if the values are in ascending order.
 *
 * If not, we need to flip the leftmost dimension (p = p(::-1,:,:) in
 * NCL-ese), *and* flip the other 3 input arrays in the same fashion.
 */
  coerce_subset_input_double(p,tmp_p_orig,0,type_p,size_cape,0,NULL,NULL);

  if(tmp_p_orig[0] > tmp_p_orig[(mkzh-1)*size_zsfc] ) {
    flip     = True;
    tmp_p    = (double *)calloc(size_cape,sizeof(double));
    tmp_t    = (double *)calloc(size_cape,sizeof(double));
    tmp_q    = (double *)calloc(size_cape,sizeof(double));
    tmp_z    = (double *)calloc(size_cape,sizeof(double));
    tmp_cape = (double *)calloc(size_cape,sizeof(double));
    tmp_cin  = (double *)calloc(size_cape,sizeof(double));

    if(tmp_p == NULL || tmp_t == NULL || tmp_q == NULL || tmp_z == NULL ||
       tmp_cape == NULL || tmp_cin == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_3d: Unable to allocate memory for flipping arrays");
      return(NhlFATAL);
    }
  }
  else  {
    flip     = False;
    tmp_p    = tmp_p_orig;
    tmp_t    = tmp_t_orig;
    tmp_q    = tmp_q_orig;
    tmp_z    = tmp_z_orig;
    tmp_cape = tmp_cape_orig;
    tmp_cin  = tmp_cin_orig;
  }

/*
 * Get path to psadilookup.dat file required by this routine. 
 */
  psa_file = get_psa_file();

/* Allocate space for errmsg*/
  errmsg = (char *) calloc(ERRLEN, sizeof(char))

/*
 * Loop through time,nz and call the Fortran routine.
 */ 
  index_cape = index_zsfc = 0;
  index_cin  = ntime * nz * size_cape;    /* Second half of output array */

  for(i = 0; i < ntime*nz; i++) {
/*
 * Coerce subset of input arrays to double if necessary.
 */
    if(i > 0) {
      coerce_subset_input_double(p,tmp_p_orig,index_cape,type_p,
                                 size_cape,0,NULL,NULL);
    }
/*
 * Multiple pressure values by 0.01 to convert from Pa to hPa.
 * The assumption is that pressure values come in as Pa.
 */
    convert_to_hPa(tmp_p_orig,size_cape);

    if(type_t != NCL_double) {
      coerce_subset_input_double(t,tmp_t_orig,index_cape,type_t,
                                 size_cape,0,NULL,NULL);
    }
    else {
/*
 * Point tmp_t_orig to appropriate location in t.
 */
      tmp_t_orig = &((double*)t)[index_cape];
    }
    if(type_q != NCL_double) {
      coerce_subset_input_double(q,tmp_q_orig,index_cape,type_q,
                                 size_cape,0,NULL,NULL);
    }
    else {
/*
 * Point tmp_q_orig to appropriate location in q.
 */
      tmp_q_orig = &((double*)q)[index_cape];
    }
    if(type_z != NCL_double) {
      coerce_subset_input_double(z,tmp_z_orig,index_cape,type_z,
                                 size_cape,0,NULL,NULL);
    }
    else {
/*
 * Point tmp_z_orig to appropriate location in z.
 */
      tmp_z_orig = &((double*)z)[index_cape];
    }

    if(type_psfc != NCL_double) {
      coerce_subset_input_double(psfc,tmp_psfc,index_zsfc,type_psfc,
                                 size_zsfc,0,NULL,NULL);
    }
    else {
/*
 * Point tmp_psfc to appropriate location in psfc.
 */
      tmp_psfc = &((double*)psfc)[index_zsfc];
    }
    if(type_zsfc != NCL_double) {
      coerce_subset_input_double(zsfc,tmp_zsfc,index_zsfc,type_zsfc,
                                 size_zsfc,0,NULL,NULL);
    }
    else {
/*
 * Point tmp_zsfc to appropriate location in zsfc.
 */
      tmp_zsfc = &((double*)zsfc)[index_zsfc];
    }
    
/*
 * If the pressure values need to be flipped, we also need to flip
 * the z, q, and t values in the same fashion.
 */
    if(flip) {
      flip_it(tmp_p_orig,tmp_p,mkzh,size_zsfc);
      flip_it(tmp_t_orig,tmp_t,mkzh,size_zsfc);
      flip_it(tmp_q_orig,tmp_q,mkzh,size_zsfc);
      flip_it(tmp_z_orig,tmp_z,mkzh,size_zsfc);
    }
    else {
      tmp_p = tmp_p_orig;
      tmp_t = tmp_t_orig;
      tmp_q = tmp_q_orig;
      tmp_z = tmp_z_orig;
    }

/*
 * Call Fortran routine.
 */
    NGCALLF(dcapecalc3d,DCAPECALC3D)(tmp_p, tmp_t, tmp_q, tmp_z, tmp_zsfc,
                                     tmp_psfc, tmp_cape_orig, tmp_cin_orig,
                                     &cmsg,&imiy, &imjx, &imkzh, &i3dflag,
                                     &iter,psa_file,&errstat,errmsg,
									 strlen(psa_file),ERRLEN);

/* Terminate if there was an error */
    if (errstat != 0) {
    	fprintf(stderr, errmsg);
    	exit(errstat);
    }

/*
 * If we flipped arrays before going into the Fortran routine, we need
 * to flip the output values as well.
 */
    if(flip) {
      flip_it(tmp_cape_orig,tmp_cape,mkzh,size_zsfc);
      flip_it(tmp_cin_orig,tmp_cin,mkzh,size_zsfc);
    }
    else {
      tmp_cape = tmp_cape_orig;
      tmp_cin  = tmp_cin_orig;
    }
/*
 * If the output is to be float, then do the coercion here.
 */
    coerce_output_float_or_double(cape,tmp_cape,type_cape,size_cape,index_cape);
    coerce_output_float_or_double(cape,tmp_cin,type_cape,size_cape,index_cin);

/*
 * Implement the pointers into the arrays.
 */
    index_cape += size_cape;
    index_cin  += size_cape;
    index_zsfc += size_zsfc;
  }
/*
 * Free memory.
 */
  NclFree(tmp_p_orig);
  NclFree(tmp_cape_orig);
  NclFree(tmp_cin_orig);
  if(type_t != NCL_double) NclFree(tmp_t_orig);
  if(type_q != NCL_double) NclFree(tmp_q_orig);
  if(type_z != NCL_double) NclFree(tmp_z_orig);
  if(type_zsfc != NCL_double) NclFree(tmp_zsfc);
  if(type_psfc != NCL_double) NclFree(tmp_psfc);
  if(flip) {
    NclFree(tmp_p);
    NclFree(tmp_t);
    NclFree(tmp_q);
    NclFree(tmp_z);
    NclFree(tmp_cape);
    NclFree(tmp_cin);
  }
  NclFree(psa_file);
  NclFree(errmsg);

/*
 * Get dimension info to see if we have named dimensions.
 * This will be used for return variable.
 */
  dim_info_t = get_wrf_dim_info(1,7,ndims_t,dsizes_t);
  if(dim_info_t != NULL) {
    dim_info = malloc(sizeof(NclDimRec)*ndims_cape);
    if(dim_info == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_3d: Unable to allocate memory for holding dimension information");
      return(NhlFATAL);
    }
    for(i = 0; i < ndims_cape; i++ ) {
      dim_info[i].dim_num  = i;
      dim_info[i].dim_size = dsizes_cape[i];
      if(i != 0) dim_info[i].dim_quark = dim_info_t[i-1].dim_quark;
      else       dim_info[0].dim_quark = NrmStringToQuark("cape_cin");
    }
  }


/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)cape,
                            &missing_cape,
                            ndims_cape,
                            dsizes_cape,
                            TEMPORARY,
                            NULL,
                            type_obj_cape
                            );
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          -1,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );

  NclFree(dsizes_cape);
  if(dim_info   != NULL) NclFree(dim_info);
  NclFree(dim_info_t);

/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);
}


/*
 * The wrf_cape_2d wrapper is for the case where I3DFLAG is set to
 * 0 in the Fortran rip_cape.f file.  In this case, 4 2D arrays
 * are returned: cape, cin, lcl, and lfc, but they are all returned 
 * in one big array whose leftmost dimension is 4:
 *
 * This wrapper is similar to rip_cape_2d except it will flip the first
 * four input arrays if the pressure values are not decreasing. 
 * It will also multiple the pressure values by 0.01 to convert 
 * from hPa to Pa.
 *
 *   index 0 = cape
 *   index 1 = cin
 *   index 2 = lcl
 *   index 3 = lfc
 */
NhlErrorTypes wrf_cape_2d_W( void )
{
/*
 * Input array variables
 */
  void *p, *t, *q, *z, *zsfc, *psfc;
  logical *ter_follow;
  double *tmp_p = NULL;
  double *tmp_t = NULL;
  double *tmp_q = NULL;
  double *tmp_z = NULL;
  double *tmp_zsfc = NULL;
  double *tmp_psfc = NULL;
  double *tmp_p_orig = NULL;
  double *tmp_t_orig = NULL;
  double *tmp_q_orig = NULL;
  double *tmp_z_orig = NULL;
  int ndims_p, ndims_t, ndims_q, ndims_z, ndims_zsfc, ndims_psfc;
  ng_size_t dsizes_p[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_t[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_q[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_z[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_zsfc[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_psfc[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_p, type_t, type_q, type_z, type_zsfc, type_psfc;

/*
 * Output array variables
 */
  void *cape;
  double *tmp_cape, *tmp_cin, cmsg;
  NclBasicDataTypes type_cape;
  NclObjClass type_obj_cape;
  int ndims_cape = 0;
  NclScalar missing_cape;
  ng_size_t *dsizes_cape;

/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info = NULL;
  NclDimRec *dim_info_t;

/*
 * Variables for returning the output array with attributes and/or
 * dimension names attached.
 */
  NclMultiDValData return_md;
  NclVar tmp_var;
  NclStackEntry return_data;
/*
 * Declare various variables for random purposes.
 */
  ng_size_t i;
  ng_size_t miy = 0;
  ng_size_t mjx = 0;
  ng_size_t mkzh = 0;
  ng_size_t ntime = 0;
  ng_size_t nz = 0;
  ng_size_t size_cape, size_output, size_zsfc;
  ng_size_t size_left_zsfc;
  ng_size_t index_cape, index_zsfc;
  ng_size_t index_output_cape, index_output_cin, index_output_lcl;
  ng_size_t index_output_lfc, mkzh0_index, mkzh1_index, mkzh2_index;
  int imiy, imjx, imkzh, iter, i3dflag=0;
  char *psa_file;
  logical flip;
  int errstat;
  char* errmsg;

/*
 * Retrieve parameters
 *
 * Note that any of the pointer parameters can be set to NULL,
 * which implies you don't care about its value.
 *
 */
  p = (void*)NclGetArgValue(
          0,
          7,
          &ndims_p,
          dsizes_p,
          NULL,
          NULL,
          &type_p,
          DONT_CARE);

  t = (void*)NclGetArgValue(
          1,
          7,
          &ndims_t,
          dsizes_t,
          NULL,
          NULL,
          &type_t,
          DONT_CARE);


  q = (void*)NclGetArgValue(
          2,
          7,
          &ndims_q,
          dsizes_q,
          NULL,
          NULL,
          &type_q,
          DONT_CARE);

  z = (void*)NclGetArgValue(
          3,
          7,
          &ndims_z,
          dsizes_z,
          NULL,
          NULL,
          &type_z,
          DONT_CARE);

  zsfc = (void*)NclGetArgValue(
          4,
          7,
          &ndims_zsfc,
          dsizes_zsfc,
          NULL,
          NULL,
          &type_zsfc,
          DONT_CARE);

  psfc = (void*)NclGetArgValue(
          5,
          7,
          &ndims_psfc,
          dsizes_psfc,
          NULL,
          NULL,
          &type_psfc,
          DONT_CARE);

  ter_follow = (logical*)NclGetArgValue(
          6,
          7,
          NULL,
          NULL,
          NULL,
          NULL,
          NULL,
          DONT_CARE);
  
  if(*ter_follow) iter = 1;
  else            iter = 0;


/*
 * Check the input dimension sizes. There are three possible cases
 * for the input dimension sizes:
 *
 *  - p,t,q,z (nz,time,lev,lat,lon) and psfc,zsfc (nz,time,lat,lon)
 *  - p,t,q,z (time,lev,lat,lon) and psfc,zsfc (time,lat,lon)
 *  - p,t,q,z (lev,lat,lon) and psfc,zsfc (lat,lon)
 */
  if(ndims_p != ndims_t || ndims_p != ndims_q || ndims_p != ndims_z) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: The p, t, q, and z arrays must all have the same number of dimensions");
    return(NhlFATAL);
  }
  if(ndims_p != 3 && ndims_p != 4 && ndims_p != 5) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: The p, t, q, and z arrays must be 3-, 4- or 5-dimensional\n");
    return(NhlFATAL);
  }
/*
 * Check zsfc and psfc dimension sizes.
 */
  if((ndims_zsfc != ndims_psfc) || (ndims_zsfc != ndims_p-1)) { 
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: The zsfc and psfc arrays must have the same number of dimensions and be one less dimension than the other input arrays");
    return(NhlFATAL);
  }

/*
 * Now check that the dimension sizes are equal to each other.
 */
  for(i = 0; i < ndims_p; i++) {
    if(dsizes_p[i] != dsizes_t[i] || dsizes_p[i] != dsizes_q[i] || 
       dsizes_p[i] != dsizes_z[i]) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: p, t, q, and z must be the same dimensionality");
    return(NhlFATAL);
    }
  }

  for(i = 0; i < ndims_psfc; i++) {
    if(dsizes_psfc[i] != dsizes_zsfc[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: psfc and zsfc must be the same dimensionality");
      return(NhlFATAL);
    }
  }
  if(ndims_p == 5) {
/*
 * Store dimension sizes.
 */
    nz    = dsizes_p[0];        /* nz */
    ntime = dsizes_p[1];        /* time */
    mkzh  = dsizes_p[2];        /* lev */
    mjx   = dsizes_p[3];        /* lat */
    miy   = dsizes_p[4];        /* lon */
    ndims_cape = 5;
    if(dsizes_psfc[0] != nz || dsizes_psfc[1] != ntime || 
       dsizes_psfc[2] != mjx || dsizes_psfc[3] != miy) { 
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: If p,q,t,z are 5-dimensional (nz x time x lev x lat x lon), psfc,zsfc must be 4-dimensional (nz x time x lat x lon)");
      return(NhlFATAL);

    }
  }
  else if(ndims_p == 4) {
/*
 * Store dimension sizes.
 */
    nz    = 1;
    ntime = dsizes_p[0];        /* time */
    mkzh  = dsizes_p[1];        /* lev */
    mjx   = dsizes_p[2];        /* lat */
    miy   = dsizes_p[3];        /* lon */
    ndims_cape = 4;
    if(dsizes_psfc[0] != ntime || dsizes_psfc[1] != mjx ||
       dsizes_psfc[2] != miy) { 
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: If p,q,t,z are 4-dimensional (time x lev x lat x lon), psfc,zsfc must be 3-dimensional (time x lat x lon)");
      return(NhlFATAL);

    }
  }
  else if(ndims_p == 3) {
/*
 * Store dimension sizes.
 */
    nz    = 1;
    ntime = 1;
    mkzh = dsizes_p[0];           /* lev */
    mjx  = dsizes_p[1];           /* lat */
    miy  = dsizes_p[2];           /* lon */
    ndims_cape = 3;
    if(dsizes_psfc[0] != mjx || dsizes_psfc[1] != miy) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: If p,q,t,z are 3-dimensional (time x lev x lat x lon), psfc,zsfc must be 2-dimensional (lat x lon)");
      return(NhlFATAL);
    }
  }
/*
 * If mkzh is not at least size 3, then this dimension won't be big 
 * enough to contain the cin, lcl, and lfc values.
 */
  if(mkzh < 3) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: The level dimension must have at least 3 elements");
    return(NhlFATAL);
  }

/*
 * Test input dimension sizes.
 */
  if((miy > INT_MAX) || (mjx > INT_MAX) || (mkzh > INT_MAX)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: one or more dimension sizes is greater than INT_MAX");
    return(NhlFATAL);
  }
  imiy = (int) miy;
  imjx = (int) mjx;
  imkzh = (int) mkzh;


/*
 * Calculate size of output array. The output array size depends on
 * the size of p,t,q,z:
 *
 *  - p,t,q,z (nz,time,lev,lat,lon) and psfc,zsfc (nz,time,lat,lon)
 *       output array: (4,nz,time,lat,lon)
 *  - p,t,q,z (time,lev,lat,lon) and psfc,zsfc (time,lat,lon)
 *       output array: (4,time,lat,lon)
 *  - p,t,q,z (lev,lat,lon) and psfc,zsfc (lat,lon)
 *       output array: (4,lat,lon)
 */
  dsizes_cape = (ng_size_t *)calloc(ndims_cape,sizeof(ng_size_t));
  if(dsizes_cape == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: Unable to allocate memory for array dimensionality");
    return(NhlFATAL);
  }

                                    /* 0=cape, 1=cin, 2=lcl, 3=lfc */
  if(ndims_cape == 5) {
    dsizes_cape[0] = 4;    /* To hold the 4 different variables. */
    dsizes_cape[1] = nz;
    dsizes_cape[2] = ntime;
    dsizes_cape[3] = mjx;
    dsizes_cape[4] = miy;
  }
  else if(ndims_cape == 4) {
    dsizes_cape[0] = 4;    /* To hold the 4 different variables. */
    dsizes_cape[1] = ntime;
    dsizes_cape[2] = mjx;
    dsizes_cape[3] = miy;
  }
  else if(ndims_cape == 3) {
    dsizes_cape[0] = 4;    /* To hold the 4 different variables. */
    dsizes_cape[1] = mjx;
    dsizes_cape[2] = miy;
  }

  size_zsfc   = mjx * miy;
  size_cape   = mkzh * size_zsfc;
  mkzh0_index = (mkzh-1) * size_zsfc;    /* Indexes into cin array for   */
  mkzh1_index = (mkzh-2) * size_zsfc;    /* returning cin, lcl, and lfc  */
  mkzh2_index = (mkzh-3) * size_zsfc;    /* respectively. */
  size_left_zsfc = size_zsfc * ntime * nz;
  size_output = 4 * size_left_zsfc;

/* 
 * Allocate space for output and temporary arrays.  Even if the input
 * arrays are already double, go ahead and allocate some space for
 * them b/c we have to copy the values back to 4 different locations.
 *
 * The addition of missing values was added in V6.1.0.
 */
  if(type_p == NCL_double || type_t == NCL_double || type_q == NCL_double ||
     type_z == NCL_double) {
    type_cape     = NCL_double;
    type_obj_cape = nclTypedoubleClass;
    cape          = (double *)calloc(size_output,sizeof(double));
    missing_cape.doubleval = ((NclTypeClass)nclTypedoubleClass)->type_class.default_mis.doubleval;
    cmsg = missing_cape.doubleval;
  }
  else {
    type_cape     = NCL_float;
    type_obj_cape = nclTypefloatClass;
    cape          = (float *)calloc(size_output,sizeof(float));
    missing_cape.floatval = ((NclTypeClass)nclTypefloatClass)->type_class.default_mis.floatval;
    cmsg = (double)missing_cape.floatval;
  }
  tmp_cape = (double *)calloc(size_cape,sizeof(double));
  tmp_cin  = (double *)calloc(size_cape,sizeof(double));
  if(cape == NULL || tmp_cape == NULL || tmp_cin == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: Unable to allocate memory for output arrays");
    return(NhlFATAL);
  }

/*
 * Allocate memory for coercing input arrays to double, if necessary.
 *
 * Force a copy of variable p, because we need to multiply it by 0.01,
 * and we don't want this to propagate back to the NCL script.
 */
  tmp_p_orig = (double *)calloc(size_cape,sizeof(double));
  if(tmp_p_orig == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: Unable to allocate memory for coercing input arrays to double");
    return(NhlFATAL);
  }

  if(type_t != NCL_double) {
    tmp_t_orig = (double *)calloc(size_cape,sizeof(double));
    if(tmp_t_orig == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  if(type_q != NCL_double) {
    tmp_q_orig = (double *)calloc(size_cape,sizeof(double));
    if(tmp_q_orig == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  if(type_z != NCL_double) {
    tmp_z_orig = (double *)calloc(size_cape,sizeof(double));
    if(tmp_z_orig == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  if(type_zsfc != NCL_double) {
    tmp_zsfc = (double *)calloc(size_zsfc,sizeof(double));
    if(tmp_zsfc == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  if(type_psfc != NCL_double) {
    tmp_psfc = (double *)calloc(size_zsfc,sizeof(double));
    if(tmp_psfc == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

/*
 * We need to coerce the pressure array once outside the loop to
 * check if the values are in ascending order.
 *
 * If not, we need to flip the leftmost dimension (p = p(::-1,:,:) in
 * NCL-ese), *and* flip the other 3 input arrays in the same fashion.
 */
  coerce_subset_input_double(p,tmp_p_orig,0,type_p,size_cape,0,NULL,NULL);

  if(tmp_p_orig[0] > tmp_p_orig[(mkzh-1)*size_zsfc] ) {
    flip  = True;
    tmp_p = (double *)calloc(size_cape,sizeof(double));
    tmp_t = (double *)calloc(size_cape,sizeof(double));
    tmp_q = (double *)calloc(size_cape,sizeof(double));
    tmp_z = (double *)calloc(size_cape,sizeof(double));

    if(tmp_p == NULL || tmp_t == NULL || tmp_q == NULL || tmp_z == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: Unable to allocate memory for flipping arrays");
      return(NhlFATAL);
    }
  }
  else  {
    flip  = False;
    tmp_p = tmp_p_orig;
    tmp_t = tmp_t_orig;
    tmp_q = tmp_q_orig;
    tmp_z = tmp_z_orig;
  }
/*
 * Get path to psadilookup.dat file required by this routine. 
 */
  psa_file = get_psa_file();

  /* Allocate space for errmsg*/
  errmsg = (char *) calloc(ERRLEN, sizeof(char))

/*
 * Loop through time,nz and call the Fortran routine.
 */ 
  index_cape        = index_zsfc = 0;
  index_output_cape = 0;
  index_output_cin  = size_left_zsfc;
  index_output_lcl  = 2 * size_left_zsfc;
  index_output_lfc  = 3 * size_left_zsfc;

  for(i = 0; i < ntime*nz; i++) {
/*
 * Coerce subset of input arrays to double if necessary.
 */
    if(i > 0) {
      coerce_subset_input_double(p,tmp_p_orig,index_cape,type_p,
                                 size_cape,0,NULL,NULL);
    }
/*
 * Multiple pressure values by 0.01 to convert from Pa to hPa.
 * The assumption is that pressure values come in as Pa.
 */
    convert_to_hPa(tmp_p_orig,size_cape);

    if(type_t != NCL_double) {
      coerce_subset_input_double(t,tmp_t_orig,index_cape,type_t,
                                 size_cape,0,NULL,NULL);
    }
    else {
/*
 * Point tmp_t_orig to appropriate location in t.
 */
      tmp_t_orig = &((double*)t)[index_cape];
    }
    if(type_q != NCL_double) {
      coerce_subset_input_double(q,tmp_q_orig,index_cape,type_q,
                                 size_cape,0,NULL,NULL);
    }
    else {
/*
 * Point tmp_q_orig to appropriate location in q.
 */
      tmp_q_orig = &((double*)q)[index_cape];
    }
    if(type_z != NCL_double) {
      coerce_subset_input_double(z,tmp_z_orig,index_cape,type_z,
                                 size_cape,0,NULL,NULL);
    }
    else {
/*
 * Point tmp_z_orig to appropriate location in z.
 */
      tmp_z_orig = &((double*)z)[index_cape];
    }

    if(type_psfc != NCL_double) {
      coerce_subset_input_double(psfc,tmp_psfc,index_zsfc,type_psfc,
                                 size_zsfc,0,NULL,NULL);
    }
    else {
/*
 * Point tmp_psfc to appropriate location in psfc.
 */
      tmp_psfc = &((double*)psfc)[index_zsfc];
    }
    if(type_zsfc != NCL_double) {
      coerce_subset_input_double(zsfc,tmp_zsfc,index_zsfc,type_zsfc,
                                 size_zsfc,0,NULL,NULL);
    }
    else {
/*
 * Point tmp_zsfc to appropriate location in zsfc.
 */
      tmp_zsfc = &((double*)zsfc)[index_zsfc];
    }
    
/*
 * If the pressure values need to be flipped, we also need to flip
 * the z, q, and t values in the same fashion.
 */
    if(flip) {
      flip_it(tmp_p_orig,tmp_p,mkzh,size_zsfc);
      flip_it(tmp_t_orig,tmp_t,mkzh,size_zsfc);
      flip_it(tmp_q_orig,tmp_q,mkzh,size_zsfc);
      flip_it(tmp_z_orig,tmp_z,mkzh,size_zsfc);
    }
    else {
      tmp_p = tmp_p_orig;
      tmp_t = tmp_t_orig;
      tmp_q = tmp_q_orig;
      tmp_z = tmp_z_orig;
    }

    errstat = 0;
    errmsg = "";
/*
 * Call Fortran routine.
 */
    NGCALLF(dcapecalc3d,DCAPECALC3D)(tmp_p, tmp_t, tmp_q, tmp_z, tmp_zsfc,
                                     tmp_psfc, tmp_cape, tmp_cin, &cmsg,
                                     &imiy, &imjx, &imkzh, &i3dflag, &iter,
                                     psa_file,&errstat,errmsg,strlen(psa_file),
									 ERRLEN);

    /* Terminate if there was an error */
    if (errstat != 0) {
    	fprintf(stderr, errmsg);
    	exit(errstat);
    }
/*
 * Even if we flipped arrays before going into the Fortran routine, do
 * NOT flip them on the output.
 *
 * Copy the values back out to the correct places in the "cape" array.
 *
 * This is a bit whacky, because the Fortran code is doing something
 * fancy to save memory. The "tmp_cin" array contains the cin values in
 * the last mkzh section, the lcl values in the 2nd-to-last mkzh
 * section, and the lfc values in the 3rd-to-last mkzh section.
 *
 * The "tmp_cape" array contains its values in the last mkzh section
 * of the tmp_cape array.
 */
    coerce_output_float_or_double(cape,&tmp_cape[mkzh0_index],type_cape,
                                  size_zsfc,index_output_cape);
    coerce_output_float_or_double(cape,&tmp_cin[mkzh0_index],type_cape,
                                  size_zsfc,index_output_cin);
    coerce_output_float_or_double(cape,&tmp_cin[mkzh1_index],type_cape,
                                  size_zsfc,index_output_lcl);
    coerce_output_float_or_double(cape,&tmp_cin[mkzh2_index],type_cape,
                                  size_zsfc,index_output_lfc);
/*
 * Implement the pointers into the arrays.
 */
    index_cape += size_cape;
    index_zsfc += size_zsfc;
    index_output_cape += size_zsfc;
    index_output_cin  += size_zsfc;
    index_output_lcl  += size_zsfc;
    index_output_lfc  += size_zsfc;
  }
/*
 * Free memory.
 */
  NclFree(tmp_p_orig);
  if(type_t != NCL_double) NclFree(tmp_t_orig);
  if(type_q != NCL_double) NclFree(tmp_q_orig);
  if(type_z != NCL_double) NclFree(tmp_z_orig);
  if(type_zsfc != NCL_double) NclFree(tmp_zsfc);
  if(type_psfc != NCL_double) NclFree(tmp_psfc);
  NclFree(tmp_cape);
  NclFree(tmp_cin);
  if(flip) {
    NclFree(tmp_p);
    NclFree(tmp_t);
    NclFree(tmp_q);
    NclFree(tmp_z);
  }
  NclFree(psa_file);
  NclFree(errmsg);

/*
 * Get dimension info to see if we have named dimensions.
 * This will be used for return variable.
 */
  dim_info_t = get_wrf_dim_info(1,7,ndims_t,dsizes_t);
  if(dim_info_t != NULL) {
    dim_info = malloc(sizeof(NclDimRec)*ndims_cape);
    if(dim_info == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_cape_2d: Unable to allocate memory for holding dimension information");
      return(NhlFATAL);
    }
    for(i = 0; i < ndims_cape; i++ ) {
      dim_info[i].dim_num  = i;
      dim_info[i].dim_size = dsizes_cape[i];
    }
    dim_info[0].dim_quark = NrmStringToQuark("mcape_mcin_lcl_lfc");
    for(i = 0; i < ndims_t-3; i++) {
      dim_info[i+1].dim_quark = dim_info_t[i].dim_quark;
    }
    dim_info[ndims_cape-2].dim_quark = dim_info_t[ndims_t-2].dim_quark;
    dim_info[ndims_cape-1].dim_quark = dim_info_t[ndims_t-1].dim_quark;

  }


/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)cape,
                            &missing_cape,
                            ndims_cape,
                            dsizes_cape,
                            TEMPORARY,
                            NULL,
                            type_obj_cape
                            );
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          -1,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );

  NclFree(dsizes_cape);
  if(dim_info   != NULL) NclFree(dim_info);
  NclFree(dim_info_t);

/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);
}


/*
 * Retrieve the dimension name info of a particular
 * input argument to a WRF NCL function. If there are
 * no named dimensions, *and* you have at least a 2D
 * array, then set the last two dimension names to
 * "south_north" x "west_east".
 */
NclDimRec *get_wrf_dim_info(int arg_num,int num_args,int ndims_arg,ng_size_t *dsizes_arg)
{
  NclDimRec *dim_info;
  int i, is_named;

  /* this is now separately malloced */
  dim_info = get_dim_info(arg_num,num_args);

  is_named = 0;
  if(ndims_arg >= 2) {
    if(dim_info != NULL) {
/*
 * Check if we actually have any named dimensions.
 */
      i = 0;
      while(i < ndims_arg && !is_named ) {
        if(dim_info[i++].dim_quark != -1) is_named = 1;
      }
    }
    if(!is_named) {
/*
 * If we are here, then we know we have no named dimensions,
 * and hence need to create some.
 */
      if(dim_info == NULL) {
        dim_info = malloc(sizeof(NclDimRec)*ndims_arg);
        if(dim_info == NULL) {
          NhlPError(NhlWARNING,NhlEUNKNOWN,"wrf_get_dim_info: Unable to allocate memory for setting dimension names");
          return(NULL);
        }
      }
      for(i = 0; i < ndims_arg; i++ ) {
        dim_info[i].dim_num   = i;
        dim_info[i].dim_quark = -1;
        dim_info[i].dim_size  = dsizes_arg[i];
      }
      dim_info[ndims_arg-2].dim_quark = NrmStringToQuark("south_north");
      dim_info[ndims_arg-1].dim_quark = NrmStringToQuark("west_east");
    }
  }
  return(dim_info);
}

NhlErrorTypes wrf_eth_W( void )
{
/*
 * Input variables
 */
/*
 * Argument # 0
 */
  void *qv;
  double *tmp_qv = NULL;
  int ndims_qv;
  ng_size_t dsizes_qv[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_qv;

/*
 * Argument # 1
 */
  void *t;
  double *tmp_t = NULL;
  int ndims_t;
  ng_size_t dsizes_t[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_t;

/*
 * Argument # 2
 */
  void *p;
  double *tmp_p = NULL;
  int ndims_p;
  ng_size_t dsizes_p[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_p;

/*
 * Return variable
 */
  void *eth;
  double *tmp_eth = NULL;
  NclBasicDataTypes type_eth;
  NclObjClass type_obj_eth;
  NclQuark *description, *units;
  char *cdescription, *cunits;

/*
 * Variables for returning the output array with dimension names attached.
 */
  int att_id;
  ng_size_t dsizes[1];
  NclMultiDValData return_md, att_md;
  NclVar tmp_var;
  NclStackEntry return_data;

/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info;

/*
 * Various
 */
  int btdim, sndim, wedim, nbtsnwe;
  ng_size_t index_eth, i, size_leftmost, size_eth;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
/*
 * Get argument # 0
 */
  qv = (void*)NclGetArgValue(
           0,
           3,
           &ndims_qv,
           dsizes_qv,
           NULL,
           NULL,
           &type_qv,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_qv < 3) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_eth: The qv array must have at least 3 dimensions");
    return(NhlFATAL);
  }
  btdim = dsizes_qv[ndims_qv-3];
  sndim = dsizes_qv[ndims_qv-2];
  wedim = dsizes_qv[ndims_qv-1];
  nbtsnwe = btdim * sndim * wedim;

/*
 * Get argument # 1
 */
  t = (void*)NclGetArgValue(
           1,
           3,
           &ndims_t,
           dsizes_t,
           NULL,
           NULL,
           &type_t,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_t != ndims_qv) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_eth: The qv and t arrays must have the same number of dimensions");
    return(NhlFATAL);
  }

/*
 * Get argument # 2
 */
  p = (void*)NclGetArgValue(
           2,
           3,
           &ndims_p,
           dsizes_p,
           NULL,
           NULL,
           &type_p,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_p != ndims_qv) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_eth: The p and t arrays must have the same number of dimensions");
    return(NhlFATAL);
  }

  for(i = 0; i < ndims_qv; i++) {
    if(dsizes_t[i] != dsizes_qv[i] || dsizes_p[i] != dsizes_qv[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_eth: The qv, t, and p arrays must have the same dimension sizes");
      return(NhlFATAL);
    }
  }

/*
 * Calculate size of leftmost dimensions.
 */
  size_leftmost = 1;
  for(i = 0; i < ndims_qv-3; i++) size_leftmost *= dsizes_qv[i];

/*
 * The output type defaults to float, unless any input arrays are double.
 */
  type_eth     = NCL_float;
  type_obj_eth = nclTypefloatClass;

/* 
 * Allocate space for coercing input arrays.  If any of the input
 * is already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 *
 * Allocate space for tmp_qv.
 */
  if(type_qv != NCL_double) {
    tmp_qv = (double *)calloc(nbtsnwe,sizeof(double));
    if(tmp_qv == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_eth: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_eth     = NCL_double;
    type_obj_eth = nclTypedoubleClass;
  }
/*
 * Allocate space for tmp_t.
 */
  if(type_t != NCL_double) {
    tmp_t = (double *)calloc(nbtsnwe,sizeof(double));
    if(tmp_t == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_eth: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_eth     = NCL_double;
    type_obj_eth = nclTypedoubleClass;
  }
/*
 * Allocate space for tmp_p.
 */
  if(type_p != NCL_double) {
    tmp_p = (double *)calloc(nbtsnwe,sizeof(double));
    if(tmp_p == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_eth: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_eth     = NCL_double;
    type_obj_eth = nclTypedoubleClass;
  }

/*
 * Calculate size of output array and allocate space for it.
 */
  size_eth = size_leftmost * nbtsnwe;

  if(type_eth != NCL_double) {
    eth = (void *)calloc(size_eth, sizeof(float));
    tmp_eth = (double *)calloc(nbtsnwe,sizeof(double));
    if(tmp_eth == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_eth: Unable to allocate memory for temporary output array");
      return(NhlFATAL);
    }
  }
  else {
    eth = (void *)calloc(size_eth, sizeof(double));
  }
  if(eth == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_eth: Unable to allocate memory for output array");
    return(NhlFATAL);
  }

/*
 * Loop across leftmost dimensions and call the Fortran routine for each
 * subsection of the input arrays.
 */
  index_eth = 0;

  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of qv (tmp_qv) to double if necessary.
 */
    if(type_qv != NCL_double) {
      coerce_subset_input_double(qv,tmp_qv,index_eth,type_qv,nbtsnwe,
                                 0,NULL,NULL);
    }
    else {
      tmp_qv = &((double*)qv)[index_eth];
    }

/*
 * Coerce subsection of t (tmp_t) to double if necessary.
 */
    if(type_t != NCL_double) {
      coerce_subset_input_double(t,tmp_t,index_eth,type_t,nbtsnwe,
                                 0,NULL,NULL);
    }
    else {
      tmp_t = &((double*)t)[index_eth];
    }

/*
 * Coerce subsection of p (tmp_p) to double if necessary.
 */
    if(type_p != NCL_double) {
      coerce_subset_input_double(p,tmp_p,index_eth,type_p,nbtsnwe,
                                 0,NULL,NULL);
    }
    else {
      tmp_p = &((double*)p)[index_eth];
    }

/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_eth == NCL_double) tmp_eth = &((double*)eth)[index_eth];

/*
 * Call the Fortran routine.
 */
    NGCALLF(deqthecalc,DEQTHECALC)(tmp_qv, tmp_t, tmp_p, tmp_eth, 
                                   &wedim, &sndim, &btdim);

/*
 * Coerce output back to float if necessary.
 */
    if(type_eth == NCL_float) {
      coerce_output_float_only(eth,tmp_eth,nbtsnwe,index_eth);
    }
    index_eth += nbtsnwe;
  }

/*
 * Free unneeded memory.
 */
  if(type_qv  != NCL_double) NclFree(tmp_qv);
  if(type_t   != NCL_double) NclFree(tmp_t);
  if(type_p   != NCL_double) NclFree(tmp_p);
  if(type_eth != NCL_double) NclFree(tmp_eth);

/*
 * Retrieve dimension names from the "t" variable, if any.
 * These dimension names will later be attached to the output variable.
 */
  dim_info = get_wrf_dim_info(1,3,ndims_t,dsizes_t);

/*
 * Set up return value.
 */
/*
 * Return value back to NCL script.
 */

  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)eth,
                            NULL,
                            ndims_t,
                            dsizes_t,
                            TEMPORARY,
                            NULL,
                            type_obj_eth
                            );
/*
 * Set up some attributes ("description" and "units") to return.
 */
  cdescription = (char *)calloc(33,sizeof(char));
  strcpy(cdescription,"Equivalent Potential Temperature");
  description  = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *description = NrmStringToQuark(cdescription);
  free(cdescription);

  cunits       = (char *)calloc(2,sizeof(char));
  strcpy(cunits,"K");
  units        = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *units       = NrmStringToQuark(cunits);
  free(cunits);

/*
 * Set up attributes to return.
 */
  att_id = _NclAttCreate(NULL,NULL,Ncl_Att,0,NULL);

  dsizes[0] = 1;
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)description,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "description",
             att_md,
             NULL
             );
    
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)units,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "units",
             att_md,
             NULL
             );
    
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          att_id,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );

  NclFree(dim_info);
/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);

}

/*
 * The wrf_wps_xxx_int functions that follow are for reading
 * WRF WPS intermediate files.  You can use wrf_wps_read_int
 * as an "all-in-one" function that creates one 3D variable
 * with all the data and attributes attached, or you can 
 * use the three individual functions, wrf_wps_open_int, 
 * wrf_wps_rdhead_int, and wrf_wps_rddata_int that allows you
 * to read the data in one 2D slab at a time.
 */

/*
 * This function simply opens the WRF/WPS intermediate file
 * and returns a status.
 */
NhlErrorTypes wrf_wps_open_int_W( void )
{

/*
 * Argument # 0
 */
  NrmQuark *filename;
  char *cfilename;
/*
 * Return variable
 */
  int *istatus;

/*
 * Various
 */
  int ret, ndims;
  ng_size_t dsizes[1];

/*
 * Get argument # 0
 */
  filename = (NrmQuark *)NclGetArgValue(
           0,
           1,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);
/*
 * Convert to character string.
 */
  cfilename = NrmQuarkToString(*filename);

/* 
 * Allocate space for output array.
 */
  ndims     = 1;
  dsizes[0] = 1;


/*  
 * Allocate return integer
 */
  istatus = (int*)calloc(1, sizeof(int));
/*
 * Call the Fortran routine.
 */
  NGCALLF(plotfmt_open,PLOTFMT_open)(cfilename, istatus,
                                     strlen(cfilename));
  if(*istatus != 0) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_wps_open_int: The input file '%s' could not be opened.\nCheck that it exists and is spelled correctly",cfilename);
    return(NhlFATAL);
  }
/*
 * Return value back to NCL script.
 */
  ret = NclReturnValue(istatus,ndims,dsizes,NULL,NCL_int,0);
  return(ret);
}

/*
 * This function simply closes a WRF/WPS intermediate file
 * that was opened with wrf_wps_open/rddata/rhead_int.
 * Note that nothing is currently done with the istatus
 * variable. This might change if we decide that 
 * istatus can contain a "unit" attribute that indicates
 * which Fortran unit was used to open the file.
 */
NhlErrorTypes wrf_wps_close_int_W( void )
{
  int *istatus;

  istatus = (int *)NclGetArgValue(
           0,
           1,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);

  NGCALLF(plotfmt_close,PLOTFMT_CLOSE)();
  return(NhlNOERROR);

}


/*
 * This function takes an open WRF/WPS intermediate file,
 * and reads header information. This header information
 * then gives you the info needed to read the slab via
 * wrf_wps_rddata_int.
 */
NhlErrorTypes wrf_wps_rdhead_int_W( void )
{

/*
 * Argument # 0
 */
  int *istatus;

/*
 * Argument # 1
 */
  void *head;
  float *rhead;
  NclBasicDataTypes type_head;
  ng_size_t dsizes_head[1];

/*
 * Argument # 2
 */
  NrmQuark *field;
  char *cfield;
  int FIELD_LEN=9;
/*
 * Argument # 3
 */
  NrmQuark *hdate;
  char *chdate;
  int HDATE_LEN=24;

/*
 * Argument # 4
 */
  NrmQuark *units;
  char *cunits;
  int UNITS_LEN=25;

/*
 * Argument # 5
 */
  NrmQuark *mapsc;
  char *cmapsc;
  int MAPSC_LEN=32;

/*
 * Argument # 6
 */
  NrmQuark *desc;
  char *cdesc;
  int DESC_LEN=46;

/*
 * Various
 */
  int i;

/*
 * Get argument #0
 */
  istatus = (int *)NclGetArgValue(
           0,
           7,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);

/*
 * Get arguments # 1-6
 */
  head = (void *)NclGetArgValue(
           1,
           7,
           NULL,
           dsizes_head,
           NULL,
           NULL,
           &type_head,
           DONT_CARE);

  rhead = coerce_input_float(head, type_head, dsizes_head[0], 0, 
                             NULL, NULL);

  field = (NrmQuark *)NclGetArgValue(
           2,
           7,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);
  hdate = (NrmQuark *)NclGetArgValue(
           3,
           7,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);
  units = (NrmQuark *)NclGetArgValue(
           4,
           7,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);
  mapsc = (NrmQuark *)NclGetArgValue(
           5,
           7,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);
  desc = (NrmQuark *)NclGetArgValue(
           6,
           7,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);

/*
 * Allocate space needed for each string
 */
  cfield = (char *)malloc((FIELD_LEN+1)*sizeof(char));
  chdate = (char *)malloc((HDATE_LEN+1)*sizeof(char));
  cunits = (char *)malloc((UNITS_LEN+1)*sizeof(char));
  cmapsc = (char *)malloc((MAPSC_LEN+1)*sizeof(char));
  cdesc  = (char *)malloc((DESC_LEN+1)*sizeof(char));
/*
 * Call the Fortran routine.
 */
  NGCALLF(plotfmt_rdhead,PLOTFMT_RDHEAD)(istatus,rhead,cfield,chdate,
                                         cunits,cmapsc,cdesc,
                                         FIELD_LEN,HDATE_LEN,
                                         UNITS_LEN,MAPSC_LEN,
                                         DESC_LEN);

/*
 * Strip off potential white space at end of each string.
 */
  i = FIELD_LEN-1;
  while(i >=0 && (cfield[i] == ' ' || cfield[i] == '\t')) i--;
  cfield[i+1] = '\0';

  i = HDATE_LEN-1;
  while(i >= 0 && (chdate[i] == ' ' || chdate[i] == '\t')) i--;
  chdate[i+1] = '\0';

  i = UNITS_LEN-1;
  while( i >= 0 && (cunits[i] == ' ' || cunits[i] == '\t')) i--;
  cunits[i+1] = '\0';

  i = MAPSC_LEN-1;
  while( i >= 0 && (cmapsc[i] == ' ' || cmapsc[i] == '\t')) i--;
  cmapsc[i+1] = '\0';

  i = DESC_LEN-1;
  while( i >= 0 && (cdesc[i] == ' ' || cdesc[i] == '\t')) i--;
  cdesc[i+1] = '\0';

  *field = NrmStringToQuark(cfield);
  *hdate = NrmStringToQuark(chdate);
  *units = NrmStringToQuark(cunits);
  *mapsc = NrmStringToQuark(cmapsc);
  *desc  = NrmStringToQuark(cdesc);

  return(NhlNOERROR);
}

/*
 * This function takes an open WRF/WPS intermediate file,
 * and reads a 2D slab. 
 */
NhlErrorTypes wrf_wps_rddata_int_W( void )
{

/*
 * Argument # 0
 */
  int *istatus;

/*
 * Arguments #1-2
 */
  void *tmp_nx, *tmp_ny;
  ng_size_t *nx, *ny;
  NclBasicDataTypes type_nx, type_ny;
  int inx, iny;

/*
 * Return
 */
  float *slab;
  ng_size_t dsizes[2];
  int ret;

/*
 * Get argument #0
 */
  istatus = (int *)NclGetArgValue(
           0,
           3,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);

/*
 * Get arguments #1-2
 */
  tmp_nx = (void *)NclGetArgValue(
           1,
           3,
           NULL,
           NULL,
           NULL,
           NULL,
           &type_nx,
           DONT_CARE);

  tmp_ny = (void *)NclGetArgValue(
           2,
           3,
           NULL,
           NULL,
           NULL,
           NULL,
           &type_ny,
           DONT_CARE);

/*
 * Convert the input dimensions to ng_size_t.
 */
  nx = get_dimensions(tmp_nx,1,type_nx,"wrf_wps_rddata_int");
  ny = get_dimensions(tmp_ny,1,type_ny,"wrf_wps_rddata_int");
  if(nx == NULL || ny == NULL) 
    return(NhlFATAL);

  if((*nx > INT_MAX) || (*ny > INT_MAX)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_wps_rddata_int: nx and/or ny is greater than INT_MAX");
    return(NhlFATAL);
  }
  inx = (int) *nx;
  iny = (int) *ny;

  slab = (float*)calloc(inx*iny,sizeof(float));
  if(slab == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_wps_rddata_int: Unable to allocate memory for output array");
    return(NhlFATAL);
  }

/*
 * Call the Fortran routine.
 */
  NGCALLF(plotfmt_rddata,PLOTFMT_RDDATA)(istatus,&inx,&iny,slab);

/*
 * Return value back to NCL script.
 */
  dsizes[0] = *ny;
  dsizes[1] = *nx;
  ret = NclReturnValue(slab,2,dsizes,NULL,NCL_float,0);
  return(ret);
}


/*
 * This function is a "3-in-1" function that does the
 * work of wrf_wps_open_int, wrf_wps_rdhead_int, and 
 * wrf_wps_rddata_int. It returns a 3D float array 
 * that contains all the necessary data and attributes.
 * 
 */
NhlErrorTypes wrf_wps_read_int_W( void )
{

/*
 * Argument # 0
 */
  NrmQuark *filename;
  char *cfilename;

/*
 * Return values.  The slab will be returned, along with a bunch
 * of attributes depending on the projection.
 */
  int NHEAD=14, FIELD_LEN=9, HDATE_LEN=24, UNITS_LEN=25;
  int MAPSC_LEN=32, DESCR_LEN=46;
  float *slab, *rhead, *slab_s;
  float rhead_s[NHEAD];
  ng_size_t dsizes_slab[3], dsizes_rhead[2], dsizes_field[1], dsizes_hdate[1];
  ng_size_t dsizes_units[1], dsizes_mapsc[1], dsizes_descr[1];
  NrmQuark *field, *hdate, *units, *mapsc, *descr;
  char *cfield, *chdate, *cunits, *cmapsc, *cdescr;
  NclScalar missing_slab;
/*
 * Various 
 */
  int i, j, n, istatus, nx, ny, nxny, max_nx, max_ny, nfields;
  int index_slab, index_rhead; 

/*
 * Attribute variables
 */
  int att_id;
  NclMultiDValData att_md, return_md;
  NclVar tmp_var;
  NclStackEntry return_data;


/*
 * Get argument #0
 */
  filename = (NrmQuark *)NclGetArgValue(
           0,
           1,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);

/*
 * Convert to character string.
 */
  cfilename = NrmQuarkToString(*filename);

/*
 * Call the Fortran routine to open the file
 */
  NGCALLF(plotfmt_open,PLOTFMT_open)(cfilename, &istatus,
                                     strlen(cfilename));
  if(istatus != 0) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_wps_read_int: The input file '%s' could not be opened.\nCheck that it exists and is spelled correctly",cfilename);
    return(NhlFATAL);
  }

  cfield = (char *)malloc((FIELD_LEN+1)*sizeof(char));
  chdate = (char *)malloc((HDATE_LEN+1)*sizeof(char));
  cunits = (char *)malloc((UNITS_LEN+1)*sizeof(char));
  cmapsc = (char *)malloc((MAPSC_LEN+1)*sizeof(char));
  cdescr = (char *)malloc((DESCR_LEN+1)*sizeof(char));

/*
 *  Read each field so we can count how many there are.
 */
  nfields = 0;
  while (istatus == 0) {
    /* Read the header */

    NGCALLF(plotfmt_rdhead,PLOTFMT_RDHEAD)(&istatus,&rhead_s[0],cfield,
                                           chdate,cunits,cmapsc,cdescr,
                                           FIELD_LEN,HDATE_LEN,
                                           UNITS_LEN,MAPSC_LEN,
                                           DESCR_LEN);

    if(istatus == 0) {
      nx = (int)rhead_s[3];
      ny = (int)rhead_s[4];
      slab_s = (float*)calloc(nx*ny,sizeof(float));
      if(slab_s == NULL) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_wps_read_int: Unable to allocate memory for output array");
        return(NhlFATAL);
      }
      if(nfields == 0) {
        max_nx = nx;
        max_ny = ny;
      }
      else {
        max_nx = max(max_nx,nx);
        max_ny = max(max_ny,ny);
      }
    /* Read the data */
      NGCALLF(plotfmt_rddata,PLOTFMT_RDDATA)(&istatus,&nx,&ny,slab_s);
      NclFree(slab_s);
      if(istatus == 0) nfields++; 
    }
  }
  nxny = max_nx * max_ny;

  /* Allocate the return arrays */
  rhead  = (float*)calloc(nfields*NHEAD,sizeof(float));
  slab   = (float*)calloc(nfields*nxny,sizeof(float));
  field  = (NclQuark*)NclMalloc(nfields*sizeof(NclQuark));
  hdate  = (NclQuark*)NclMalloc(nfields*sizeof(NclQuark));
  units  = (NclQuark*)NclMalloc(nfields*sizeof(NclQuark));
  mapsc  = (NclQuark*)NclMalloc(nfields*sizeof(NclQuark));
  descr  = (NclQuark*)NclMalloc(nfields*sizeof(NclQuark));

  if(rhead == NULL || slab == NULL || field == NULL || hdate == NULL ||
     units == NULL || mapsc == NULL || descr == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_wps_read_int: Unable to allocate memory for output array and attributes");
    return(NhlFATAL);
  }

  /* Close and reopen file. */
  NGCALLF(plotfmt_close,PLOTFMT_CLOSE)();
  NGCALLF(plotfmt_open,PLOTFMT_open)(cfilename, &istatus,
                                     strlen(cfilename));

  missing_slab.floatval = -1e30;
  index_slab = index_rhead = 0;

  /* Now loop through again and fill up the output arrays */

  for(n = 0; n < nfields; n++) {

    /* Read the header into allocated arrays */
    NGCALLF(plotfmt_rdhead,PLOTFMT_RDHEAD)(&istatus,&rhead[index_rhead],
                                           cfield,chdate,cunits,
                                           cmapsc,cdescr,
                                           FIELD_LEN,HDATE_LEN,
                                           UNITS_LEN,MAPSC_LEN,
                                           DESCR_LEN);
/*
 * Strip off potential white space at end of each string.
 */
    i = FIELD_LEN-1;
    while(i >=0 && (cfield[i] == ' ' || 
                    cfield[i] == '\t')) i--;
    cfield[i+1] = '\0';
    
    i = HDATE_LEN-1;
    while(i >= 0 && (chdate[i] == ' ' || 
                     chdate[i] == '\t')) i--;
    chdate[i+1] = '\0';
    
    i = UNITS_LEN-1;
    while( i >= 0 && (cunits[i] == ' ' || 
                      cunits[i] == '\t')) i--;
    cunits[i+1] = '\0';
    
    i = MAPSC_LEN-1;
    while( i >= 0 && (cmapsc[i] == ' ' || 
                      cmapsc[i] == '\t')) i--;
    cmapsc[i+1] = '\0';
    
    i = DESCR_LEN-1;
    while( i >= 0 && (cdescr[i] == ' ' || 
                      cdescr[i] == '\t')) i--;
    cdescr[i+1] = '\0';
    
    field[n] = NrmStringToQuark(cfield);
    hdate[n] = NrmStringToQuark(chdate);
    units[n] = NrmStringToQuark(cunits);
    mapsc[n] = NrmStringToQuark(cmapsc);
    descr[n] = NrmStringToQuark(cdescr);

/*
 * We have to get nx and ny every time in the loop, because they 
 * can be different sizes.
 */
    nx = (int)rhead[index_rhead+3];
    ny = (int)rhead[index_rhead+4];
    slab_s = (float*)calloc(nx*ny,sizeof(float));
    if(slab_s == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_wps_read_int: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
    /* Read data into temporary array */
    NGCALLF(plotfmt_rddata,PLOTFMT_RDDATA)(&istatus,&nx,&ny,&slab_s[0]);

    /* Copy slab subset back out to big slab array */
    for(i = 0; i < ny; i++) {
      for(j = 0; j < nx; j++) {
        slab[index_slab+(i*max_nx)+j] = slab_s[(i*nx)+j];
      }
    /* Fill rest with missing values */
      for(j = nx; j < max_nx; j++) {
        slab[index_slab+(i*max_nx)+j] = missing_slab.floatval;
      }
    }

    /* For next time through loop */
    index_slab  += nxny;
    index_rhead += NHEAD;
  }

  /* Close the file. */
  NGCALLF(plotfmt_close,PLOTFMT_CLOSE)();

  /* Free memory */
  free(slab_s);
  free(cfield);
  free(chdate);
  free(cunits);
  free(cmapsc);
  free(cdescr);

/*
 * Return slab and attributes back to NCL script.
 */
  dsizes_slab[0]  = nfields;
  dsizes_slab[1]  = max_ny;
  dsizes_slab[2]  = max_nx;
  dsizes_rhead[0] = nfields;
  dsizes_rhead[1] = NHEAD;
  dsizes_field[0] = nfields;
  dsizes_hdate[0] = nfields;
  dsizes_units[0] = nfields;
  dsizes_mapsc[0] = nfields;
  dsizes_descr[0] = nfields;

  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)slab,
                            &missing_slab,
                            3,
                            dsizes_slab,
                            TEMPORARY,
                            NULL,
                            nclTypefloatClass
                            );
/*
 * Set up attributes to return.
 */
  att_id = _NclAttCreate(NULL,NULL,Ncl_Att,0,NULL);

  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)field,
                         NULL,
                         1,
                         dsizes_field,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "field",
             att_md,
             NULL
             );
    
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)hdate,
                         NULL,
                         1,
                         dsizes_hdate,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "hdate",
             att_md,
             NULL
             );
    
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)units,
                         NULL,
                         1,
                         dsizes_units,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "units",
             att_md,
             NULL
             );
    
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)mapsc,
                         NULL,
                         1,
                         dsizes_mapsc,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "map_source",
             att_md,
             NULL
             );
    
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)descr,
                         NULL,
                         1,
                         dsizes_descr,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "description",
             att_md,
             NULL
             );
    

  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)rhead,
                         NULL,
                         2,
                         dsizes_rhead,
                         TEMPORARY,
                         NULL,
                         nclTypefloatClass
                         );

  _NclAddAtt(
             att_id,
             "rhead",
             att_md,
             NULL
             );

  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          NULL,
                          att_id,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );

/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);

  return(NhlNOERROR);
}


/*
 * This next section contains wrappers for wrf_bint and wrf_iclw.
 * These wrappers have never been tested or documented, and are
 * not registered in wrapper.c. I've let them here in case we need
 * to resuscitate them.
 */


NhlErrorTypes wrf_wps_read_nml_W( void )
{

/*
 * Argument # 0
 */
  NrmQuark *namelist;
  char *cnamelist;
/*
 * Return variable
 */
  float *pgrids_var;
  int size_output, ndims_output;
  ng_size_t dsizes_output[2];
  NclScalar missing_output;

/*
 * Various
 */
  int NVAR=19, MAX_DOMAINS=21;
  int ret;
/*
 * Get argument # 0
 */
  namelist = (NrmQuark *)NclGetArgValue(
           0,
           1,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);
/*
 * Convert to character string.
 */
  cnamelist = NrmQuarkToString(*namelist);

/* 
 * Allocate space for output array.
 */
  ndims_output     = 2;
  dsizes_output[0] = MAX_DOMAINS; 
  dsizes_output[1] = NVAR;
  size_output      = NVAR*MAX_DOMAINS;
  pgrids_var       = (float*)calloc(size_output, sizeof(float));
  if(pgrids_var == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_wps_read_nml: Unable to allocate memory for output array");
    return(NhlFATAL);
  }
  missing_output.floatval = -999.;

/*
 * Call the Fortran routine.
 */
  NGCALLF(plotgrids_var,PLOTGRIDS_VAR)(cnamelist, pgrids_var, 
                                       strlen(cnamelist));
/*
 * Return value back to NCL script.
 */
  ret = NclReturnValue(pgrids_var,ndims_output,dsizes_output,
                       &missing_output,NCL_float,0);
  return(ret);
}

NhlErrorTypes wrf_bint_W( void )
{
/*
 * Input array variables
 */
  void *data_in, *obsii, *obsjj;
  double *tmp_data_in = NULL;
  double *tmp_obsii = NULL;
  double *tmp_obsjj = NULL;

  int *icrs, *jcrs;
  int ndims_data_in, ndims_obsii, ndims_obsjj;
  ng_size_t dsizes_data_in[NCL_MAX_DIMENSIONS]; 
  ng_size_t dsizes_obsii[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_obsjj[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_data_in, type_obsii, type_obsjj;
  
/*
 * Output variable.
 */
  void *data_out;
  double *tmp_data_out = NULL;
  ng_size_t *dsizes_data_out, size_data_out;
  NclBasicDataTypes type_data_out;
/*
 * Various
 */
  int ret;
  ng_size_t i, nx, ny, nz, nobsicrs, nobsjcrs, size_leftmost;
  ng_size_t nxyz, nobsij, nobsijz, index_data_in, index_data_out, index_nobsij;
  int inx, iny, inz, inobsicrs, inobsjcrs;
/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
  data_in = (void*)NclGetArgValue(
           0,
           5,
           &ndims_data_in,
           dsizes_data_in,
           NULL,
           NULL,
           &type_data_in,
           DONT_CARE);

  obsii = (void*)NclGetArgValue(
           1,
           5,
           &ndims_obsii,
           dsizes_obsii,
           NULL,
           NULL,
           &type_obsii,
           DONT_CARE);

  obsjj = (void*)NclGetArgValue(
           2,
           5,
           &ndims_obsjj,
           dsizes_obsjj,
           NULL,
           NULL,
           &type_obsjj,
           DONT_CARE);

  icrs = (int*)NclGetArgValue(
           3,
           5,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);

  jcrs = (int*)NclGetArgValue(
           4,
           5,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);

/*
 * Error checking.
 */
  if(ndims_data_in < 2 || ndims_obsii < 2 || ndims_obsjj < 2) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_bint: The data_in, obsii, and obsjj arrays must have at least two dimensions");
    return(NhlFATAL);
  }
  if(ndims_obsii != ndims_obsjj) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_bint: The obsii and obsjj arrays must have the same number of dimensions");
    return(NhlFATAL);
  }
  if((ndims_data_in == 2 && ndims_obsii != 2) || 
     (ndims_data_in  > 2 && ndims_data_in != (ndims_obsii+1))) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_bint: The data_in, obsii, and obsjj arrays must all be two-dimensional, or data_in must be greater than two dimensions and have one more dimension than obsii and obsjj");
    return(NhlFATAL);
  }
  for(i = 0; i < ndims_obsii; i++) {
    if(dsizes_obsii[i] != dsizes_obsjj[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_bint: The obsii and obsjj arrays must be the same dimension sizes");
      return(NhlFATAL);
    }
  }

/*
 * If data_in is greater than 3 dimensions, then check that these
 * extra dimensions are all the same length in the three input 
 * arrays.
 *
 * While we're here, calculate the size of the leftmost dimensions.
 */
  size_leftmost = 1;
  if(ndims_data_in > 3) {
    for(i = 0; i < ndims_data_in-3; i++) {
      if(dsizes_data_in[i] != dsizes_obsii[i]) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_bint: The rightmost dimensions of data_in, obsii, obsjj must be the same");
        return(NhlFATAL);
      }
      size_leftmost *= dsizes_data_in[i];
    }
  }

/*
 * Store some dimension sizes and output data array sizes.
 */
  nx = dsizes_data_in[ndims_data_in-1];
  ny = dsizes_data_in[ndims_data_in-2];
  if(ndims_data_in > 2) {
    nz = dsizes_data_in[ndims_data_in-3];
  }
  else {
    nz = 1;
  }
  nobsicrs = dsizes_obsii[ndims_obsii-1];
  nobsjcrs = dsizes_obsii[ndims_obsii-2];
  nxyz     = nx * ny * nz;
  nobsij   = nobsicrs * nobsjcrs;
  nobsijz  = nobsij * nz;

  size_data_out = size_leftmost * nobsijz;

/*
 * Test input dimension sizes.
 */
  if((nx > INT_MAX) || (ny > INT_MAX) || (nobsicrs > INT_MAX) || 
     (nobsjcrs > INT_MAX) || (nz > INT_MAX)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_bint: one or more dimension sizes is greater than INT_MAX");
    return(NhlFATAL);
    }
    inx = (int) nx;
    iny = (int) ny;
    inz = (int) nz;
    inobsicrs = (int) nobsicrs;
    inobsjcrs = (int) nobsjcrs;

/* 
 * Allocate space for coercing input arrays.  If the input data_in, obsii,
 * or obsjj are already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 *
 * The output type defaults to float, unless any of the input arrays
 * are double.
 */
  type_data_out = NCL_float;
  if(type_data_in != NCL_double) {
    tmp_data_in = (double *)calloc(nxyz,sizeof(double));
    if(tmp_data_in == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_bint: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_data_out = NCL_double;
  }

  if(type_obsii != NCL_double) {
    tmp_obsii = (double *)calloc(nobsij,sizeof(double));
    if(tmp_obsii == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_bint: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_data_out = NCL_double;
  }

  if(type_obsjj != NCL_double) {
    tmp_obsjj = (double *)calloc(nobsij,sizeof(double));
    if(tmp_obsjj == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_bint: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_data_out = NCL_double;
  }

/*
 * Allocate space for output array.
 */ 
  if(type_data_out == NCL_double) {
    data_out = (double *)calloc(size_data_out,sizeof(double));
    if(data_out == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_bint: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }
  else {
    data_out     = (float *)calloc(size_data_out,sizeof(float));
    tmp_data_out = (double *)calloc(nobsijz,sizeof(double));
    if(tmp_data_out == NULL || data_out == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_bint: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }

/*
 * Create dimension sizes for output array.
 */
  dsizes_data_out = (ng_size_t*)calloc(ndims_data_in,sizeof(ng_size_t));  
  if( dsizes_data_out == NULL ) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_bint: Unable to allocate memory for holding dimension sizes");
    return(NhlFATAL);
  }
  for(i = 1; i < ndims_data_in-2; i++) dsizes_data_out[i] = dsizes_data_in[i];
  dsizes_data_out[ndims_data_in-2] = nobsjcrs;
  dsizes_data_out[ndims_data_in-1] = nobsicrs;
/*
 * Loop across leftmost dimensions and call the Fortran routine for each
 * one-dimensional subsection.
 */
  index_data_in = index_data_out = index_nobsij = 0;
  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of data_in (tmp_data_in) to double if necessary.
 */
    if(type_data_in != NCL_double) {
      coerce_subset_input_double(data_in,tmp_data_in,index_data_in,
                                 type_data_in,nxyz,0,NULL,NULL);
    }
    else {
      tmp_data_in = &((double*)data_in)[index_data_in];
    }
/*
 * Coerce subsection of obsii (tmp_obsii) to double if ncessary.
 */
    if(type_obsii != NCL_double) {
      coerce_subset_input_double(obsii,tmp_obsii,index_nobsij,type_obsii,
                                 nobsij,0,NULL,NULL);
    }
    else {
      tmp_obsii = &((double*)obsii)[index_nobsij];
    }

/*
 * Coerce subsection of obsjj (tmp_obsjj) to double if ncessary.
 */
    if(type_obsjj != NCL_double) {
      coerce_subset_input_double(obsjj,tmp_obsjj,index_nobsij,type_obsjj,
                                 nobsij,0,NULL,NULL);
    }
    else {
      tmp_obsjj = &((double*)obsjj)[index_nobsij];
    }

/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_data_out == NCL_double) {
      tmp_data_out = &((double*)data_out)[index_data_out];
    }
/*
 * Call Fortran routine.
 */
    NGCALLF(dbint3d,DBINT3D)(tmp_data_out,tmp_obsii,tmp_obsjj,tmp_data_in,
                             &inx,&iny,&inz,&inobsicrs,&inobsjcrs,icrs,jcrs);
/*
 * Coerce output back to float if necessary.
 */
    if(type_data_out == NCL_float) {
      coerce_output_float_only(data_out,tmp_data_out,nobsijz,index_data_out);
    }

/*
 * Increment indices.
 */
    index_data_in  += nxyz;
    index_data_out += nobsijz;
    index_nobsij   += nobsij;
  }
/*
 * Free up memory.
 */
  if(type_data_in  != NCL_double) NclFree(tmp_data_in);
  if(type_obsii    != NCL_double) NclFree(tmp_obsii);
  if(type_obsjj    != NCL_double) NclFree(tmp_obsjj);
  if(type_data_out != NCL_double) NclFree(tmp_data_out);

  ret = NclReturnValue(data_out,ndims_data_in,dsizes_data_out,NULL,
                        type_data_out,0);
  NclFree(dsizes_data_out);
  return(ret);
}


NhlErrorTypes wrf_iclw_W( void )
{

/*
 * Input variables
 */
/*
 * Argument # 0
 */
  void *p;
  double *tmp_p = NULL;
  int ndims_p;
  ng_size_t dsizes_p[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_p;

/*
 * Argument # 1
 */
  void *qc;
  double *tmp_qc = NULL;
  int ndims_qc;
  ng_size_t dsizes_qc[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_qc;

/*
 * Return variable
 */
  void *iclw;
  NclQuark *description, *units;
  char *cdescription, *cunits;
  double *tmp_iclw = NULL;
  int ndims_iclw;
  ng_size_t *dsizes_iclw;
  NclBasicDataTypes type_iclw;
  NclObjClass type_obj_iclw;
/*
 * Variables for returning the output array with attributes attached.
 */
  int att_id;
  ng_size_t dsizes[1];
  NclMultiDValData att_md, return_md;
  NclVar tmp_var;
  NclStackEntry return_data;

/*
 * Various
 */
  ng_size_t nz, ny, nx, nznynx, nynx;
  ng_size_t index_p, index_iclw;
  ng_size_t i, ndims_leftmost, size_leftmost, size_output;
  int inx, iny, inz;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
/*
 * Get argument # 0
 */
  p = (void*)NclGetArgValue(
           0,
           2,
           &ndims_p,
           dsizes_p,
           NULL,
           NULL,
           &type_p,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_p < 3) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_iclw: The p array must have at least 3 dimensions");
    return(NhlFATAL);
  }
  nz = dsizes_p[ndims_p-3];
  ny = dsizes_p[ndims_p-2];
  nx = dsizes_p[ndims_p-1];
  nynx   = ny * nx;
  nznynx = nz * nynx;

/*
 * Test dimension sizes.
 */
  if((nx > INT_MAX) || (ny > INT_MAX) || (nz > INT_MAX)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_iclw: nx, ny and/or is greater than INT_MAX");
    return(NhlFATAL);
  }
  inx = (int) nx;
  iny = (int) ny;
  inz = (int) nz;

/*
 * Get argument # 1
 */
  qc = (void*)NclGetArgValue(
           1,
           2,
           &ndims_qc,
           dsizes_qc,
           NULL,
           NULL,
           &type_qc,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_qc < 3) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_iclw: The qc array must have at least 3 dimensions");
    return(NhlFATAL);
  }

  if(dsizes_qc[ndims_qc-3] != nz || 
     dsizes_qc[ndims_qc-2] != ny || 
     dsizes_qc[ndims_qc-1] != nx) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_iclw: The rightmost dimensions of qc must be nz x ny x nx");
    return(NhlFATAL);
  }

/*
 * Calculate size of leftmost dimensions.
 */
  size_leftmost  = 1;
  ndims_leftmost = ndims_p-3;
  for(i = 0; i < ndims_leftmost; i++) {
    if(dsizes_qc[i] != dsizes_p[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_iclw: The leftmost dimensions of p and qc must be the same");
      return(NhlFATAL);
    }
    size_leftmost *= dsizes_p[i];
  }

/*
 * The output type defaults to float, unless either input array is double.
 */
  type_iclw     = NCL_float;
  type_obj_iclw = nclTypefloatClass;

/* 
 * Allocate space for coercing input arrays.  If any of the input
 * is already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 */
/*
 * Allocate space for tmp_p.
 */
  if(type_p != NCL_double) {
    tmp_p = (double *)calloc(nznynx,sizeof(double));
    if(tmp_p == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_iclw: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_iclw     = NCL_double;
    type_obj_iclw = nclTypedoubleClass;
  }
/*
 * Allocate space for tmp_qc.
 */
  if(type_qc != NCL_double) {
    tmp_qc = (double *)calloc(nznynx,sizeof(double));
    if(tmp_qc == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_iclw: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_iclw     = NCL_double;
    type_obj_iclw = nclTypedoubleClass;
  }

/*
 * Calculate size of output array.
 */
  size_output = size_leftmost * nynx;

/* 
 * Allocate space for output array.
 */
  if(type_iclw != NCL_double) {
    iclw = (void *)calloc(size_output, sizeof(float));
    tmp_iclw = (double *)calloc(nynx,sizeof(double));
    if(tmp_iclw == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_iclw: Unable to allocate memory for temporary output array");
      return(NhlFATAL);
    }
  }
  else {
    iclw = (void *)calloc(size_output, sizeof(double));
  }
  if(iclw == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_iclw: Unable to allocate memory for output array");
    return(NhlFATAL);
  }

/* 
 * Allocate space for output dimension sizes and set them.
 */
  ndims_iclw = ndims_leftmost + 2;
  dsizes_iclw = (ng_size_t*)calloc(ndims_iclw,sizeof(ng_size_t));  
  if( dsizes_iclw == NULL ) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_iclw: Unable to allocate memory for holding dimension sizes");
    return(NhlFATAL);
  }
  for(i = 0; i < ndims_iclw-2; i++) dsizes_iclw[i] = dsizes_p[i];
  dsizes_iclw[ndims_iclw-2] = ny;
  dsizes_iclw[ndims_iclw-1] = nx;

/*
 * Loop across leftmost dimensions and call the Fortran routine for each
 * subsection of the input arrays..
 */
  index_p = index_iclw = 0;

  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of p (tmp_p) to double if necessary.
 */
    if(type_p != NCL_double) {
      coerce_subset_input_double(p,tmp_p,index_p,type_p,nznynx,0,NULL,NULL);
    }
    else {
      tmp_p = &((double*)p)[index_p];
    }

/*
 * Coerce subsection of qc (tmp_qc) to double if necessary.
 */
    if(type_qc != NCL_double) {
      coerce_subset_input_double(qc,tmp_qc,index_p,type_qc,nznynx,0,NULL,NULL);
    }
    else {
      tmp_qc = &((double*)qc)[index_p];
    }
/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_iclw == NCL_double) tmp_iclw = &((double*)iclw)[index_iclw];

/*
 * Call the Fortran routine.
 */
    NGCALLF(dcomputeiclw,DCOMPUTEICLW)(tmp_iclw, tmp_p, tmp_qc, &inx, &iny, &inz);
/*
 * Coerce output back to float if necessary.
 */
    if(type_iclw == NCL_float) {
      coerce_output_float_only(iclw,tmp_iclw,nynx,index_iclw);
    }
    index_p    += nznynx;
    index_iclw += nynx;
  }

/*
 * Free unneeded memory.
 */
  if(type_p    != NCL_double) NclFree(tmp_p);
  if(type_qc   != NCL_double) NclFree(tmp_qc);
  if(type_iclw != NCL_double) NclFree(tmp_iclw);

/*
 * Set up some attributes ("description" and "units") to return.
 */
  cdescription = (char *)calloc(16,sizeof(char));
  cunits       = (char *)calloc(3,sizeof(char));
  strcpy(cdescription,"Int Cloud Water");
  strcpy(cunits,"mm");
  description = (NclQuark*)NclMalloc(sizeof(NclQuark));
  units       = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *description = NrmStringToQuark(cdescription);
  *units       = NrmStringToQuark(cunits);
  free(cdescription);
  free(cunits);

/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)iclw,
                            NULL,
                            ndims_iclw,
                            dsizes_iclw,
                            TEMPORARY,
                            NULL,
                            type_obj_iclw
                            );

/*
 * Set up attributes to return.
 */
  att_id = _NclAttCreate(NULL,NULL,Ncl_Att,0,NULL);

  dsizes[0] = 1;
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)description,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "description",
             att_md,
             NULL
             );
    
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)units,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "units",
             att_md,
             NULL
             );
    
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          NULL,
                          att_id,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );
  NclFree(dsizes_iclw);
/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);

}

NhlErrorTypes wrf_wetbulb_W( void )
{
/*
 * Input array variables
 */
  void *prs, *tmk, *qvp;
  double *tmp_prs = NULL;
  double *tmp_tmk = NULL;
  double *tmp_qvp = NULL;
  int ndims_prs, ndims_tmk, ndims_qvp;
  ng_size_t dsizes_prs[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_tmk[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_qvp[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_prs, type_tmk, type_qvp;
/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info;

/*
 * Output variable and attributes.
 */
  void *twb;
  NclQuark *description, *units;
  char *cdescription, *cunits;
  double *tmp_twb = NULL;
  NclBasicDataTypes type_twb;
  NclObjClass type_obj_twb;
/*
 * Various
 */
  ng_size_t i, nx, ny, nz, nxyz, size_leftmost, index_prs, size_twb;
  int inx, iny, inz;
  char *psa_file;

/*
 * Variables for returning the output array with attributes attached.
 */
  int att_id;
  ng_size_t dsizes[1];
  NclMultiDValData att_md, return_md;
  NclVar tmp_var;
  NclStackEntry return_data;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
  prs = (void*)NclGetArgValue(
           0,
           3,
           &ndims_prs,
           dsizes_prs,
           NULL,
           NULL,
           &type_prs,
           DONT_CARE);

  tmk = (void*)NclGetArgValue(
           1,
           3,
           &ndims_tmk,
           dsizes_tmk,
           NULL,
           NULL,
           &type_tmk,
           DONT_CARE);

  qvp = (void*)NclGetArgValue(
           2,
           3,
           &ndims_qvp,
           dsizes_qvp,
           NULL,
           NULL,
           &type_qvp,
           DONT_CARE);
/*
 * Error checking. Input variables must be same size.
 */
  if(ndims_prs < 3 || ndims_prs != ndims_tmk || ndims_prs != ndims_qvp ) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_wetbulb: The prs, tmk, and qvp arrays must have at least 3 dimensions and have the same number of dimensions as each other");
    return(NhlFATAL);
  }
  for(i = 0; i < ndims_prs; i++) {
    if(dsizes_prs[i] != dsizes_tmk[i] || dsizes_prs[i] != dsizes_qvp[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_wetbulb: prs, tmk, and qvp must be the same dimensionality");
      return(NhlFATAL);
    }
  }

/*
 * Test dimension sizes.
 */
  nz = dsizes_prs[ndims_prs-1];
  ny = dsizes_prs[ndims_prs-2];
  nx = dsizes_prs[ndims_prs-3];

  if(nx > INT_MAX) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_wetbulb: nx = %ld is greater than INT_MAX", nx);
    return(NhlFATAL);
  }
  if(ny > INT_MAX) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_wetbulb: ny = %ld is greater than INT_MAX", ny);
    return(NhlFATAL);
  }
  if(nz > INT_MAX) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_wetbulb: nz = %ld is greater than INT_MAX", nz);
    return(NhlFATAL);
  }
  inx = (int) nx;
  iny = (int) ny;
  inz = (int) nz;

/*
 * Retrieve dimension names from the "tmk" variable, if any.
 * These dimension names will later be attached to the output variable.
 */
  dim_info = get_wrf_dim_info(1,3,ndims_tmk,dsizes_tmk);

/*
 * Calculate size of leftmost dimensions.
 */
  size_leftmost = 1;
  for(i = 0; i < ndims_prs-3; i++) size_leftmost *= dsizes_prs[i];
  nxyz = nx * ny * nz;
  size_twb = size_leftmost * nxyz;

/* 
 * Allocate space for coercing input arrays.  If the input prs, tmk,
 * or qvp are already double, then we don't need to allocate space
 * for temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 *
 * The output type defaults to float, unless any of the input arrays
 * are double.
 */
  type_twb     = NCL_float;
  type_obj_twb = nclTypefloatClass;
  if(type_prs != NCL_double) {
    tmp_prs = (double *)calloc(nxyz,sizeof(double));
    if(tmp_prs == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_wetbulb: Unable to allocate memory for coercing 'prs' array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_twb     = NCL_double;
    type_obj_twb = nclTypedoubleClass;
  }

  if(type_tmk != NCL_double) {
    tmp_tmk = (double *)calloc(nxyz,sizeof(double));
    if(tmp_tmk == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_wetbulb: Unable to allocate memory for coercing 'tmk' to double");
      return(NhlFATAL);
    }
  }
  else {
    type_twb     = NCL_double;
    type_obj_twb = nclTypedoubleClass;
  }

  if(type_qvp != NCL_double) {
    tmp_qvp = (double *)calloc(nxyz,sizeof(double));
    if(tmp_qvp == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_wetbulb: Unable to allocate memory for coercing 'qvp' to double");
      return(NhlFATAL);
    }
  }
  else {
    type_twb     = NCL_double;
    type_obj_twb = nclTypedoubleClass;
  }

/*
 * Allocate space for output array.
 */ 
  if(type_twb == NCL_double) {
    twb = (double *)calloc(size_twb,sizeof(double));
    if(twb == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_wetbulb: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }
  else {
    twb     = (float *)calloc(size_twb,sizeof(float));
    tmp_twb = (double *)calloc(nxyz,sizeof(double));
    if(tmp_twb == NULL || twb == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_wetbulb: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }
/*
 * Get path to psadilookup.dat file required by this routine. 
 */
  psa_file = get_psa_file();

/*
 * Loop across leftmost dimensions and call the Fortran routine for each
 * three-dimensional subsection.
 */
  index_prs = 0;
  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of p (tmp_prs) to double if necessary.
 */
    if(type_prs != NCL_double) {
      coerce_subset_input_double(prs,tmp_prs,index_prs,type_prs,
                                 nxyz,0,NULL,NULL);
    }
    else {
      tmp_prs = &((double*)prs)[index_prs];
    }
/*
 * Coerce subsection of tmk (tmp_tmk) to double if ncessary.
 */
    if(type_tmk != NCL_double) {
      coerce_subset_input_double(tmk,tmp_tmk,index_prs,type_tmk,
                                 nxyz,0,NULL,NULL);
    }
    else {
      tmp_tmk = &((double*)tmk)[index_prs];
    }

/*
 * Coerce subsection of qvp (tmp_qvp) to double if ncessary.
 */
    if(type_qvp != NCL_double) {
      coerce_subset_input_double(qvp,tmp_qvp,index_prs,type_qvp,
                                 nxyz,0,NULL,NULL);
    }
    else {
      tmp_qvp = &((double*)qvp)[index_prs];
    }

/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_twb == NCL_double) tmp_twb = &((double*)twb)[index_prs];
/*
 * Call Fortran routine.
 */
    NGCALLF(wetbulbcalc,WETBULBCALC)(tmp_prs,tmp_tmk,tmp_qvp,tmp_twb,
				     &inx,&iny,&inz,psa_file,
				     strlen(psa_file));

/*
 * Coerce output back to float if necessary.
 */
    if(type_twb == NCL_float) {
      coerce_output_float_only(twb,tmp_twb,nxyz,index_prs);
    }

    index_prs += nxyz;    /* Increment index */
  }
/*
 * Free up memory.
 */
  if(type_prs != NCL_double) NclFree(tmp_prs);
  if(type_tmk != NCL_double) NclFree(tmp_tmk);
  if(type_qvp != NCL_double) NclFree(tmp_qvp);
  if(type_twb != NCL_double) NclFree(tmp_twb);
  NclFree(psa_file);

/*
 * Set up some attributes ("description" and "units") to return.
 */
  cdescription = (char *)calloc(21,sizeof(char));
  strcpy(cdescription,"Wet Bulb Temperature");
  cunits       = (char *)calloc(2,sizeof(char));
  strcpy(cunits,"C");
  description = (NclQuark*)NclMalloc(sizeof(NclQuark));
  units       = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *description = NrmStringToQuark(cdescription);
  *units       = NrmStringToQuark(cunits);
  free(cunits);
  free(cdescription);
/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)twb,
                            NULL,
                            ndims_prs,
                            dsizes_prs,
                            TEMPORARY,
                            NULL,
                            type_obj_twb
                            );
/*
 * Set up attributes to return.
 */
  att_id = _NclAttCreate(NULL,NULL,Ncl_Att,0,NULL);

  dsizes[0] = 1;
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)description,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "description",
             att_md,
             NULL
             );
    
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)units,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "units",
             att_md,
             NULL
             );
    
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          att_id,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );

  NclFree(dim_info);
/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);
}


NhlErrorTypes wrf_omega_W( void )
{
/*
 * Input array variables
 */
  void *qvp, *tmk, *www, *prs;
  double *tmp_qvp = NULL;
  double *tmp_tmk = NULL;
  double *tmp_www = NULL;
  double *tmp_prs = NULL;
  int ndims_qvp, ndims_tmk, ndims_www, ndims_prs; 
  ng_size_t dsizes_qvp[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_tmk[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_www[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_prs[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_qvp, type_tmk, type_www, type_prs;
/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info;

/*
 * Output variable and attributes.
 */
  void *omg;
  NclQuark *description, *units;
  char *cdescription, *cunits;
  double *tmp_omg = NULL;
  NclBasicDataTypes type_omg;
  NclObjClass type_obj_omg;
/*
 * Various
 */
  ng_size_t i, nx, ny, nz, nxyz, size_leftmost, index_qvp, size_omg;
  int inx, iny, inz;

/*
 * Variables for returning the output array with attributes attached.
 */
  int att_id;
  ng_size_t dsizes[1];
  NclMultiDValData att_md, return_md;
  NclVar tmp_var;
  NclStackEntry return_data;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
  qvp = (void*)NclGetArgValue(
           0,
           4,
           &ndims_qvp,
           dsizes_qvp,
           NULL,
           NULL,
           &type_qvp,
           DONT_CARE);

  tmk = (void*)NclGetArgValue(
           1,
           4,
           &ndims_tmk,
           dsizes_tmk,
           NULL,
           NULL,
           &type_tmk,
           DONT_CARE);

  www = (void*)NclGetArgValue(
           2,
           4,
           &ndims_www,
           dsizes_www,
           NULL,
           NULL,
           &type_www,
           DONT_CARE);

  prs = (void*)NclGetArgValue(
           3,
           4,
           &ndims_prs,
           dsizes_prs,
           NULL,
           NULL,
           &type_prs,
           DONT_CARE);

/*
 * Error checking. Input variables must be same size.
 */
  if(ndims_qvp < 3 || ndims_qvp != ndims_tmk || ndims_qvp != ndims_www ||
     ndims_qvp != ndims_prs) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_omega: qvp, tmk, www, and prs must have at least 3 dimensions and have the same number of dimensions as each other");
    return(NhlFATAL);
  }
  for(i = 0; i < ndims_qvp; i++) {
    if(dsizes_qvp[i] != dsizes_tmk[i] || dsizes_qvp[i] != dsizes_prs[i] || dsizes_qvp[i] != dsizes_www[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_omega: qvp, tmk, www, and prs must be the same dimensionality");
      return(NhlFATAL);
    }
  }

/*
 * Test dimension sizes.
 */
  nz = dsizes_qvp[ndims_qvp-1];
  ny = dsizes_qvp[ndims_qvp-2];
  nx = dsizes_qvp[ndims_qvp-3];

  if(nx > INT_MAX) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_omega: nx = %ld is greater than INT_MAX", nx);
    return(NhlFATAL);
  }
  if(ny > INT_MAX) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_omega: ny = %ld is greater than INT_MAX", ny);
    return(NhlFATAL);
  }
  if(nz > INT_MAX) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_omega: nz = %ld is greater than INT_MAX", nz);
    return(NhlFATAL);
  }
  inx = (int) nx;
  iny = (int) ny;
  inz = (int) nz;

/*
 * Retrieve dimension names from the "tmk" variable, if any.
 * These dimension names will later be attached to the output variable.
 */
  dim_info = get_wrf_dim_info(1,3,ndims_tmk,dsizes_tmk);

/*
 * Calculate size of leftmost dimensions.
 */
  size_leftmost = 1;
  for(i = 0; i < ndims_qvp-3; i++) size_leftmost *= dsizes_qvp[i];
  nz = dsizes_qvp[ndims_qvp-1];
  ny = dsizes_qvp[ndims_qvp-2];
  nx = dsizes_qvp[ndims_qvp-3];
  nxyz = nx * ny * nz;
  size_omg = size_leftmost * nxyz;

/* 
 * Allocate space for coercing input arrays.  If the input qvp, tmk,
 * www, prs are already double, then we don't need to allocate space
 * for temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 *
 * The output type defaults to float, unless any of the input arrays
 * are double.
 */
  type_omg     = NCL_float;
  type_obj_omg = nclTypefloatClass;

  if(type_qvp != NCL_double) {
    tmp_qvp = (double *)calloc(nxyz,sizeof(double));
    if(tmp_qvp == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_omega: Unable to allocate memory for coercing 'qvp' to double");
      return(NhlFATAL);
    }
  }
  else {
    type_omg     = NCL_double;
    type_obj_omg = nclTypedoubleClass;
  }

  if(type_tmk != NCL_double) {
    tmp_tmk = (double *)calloc(nxyz,sizeof(double));
    if(tmp_tmk == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_omega: Unable to allocate memory for coercing 'tmk' to double");
      return(NhlFATAL);
    }
  }
  else {
    type_omg     = NCL_double;
    type_obj_omg = nclTypedoubleClass;
  }

  if(type_www != NCL_double) {
    tmp_www = (double *)calloc(nxyz,sizeof(double));
    if(tmp_www == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_omega: Unable to allocate memory for coercing 'www' to double");
      return(NhlFATAL);
    }
  }
  else {
    type_omg     = NCL_double;
    type_obj_omg = nclTypedoubleClass;
  }

  if(type_prs != NCL_double) {
    tmp_prs = (double *)calloc(nxyz,sizeof(double));
    if(tmp_prs == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_omega: Unable to allocate memory for coercing 'prs' array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_omg     = NCL_double;
    type_obj_omg = nclTypedoubleClass;
  }
/*
 * Allocate space for output array.
 */ 
  if(type_omg == NCL_double) {
    omg = (double *)calloc(size_omg,sizeof(double));
    if(omg == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_omega: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }
  else {
    omg     = (float *)calloc(size_omg,sizeof(float));
    tmp_omg = (double *)calloc(nxyz,sizeof(double));
    if(tmp_omg == NULL || omg == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_omega: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }

/*
 * Loop across leftmost dimensions and call the Fortran routine for each
 * three-dimensional subsection.
 */
  index_qvp = 0;
  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of qvp (tmp_qvp) to double if ncessary.
 */
    if(type_qvp != NCL_double) {
      coerce_subset_input_double(qvp,tmp_qvp,index_qvp,type_qvp,
                                 nxyz,0,NULL,NULL);
    }
    else {
      tmp_qvp = &((double*)qvp)[index_qvp];
    }

/*
 * Coerce subsection of tmk (tmp_tmk) to double if ncessary.
 */
    if(type_tmk != NCL_double) {
      coerce_subset_input_double(tmk,tmp_tmk,index_qvp,type_tmk,
                                 nxyz,0,NULL,NULL);
    }
    else {
      tmp_tmk = &((double*)tmk)[index_qvp];
    }

/*
 * Coerce subsection of www (tmp_www) to double if ncessary.
 */
    if(type_www != NCL_double) {
      coerce_subset_input_double(www,tmp_www,index_qvp,type_www,
                                 nxyz,0,NULL,NULL);
    }
    else {
      tmp_www = &((double*)www)[index_qvp];
    }

/*
 * Coerce subsection of p (tmp_prs) to double if necessary.
 */
    if(type_prs != NCL_double) {
      coerce_subset_input_double(prs,tmp_prs,index_qvp,type_prs,
                                 nxyz,0,NULL,NULL);
    }
    else {
      tmp_prs = &((double*)prs)[index_qvp];
    }
/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_omg == NCL_double) tmp_omg = &((double*)omg)[index_qvp];
/*
 * Call Fortran routine.
 */
    NGCALLF(omgcalc,OMGCALC)(tmp_qvp,tmp_tmk,tmp_www,tmp_prs,tmp_omg,
                             &inx,&iny,&inz);

/*
 * Coerce output back to float if necessary.
 */
    if(type_omg == NCL_float) {
      coerce_output_float_only(omg,tmp_omg,nxyz,index_qvp);
    }

    index_qvp += nxyz;    /* Increment index */
  }

/*
 * Free up memory.
 */
  if(type_qvp != NCL_double) NclFree(tmp_qvp);
  if(type_tmk != NCL_double) NclFree(tmp_tmk);
  if(type_www != NCL_double) NclFree(tmp_www);
  if(type_prs != NCL_double) NclFree(tmp_prs);
  if(type_omg != NCL_double) NclFree(tmp_omg);

/*
 * Set up some attributes ("description" and "units") to return.
 */
  cdescription = (char *)calloc(6,sizeof(char));
  strcpy(cdescription,"Omega");
  cunits       = (char *)calloc(5,sizeof(char));
  strcpy(cunits,"Pa/s");
  description = (NclQuark*)NclMalloc(sizeof(NclQuark));
  units       = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *description = NrmStringToQuark(cdescription);
  *units       = NrmStringToQuark(cunits);
  free(cunits);
  free(cdescription);

/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)omg,
                            NULL,
                            ndims_qvp,
                            dsizes_qvp,
                            TEMPORARY,
                            NULL,
                            type_obj_omg
                            );
/*
 * Set up attributes to return.
 */
  att_id = _NclAttCreate(NULL,NULL,Ncl_Att,0,NULL);

  dsizes[0] = 1;
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)description,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "description",
             att_md,
             NULL
             );
    
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)units,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "units",
             att_md,
             NULL
             );
    
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          att_id,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );

  NclFree(dim_info);
/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);
         }


NhlErrorTypes wrf_virtual_temp_W( void )
{
/*
 * Input array variables
 */
  void *temp, *ratmx;
  double *tmp_temp = NULL;
  double *tmp_ratmx = NULL;
  int ndims_temp, ndims_ratmx;
  ng_size_t dsizes_temp[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_ratmx[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_temp, type_ratmx;
/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info;

/*
 * Output variable and attributes.
 */
  void *tv;
  NclQuark *description, *units;
  char *cdescription, *cunits;
  double *tmp_tv = NULL;
  NclBasicDataTypes type_tv;
  NclObjClass type_obj_tv;
/*
 * Various
 */
  ng_size_t i, nx, ny, nz, nxyz, size_leftmost, index_temp, size_tv;
  int inx, iny, inz;

/*
 * Variables for returning the output array with attributes attached.
 */
  int att_id;
  ng_size_t dsizes[1];
  NclMultiDValData att_md, return_md;
  NclVar tmp_var;
  NclStackEntry return_data;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
  temp = (void*)NclGetArgValue(
           0,
           2,
           &ndims_temp,
           dsizes_temp,
           NULL,
           NULL,
           &type_temp,
           DONT_CARE);

  ratmx = (void*)NclGetArgValue(
           1,
           2,
           &ndims_ratmx,
           dsizes_ratmx,
           NULL,
           NULL,
           &type_ratmx,
           DONT_CARE);
/*
 * Error checking. Input variables must be same size.
 */
  if(ndims_temp < 3 || ndims_temp != ndims_ratmx) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_virtual_temp: The temp and ratmx arrays must have at least 3 dimensions and have the same number of dimensions as each other");
    return(NhlFATAL);
  }
  for(i = 0; i < ndims_temp; i++) {
    if(dsizes_temp[i] != dsizes_ratmx[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_virtual_temp: temp and ratmx must be the same dimensionality");
      return(NhlFATAL);
    }
  }

/*
 * Test dimension sizes.
 */
  nx = dsizes_temp[ndims_temp-1];
  ny = dsizes_temp[ndims_temp-2];
  nz = dsizes_temp[ndims_temp-3];

  if(nx > INT_MAX) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_virtual_temp: nx = %ld is greater than INT_MAX", nx);
    return(NhlFATAL);
  }
  if(ny > INT_MAX) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_virtual_temp: ny = %ld is greater than INT_MAX", ny);
    return(NhlFATAL);
  }
  if(nz > INT_MAX) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_virtual_temp: nz = %ld is greater than INT_MAX", nz);
    return(NhlFATAL);
  }
  inx = (int) nx;
  iny = (int) ny;
  inz = (int) nz;

/*
 * Retrieve dimension names from the "temp" variable, if any.
 * These dimension names will later be attached to the output variable.
 */
  dim_info = get_wrf_dim_info(1,3,ndims_temp,dsizes_temp);

/*
 * Calculate size of leftmost dimensions.
 */
  size_leftmost = 1;
  for(i = 0; i < ndims_temp-3; i++) size_leftmost *= dsizes_temp[i];
  nxyz = nx * ny * nz;
  size_tv = size_leftmost * nxyz;

/* 
 * Allocate space for coercing input arrays.  If the input temp
 * or ratmx are already double, then we don't need to allocate space
 * for temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 *
 * The output type defaults to float, unless any of the input arrays
 * are double.
 */
  type_tv     = NCL_float;
  type_obj_tv = nclTypefloatClass;

  if(type_temp != NCL_double) {
    tmp_temp = (double *)calloc(nxyz,sizeof(double));
    if(tmp_temp == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_virtual_temp: Unable to allocate memory for coercing 'temp' to double");
      return(NhlFATAL);
    }
  }
  else {
    type_tv     = NCL_double;
    type_obj_tv = nclTypedoubleClass;
  }

  if(type_ratmx != NCL_double) {
    tmp_ratmx = (double *)calloc(nxyz,sizeof(double));
    if(tmp_ratmx == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_virtual_temp: Unable to allocate memory for coercing 'ratmx' to double");
      return(NhlFATAL);
    }
  }
  else {
    type_tv     = NCL_double;
    type_obj_tv = nclTypedoubleClass;
  }

/*
 * Allocate space for output array.
 */ 
  if(type_tv == NCL_double) {
    tv = (double *)calloc(size_tv,sizeof(double));
    if(tv == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_virtual_temp: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }
  else {
    tv     = (float *)calloc(size_tv,sizeof(float));
    tmp_tv = (double *)calloc(nxyz,sizeof(double));
    if(tmp_tv == NULL || tv == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_virtual_temp: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }
/*
 * Loop across leftmost dimensions and call the Fortran routine for each
 * three-dimensional subsection.
 */
  index_temp = 0;
  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of temp (tmp_temp) to double if ncessary.
 */
    if(type_temp != NCL_double) {
      coerce_subset_input_double(temp,tmp_temp,index_temp,type_temp,
                                 nxyz,0,NULL,NULL);
    }
    else {
      tmp_temp = &((double*)temp)[index_temp];
    }

/*
 * Coerce subsection of ratmx (tmp_ratmx) to double if ncessary.
 */
    if(type_ratmx != NCL_double) {
      coerce_subset_input_double(ratmx,tmp_ratmx,index_temp,type_ratmx,
                                 nxyz,0,NULL,NULL);
    }
    else {
      tmp_ratmx = &((double*)ratmx)[index_temp];
    }

/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_tv == NCL_double) tmp_tv = &((double*)tv)[index_temp];
/*
 * Call Fortran routine.
 */
    NGCALLF(virtual_temp,VIRTUAL_TEMP)(tmp_temp,tmp_ratmx,tmp_tv,
				       &inx,&iny,&inz);

/*
 * Coerce output back to float if necessary.
 */
    if(type_tv == NCL_float) {
      coerce_output_float_only(tv,tmp_tv,nxyz,index_temp);
    }

    index_temp += nxyz;    /* Increment index */
  }
/*
 * Free up memory.
 */
  if(type_temp  != NCL_double) NclFree(tmp_temp);
  if(type_ratmx != NCL_double) NclFree(tmp_ratmx);
  if(type_tv    != NCL_double) NclFree(tmp_tv);

/*
 * Set up some attributes ("description" and "units") to return.
 */
  cdescription = (char *)calloc(20,sizeof(char));
  strcpy(cdescription,"Virtual Temperature");
  cunits       = (char *)calloc(2,sizeof(char));
  strcpy(cunits,"K");
  description = (NclQuark*)NclMalloc(sizeof(NclQuark));
  units       = (NclQuark*)NclMalloc(sizeof(NclQuark));
  *description = NrmStringToQuark(cdescription);
  *units       = NrmStringToQuark(cunits);
  free(cunits);
  free(cdescription);

/*
 * Set up return value.
 */
  return_md = _NclCreateVal(
                            NULL,
                            NULL,
                            Ncl_MultiDValData,
                            0,
                            (void*)tv,
                            NULL,
                            ndims_temp,
                            dsizes_temp,
                            TEMPORARY,
                            NULL,
                            type_obj_tv
                            );
/*
 * Set up attributes to return.
 */
  att_id = _NclAttCreate(NULL,NULL,Ncl_Att,0,NULL);

  dsizes[0] = 1;
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)description,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "description",
             att_md,
             NULL
             );
    
  att_md = _NclCreateVal(
                         NULL,
                         NULL,
                         Ncl_MultiDValData,
                         0,
                         (void*)units,
                         NULL,
                         1,
                         dsizes,
                         TEMPORARY,
                         NULL,
                         (NclObjClass)nclTypestringClass
                         );
  _NclAddAtt(
             att_id,
             "units",
             att_md,
             NULL
             );
    
  tmp_var = _NclVarCreate(
                          NULL,
                          NULL,
                          Ncl_Var,
                          0,
                          NULL,
                          return_md,
                          dim_info,
                          att_id,
                          NULL,
                          RETURNVAR,
                          NULL,
                          TEMPORARY
                          );

  NclFree(dim_info);
/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);
}



/*
 * This routine gets the path to the psadilookup.dat file
 * required by some WRF routines. 
 *
 * The default is to use $NCARG_ROOT/lib/ncarg/data/asc/psadilookup.dat
 * for the input data file, unless PSADILOOKUP_PATH is set by the
 * user, then it will try to use this path. 
 */
char *get_psa_file()
{
  const char *path = NULL;
  char *psa_file;
  int path_len;

  path = getenv("PSADILOOKUP_PATH");
  if ((void *)path == (void *)NULL) {
   path = _NGGetNCARGEnv("data");
    if ((void *)path != (void *)NULL) {
      path_len = strlen(path) + 21;   /* 21 = "/asc/psadilookup.dat\0" */
      psa_file = malloc(path_len*sizeof(char));
      strcpy(psa_file,path);
      strcat(psa_file,_NhlPATHDELIMITER);
      strcat(psa_file,"asc");
    }
  }
  else {
    strcpy(psa_file,path);
    path_len = strlen(path) + 17;   /* 17 = "/psadilookup.dat\0" */
    psa_file = malloc(path_len*sizeof(char));
  }
  strcat(psa_file,_NhlPATHDELIMITER);
  strcat(psa_file,"psadilookup.dat");
  strcat(psa_file,"\0");
  return(psa_file);
}

/*
 * This routine sets all values of var < 0 to 0.0. This is
 * so you don't have to do this in the NCL script. It's the
 * equivalent of:
 *
 * tmp_var = tmp_var > 0.0
 *
 */
void var_zero(double *tmp_var, ng_size_t n)
{
  ng_size_t i;

  for(i = 0; i < n; i++) {
    if(tmp_var[i] < 0.0) tmp_var[i] = 0.0;
  }
}


/* Converts from hPa to Pa. */

void convert_to_hPa(double *pp, ng_size_t np)
{
  ng_size_t i;
  for(i = 0; i < np; i++) pp[i] *= 0.01;
}


/*
 * This procedure flips the given double array in the
 * leftmost dimension, given the size of the leftmost
 * dimension, and the product of the rightmost two dimensions.
 */
void flip_it(double *tmp_from, double *tmp_to, ng_size_t nz, ng_size_t nynx)
{
  ng_size_t i, index_from, index_to, size_copy;

  size_copy = nynx*sizeof(double);
  for(i = 0; i < nz; i++) {
    index_from = (i * nynx) * sizeof(double);
    index_to   = ((nz-1-i) * nynx) * sizeof(double);
    (void *)memcpy((void*)((char*)tmp_to)   + index_to,
                   (void*)((char*)tmp_from) + index_from,size_copy);
  }
}
