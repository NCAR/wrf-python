#include <stdio.h>
#include <math.h>
#include "wrapper.h"

#define ERRLEN 512

extern void NGCALLF(dcapecalc3d,DCAPECALC3D)(double *prs, double *tmk, 
                                             double *qvp, double *ght,
                                             double *ter, double *sfp, 
                                             double *cape, double *cin, 
                                             double *cmsg,
                                             int *miy, int *mjx, int *mkzh, 
                                             int *i3dflag, int *ter_follow,
                                             char *,int);


/*
 * Function for calculating cape (from the RIP code). This function
 * depends on the "psadilookup.dat" file, which by default will be
 * searched for in $NCARG_ROOT/lib/ncarg/data/asc/), unless
 * NCARG_PSADILOOKUP is set to the location of this file.
 */

/*
 * The rip_cape_3d wrapper is for the case where I3DFLAG is set to
 * 1 in the Fortran rip_cape.f file.
 */
NhlErrorTypes rip_cape_3d_W( void )
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
  int ndims_p, ndims_t, ndims_q, ndims_z, ndims_zsfc, ndims_psfc;
  ng_size_t dsizes_p[NCL_MAX_DIMENSIONS], dsizes_t[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_q[NCL_MAX_DIMENSIONS], dsizes_z[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_zsfc[NCL_MAX_DIMENSIONS], dsizes_psfc[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_p, type_t, type_q, type_z, type_zsfc, type_psfc;

/*
 * Output array variables
 */
  void *cape;
  double *tmp_cape = NULL, cmsg;
  double *tmp_cin = NULL;
  NclBasicDataTypes type_cape;
  int ndims_cape;
  ng_size_t *dsizes_cape;
  NclScalar missing_cape;
/*
 * File input variables.
 */
  const char *path = NULL;
  char psa_file[_NhlMAXFNAMELEN];
  int errstat;
  char *errmsg;

/*
 * Declare various variables for random purposes.
 */
  ng_size_t i;
  ng_size_t miy = 0;
  ng_size_t mjx = 0;
  ng_size_t mkzh = 0;
  ng_size_t ntime = 0;
  ng_size_t size_cape, size_output, size_zsfc;
  ng_size_t index_cape, index_zsfc, index_cin;
  int i3dflag=1, scalar_zsfc;
  int iter, ret;
  int imiy, imjx, imkzh;

/*
 * The default is to use $NCARG_ROOT/lib/ncarg/data/asc/psadilookup.dat
 * for the input data file, unless PSADILOOKUP_PATH is set by the
 * user, then it will try to use this path. 
 */
  path = getenv("PSADILOOKUP_PATH");
  if ((void *)path == (void *)NULL) {
    path = _NGGetNCARGEnv("data");
    if ((void *)path != (void *)NULL) {
      strcpy(psa_file,path);
      strcat(psa_file,_NhlPATHDELIMITER);
      strcat(psa_file,"asc");
      strcat(psa_file,_NhlPATHDELIMITER);
      strcat(psa_file,"psadilookup.dat");
    }
  }
  else {
    strcpy(psa_file,path);
    strcat(psa_file,_NhlPATHDELIMITER);
    strcat(psa_file,"psadilookup.dat");
  }

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
 *  - p,t,q,z (time,lev,lat,lon) and psfc,zsfc (time,lat,lon)
 *  - p,t,q,z (lev,lat,lon) and psfc,zsfc (lat,lon)
 *  - p,t,q,z (lev) and psfc,zsfc (scalars)
 */
  if(ndims_p != ndims_t || ndims_p != ndims_q || ndims_p != ndims_z) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_3d: The p, t, q, and z arrays must all have the same number of dimensions");
    return(NhlFATAL);
  }
  if(ndims_p != 1 && ndims_p != 3 && ndims_p != 4) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_3d: The p, t, q, and z arrays must be 1-, 3-, or 4-dimensional\n");
    return(NhlFATAL);
  }
/*
 * zsfc and psfc can be scalars, if the other input arrays are 1D.
 */
  scalar_zsfc = is_scalar(ndims_zsfc,dsizes_zsfc);

  if((ndims_zsfc != ndims_psfc) || (scalar_zsfc && ndims_p != 1) || 
     (!scalar_zsfc && ndims_zsfc != ndims_p-1)) { 
    NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_3d: The zsfc and psfc arrays must have the same number of dimensions, and either be scalars or one less dimension than the other input arrays");
    return(NhlFATAL);
  }

/*
 * Now check that the dimension sizes are equal to each other.
 */
  for(i = 0; i < ndims_p; i++) {
    if(dsizes_p[i] != dsizes_t[i] || dsizes_p[i] != dsizes_q[i] || 
       dsizes_p[i] != dsizes_z[i]) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_3d: p, t, q, and z must be the same dimensionality");
    return(NhlFATAL);
    }
  }

  for(i = 0; i < ndims_psfc; i++) {
    if(dsizes_psfc[i] != dsizes_zsfc[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_3d: psfc and zsfc must be the same dimensionality");
      return(NhlFATAL);
    }
  }
/*
 * Get sizes of input arrays.
 */
  if(ndims_p == 4) {
    ntime = dsizes_p[0];          /* time, serves as a leftmost dimension */
    mkzh  = dsizes_p[1];          /* lev */
    mjx   = dsizes_p[2];          /* lat */
    miy   = dsizes_p[3];          /* lon */
  }
  else if(ndims_p == 3) {
    ntime = 1;
    mkzh = dsizes_p[0];           /* lev */
    mjx  = dsizes_p[1];           /* lat */
    miy  = dsizes_p[2];           /* lon */
  }
  else if(ndims_p == 1) {
    ntime = 1;
    mkzh = dsizes_p[0];           /* lev */
    mjx  = 1;                     /* lat */
    miy  = 1;                     /* lon */
  }

/*
 * Test input dimension sizes.
 */
  if((miy > INT_MAX) || (mjx > INT_MAX) || (mkzh > INT_MAX)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_3d: one or more input dimension sizes is greater than INT_MAX");
    return(NhlFATAL);
  }
  imiy = (int) miy;
  imjx = (int) mjx;
  imkzh = (int) mkzh;

/*
 * Check some more dimension sizes.
 */
  if(ndims_p == 4) {
    if(dsizes_psfc[0] != ntime || dsizes_psfc[1] != mjx || 
       dsizes_psfc[2] != miy) { 
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_3d: If p,q,t,z are 4-dimensional (time x lev x lat x lon), psfc,zsfc must be 3-dimensional (time x lat x lon)");
      return(NhlFATAL);
    }
  }
  if(ndims_p == 3) {
    if(dsizes_psfc[0] != mjx || dsizes_psfc[1] != miy) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_3d: If p,q,t,z are 3-dimensional (time x lev x lat x lon), psfc,zsfc must be 2-dimensional (lat x lon)");
      return(NhlFATAL);
    }
  }
/*
 * Calculate size of output array. The output array size depends on
 * the size of p,t,q,z:
 *
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
    NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_3d: Unable to allocate memory for array dimensionality");
    return(NhlFATAL);
  }

  dsizes_cape[0] = 2;                /* 0 = cape, 1 = cin */
  for(i = 0; i < ndims_p; i++ ) {
    dsizes_cape[i+1] = dsizes_p[i];
  }
  size_zsfc   = mjx * miy;
  size_cape   = mkzh * size_zsfc;       /* Also size of cin array */
  size_output = 2 * size_cape * ntime;

/* 
 * Allocate space for output arrays.  If any of the input is already double,
 * then we don't need to allocate space for temporary arrays, because
 * we'll just change the pointer into the void array appropriately.
 */
  if(type_p == NCL_double || type_t == NCL_double || type_q == NCL_double ||
     type_z == NCL_double) {
    type_cape = NCL_double;
    cape = (double *)calloc(size_output,sizeof(double));
    missing_cape.doubleval = ((NclTypeClass)nclTypedoubleClass)->type_class.default_mis.doubleval;
    cmsg = missing_cape.doubleval;
    if(cape == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_3d: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }
  else {
    type_cape = NCL_float;
    cape      = (float *)calloc(size_output,sizeof(float));
    tmp_cape  = (double *)calloc(size_cape,sizeof(double));
    tmp_cin   = (double *)calloc(size_cape,sizeof(double));
    if(cape == NULL || tmp_cape == NULL || tmp_cin == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_3d: Unable to allocate memory for output arrays");
      return(NhlFATAL);
    }
    missing_cape.floatval = ((NclTypeClass)nclTypefloatClass)->type_class.default_mis.floatval;
    cmsg = (double)missing_cape.floatval;
  }

/*
 * Allocate memory for allocating input arrays to double, if necessary.
 */
  if(type_p != NCL_double) {
    tmp_p = (double *)calloc(size_cape,sizeof(double));
    if(tmp_p == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_3d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  if(type_t != NCL_double) {
    tmp_t = (double *)calloc(size_cape,sizeof(double));
    if(tmp_t == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_3d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  if(type_q != NCL_double) {
    tmp_q = (double *)calloc(size_cape,sizeof(double));
    if(tmp_q == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_3d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  if(type_z != NCL_double) {
    tmp_z = (double *)calloc(size_cape,sizeof(double));
    if(tmp_z == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_3d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  if(type_zsfc != NCL_double) {
    tmp_zsfc = (double *)calloc(size_zsfc,sizeof(double));
    if(tmp_zsfc == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_3d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  if(type_psfc != NCL_double) {
    tmp_psfc = (double *)calloc(size_zsfc,sizeof(double));
    if(tmp_psfc == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_3d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  /* Allocate space for errmsg*/
  errmsg = (char *) calloc(ERRLEN, sizeof(char))

/*
 * Call the Fortran routine.
 */ 
  index_cape = index_zsfc = 0;
  index_cin = ntime * size_cape;    /* Second half of output array */

  for(i = 0; i < ntime; i++) {
/*
 * Coerce subset of input arrays to double if necessary.
 */
    if(type_p != NCL_double) {
      coerce_subset_input_double(p,tmp_p,index_cape,type_p,size_cape,0,NULL,NULL);
    }
    else {
/*
 * Point tmp_p to appropriate location in p.
 */
      tmp_p = &((double*)p)[index_cape];
    }
    if(type_t != NCL_double) {
      coerce_subset_input_double(t,tmp_t,index_cape,type_t,size_cape,0,NULL,NULL);
    }
    else {
/*
 * Point tmp_t to appropriate location in t.
 */
      tmp_t = &((double*)t)[index_cape];
    }
    if(type_q != NCL_double) {
      coerce_subset_input_double(q,tmp_q,index_cape,type_q,size_cape,0,NULL,NULL);
    }
    else {
/*
 * Point tmp_q to appropriate location in q.
 */
      tmp_q = &((double*)q)[index_cape];
    }
    if(type_z != NCL_double) {
      coerce_subset_input_double(z,tmp_z,index_cape,type_z,size_cape,0,NULL,NULL);
    }
    else {
/*
 * Point tmp_z to appropriate location in z.
 */
      tmp_z = &((double*)z)[index_cape];
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
 * Point tmp_cape and tmp_cin to appropriate location in cape
 * if necessary
 */
    if(type_cape == NCL_double) {
      tmp_cape = &((double*)cape)[index_cape];
      tmp_cin  = &((double*)cape)[index_cin];
    }
    

   errstat = 0;
   errmsg = "";
/*
 * Call Fortran routine.
 */
    NGCALLF(dcapecalc3d,DCAPECALC3D)(tmp_p, tmp_t, tmp_q, tmp_z, tmp_zsfc,
                                     tmp_psfc, tmp_cape, tmp_cin, &cmsg, 
                                     &imiy, &imjx, &imkzh, &i3dflag, &iter,
                                     psa_file,errstat,errmsg,strlen(psa_file),
									 ERRLEN);

/* Terminate if there was an error */
	if (errstat != 0) {
		fprintf(stderr, errmsg);
		exit(errstat);
	}

/*
 * If the output is to be float, then do the coercion here.
 */
    if(type_cape == NCL_float) {
      coerce_output_float_only(cape,tmp_cape,size_cape,index_cape);
      coerce_output_float_only(cape,tmp_cin,size_cape,index_cin);
    }
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
  if(type_p != NCL_double) NclFree(tmp_p);
  if(type_t != NCL_double) NclFree(tmp_t);
  if(type_q != NCL_double) NclFree(tmp_q);
  if(type_z != NCL_double) NclFree(tmp_z);
  if(type_zsfc != NCL_double) NclFree(tmp_zsfc);
  if(type_psfc != NCL_double) NclFree(tmp_psfc);
  if(type_cape != NCL_double) NclFree(tmp_cape);
  if(type_cape != NCL_double) NclFree(tmp_cin);
  NclFree(errmsg);
/*
 * Set up variable to return.
 */
  ret = NclReturnValue(cape,ndims_cape,dsizes_cape,&missing_cape,type_cape,0);
  NclFree(dsizes_cape);
  return(ret);
}


/*
 * The rip_cape_2d wrapper is for the case where I3DFLAG is set to
 * 0 in the Fortran rip_cape.f file.  In this case, 4 2D arrays
 * are returned: cape, cin, lcl, and lfc, but they are all returned 
 * in one big array whose leftmost dimension is 4:
 *
 *   index 0 = cape
 *   index 1 = cin
 *   index 2 = lcl
 *   index 3 = lfc
 */
NhlErrorTypes rip_cape_2d_W( void )
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
  int ndims_p, ndims_t, ndims_q, ndims_z, ndims_zsfc, ndims_psfc;
  ng_size_t dsizes_p[NCL_MAX_DIMENSIONS], dsizes_t[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_q[NCL_MAX_DIMENSIONS], dsizes_z[NCL_MAX_DIMENSIONS];
  ng_size_t dsizes_zsfc[NCL_MAX_DIMENSIONS], dsizes_psfc[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_p, type_t, type_q, type_z, type_zsfc, type_psfc;

/*
 * Output array variables
 */
  void *cape;
  double *tmp_cape = NULL, cmsg;
  double *tmp_cin = NULL;
  NclBasicDataTypes type_cape;
  int ndims_cape = 0;
  NclScalar missing_cape;
  ng_size_t *dsizes_cape;
/*
 * File input variables.
 */
  const char *path = NULL;
  char psa_file[_NhlMAXFNAMELEN];
  int errstat;
  char *errmsg;

/*
 * Declare various variables for random purposes.
 */
  ng_size_t i;
  ng_size_t miy = 0;
  ng_size_t mjx = 0;
  ng_size_t mkzh = 0;
  ng_size_t ntime = 0;
  ng_size_t size_cape, size_output, size_zsfc;
  ng_size_t size_left_zsfc;
  int i3dflag=0;
  ng_size_t index_cape, index_zsfc;
  ng_size_t index_output_cape, index_output_cin, index_output_lcl;
  ng_size_t index_output_lfc, mkzh0_index, mkzh1_index, mkzh2_index;
  int iter, ret;
  int imiy, imjx, imkzh;

/*
 * The default is to use $NCARG_ROOT/lib/ncarg/data/asc/psadilookup.dat
 * for the input data file, unless PSADILOOKUP_PATH is set by the
 * user, then it will try to use this path. 
 */
  path = getenv("PSADILOOKUP_PATH");
  if ((void *)path == (void *)NULL) {
    path = _NGGetNCARGEnv("data");
    if ((void *)path != (void *)NULL) {
      strcpy(psa_file,path);
      strcat(psa_file,_NhlPATHDELIMITER);
      strcat(psa_file,"asc");
      strcat(psa_file,_NhlPATHDELIMITER);
      strcat(psa_file,"psadilookup.dat");
    }
  }
  else {
    strcpy(psa_file,path);
    strcat(psa_file,_NhlPATHDELIMITER);
    strcat(psa_file,"psadilookup.dat");
  }

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
 * Check the input dimension sizes. There are two possible cases
 * for the input dimension sizes:
 *
 *  - p,t,q,z (time,lev,lat,lon) and psfc,zsfc (time,lat,lon)
 *  - p,t,q,z (lev,lat,lon) and psfc,zsfc (lat,lon)
 */
  if(ndims_p != ndims_t || ndims_p != ndims_q || ndims_p != ndims_z) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_2d: The p, t, q, and z arrays must all have the same number of dimensions");
    return(NhlFATAL);
  }
  if(ndims_p != 3 && ndims_p != 4) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_2d: The p, t, q, and z arrays must be 3 or 4-dimensional\n");
    return(NhlFATAL);
  }
/*
 * Check zsfc and psfc dimension sizes.
 */
  if((ndims_zsfc != ndims_psfc) || (ndims_zsfc != ndims_p-1)) { 
    NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_2d: The zsfc and psfc arrays must have the same number of dimensions and be one less dimension than the other input arrays");
    return(NhlFATAL);
  }

/*
 * Now check that the dimension sizes are equal to each other.
 */
  for(i = 0; i < ndims_p; i++) {
    if(dsizes_p[i] != dsizes_t[i] || dsizes_p[i] != dsizes_q[i] || 
       dsizes_p[i] != dsizes_z[i]) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_2d: p, t, q, and z must be the same dimensionality");
    return(NhlFATAL);
    }
  }

  for(i = 0; i < ndims_psfc; i++) {
    if(dsizes_psfc[i] != dsizes_zsfc[i]) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_2d: psfc and zsfc must be the same dimensionality");
      return(NhlFATAL);
    }
  }
  if(ndims_p == 4) {
/*
 * Store dimension sizes.
 */
    ntime = dsizes_p[0];       /* time */
    mkzh = dsizes_p[1];        /* lev */
    mjx  = dsizes_p[2];        /* lat */
    miy  = dsizes_p[3];        /* lon */
    ndims_cape = 4;
    if(dsizes_psfc[0] != ntime || dsizes_psfc[1] != mjx ||
       dsizes_psfc[2] != miy) { 
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_2d: If p,q,t,z are 4-dimensional (time x lev x lat x lon), psfc,zsfc must be 3-dimensional (time x lat x lon)");
      return(NhlFATAL);

    }
  }
  else if(ndims_p == 3) {
/*
 * Store dimension sizes.
 */
    ntime = 1;
    mkzh = dsizes_p[0];           /* lev */
    mjx  = dsizes_p[1];           /* lat */
    miy  = dsizes_p[2];           /* lon */
    ndims_cape = 3;
    if(dsizes_psfc[0] != mjx || dsizes_psfc[1] != miy) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_2d: If p,q,t,z are 3-dimensional (time x lev x lat x lon), psfc,zsfc must be 2-dimensional (lat x lon)");
      return(NhlFATAL);
    }
  }
/*
 * If mkzh is not at least size 3, then this dimension won't be big 
 * enough to contain the cin, lcl, and lfc values.
 */
  if(mkzh < 3) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_2d: The level dimension must have at least 3 elements");
    return(NhlFATAL);
  }

/*
 * Test input dimension sizes.
 */
  if((miy > INT_MAX) || (mjx > INT_MAX) || (mkzh > INT_MAX)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_2d: one or more input dimension sizes is greater than INT_MAX");
    return(NhlFATAL);
  }
  imiy = (int) miy;
  imjx = (int) mjx;
  imkzh = (int) mkzh;


/*
 * Calculate size of output array. The output array size depends on
 * the size of p,t,q,z:
 *
 *  - p,t,q,z (time,lev,lat,lon) and psfc,zsfc (time,lat,lon)
 *       output array: (4,time,lat,lon)
 *  - p,t,q,z (lev,lat,lon) and psfc,zsfc (lat,lon)
 *       output array: (4,lat,lon)
 */
  dsizes_cape = (ng_size_t *)calloc(ndims_cape,sizeof(ng_size_t));
  if(dsizes_cape == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_2d: Unable to allocate memory for array dimensionality");
    return(NhlFATAL);
  }

  dsizes_cape[0]            = 4;    /* To hold the 4 different variables. */
                                    /* 0=cape, 1=cin, 2=lcl, 3=lfc */
  dsizes_cape[ndims_cape-1] = miy;
  dsizes_cape[ndims_cape-2] = mjx;
  if(ndims_cape == 4) dsizes_cape[1] = ntime;

  size_zsfc   = mjx * miy;
  size_cape   = mkzh * size_zsfc;
  mkzh0_index = (mkzh-1) * size_zsfc;    /* Indexes into cin array for   */
  mkzh1_index = (mkzh-2) * size_zsfc;    /* returning cin, lcl, and lfc  */
  mkzh2_index = (mkzh-3) * size_zsfc;    /* respectively. */
  size_left_zsfc = size_zsfc * ntime;
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
    type_cape = NCL_double;
    cape      = (double *)calloc(size_output,sizeof(double));
    missing_cape.doubleval = ((NclTypeClass)nclTypedoubleClass)->type_class.default_mis.doubleval;
    cmsg = missing_cape.doubleval;
  }
  else {
    type_cape = NCL_float;
    cape      = (float *)calloc(size_output,sizeof(float));
    missing_cape.floatval = ((NclTypeClass)nclTypefloatClass)->type_class.default_mis.floatval;
    cmsg = (double)missing_cape.floatval;
  }
  tmp_cape = (double *)calloc(size_cape,sizeof(double));
  tmp_cin  = (double *)calloc(size_cape,sizeof(double));
  if(cape == NULL || tmp_cape == NULL || tmp_cin == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_2d: Unable to allocate memory for output arrays");
    return(NhlFATAL);
  }

/*
 * Allocate memory for allocating input arrays to double, if necessary.
 */
  if(type_p != NCL_double) {
    tmp_p = (double *)calloc(size_cape,sizeof(double));
    if(tmp_p == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_2d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  if(type_t != NCL_double) {
    tmp_t = (double *)calloc(size_cape,sizeof(double));
    if(tmp_t == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_2d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  if(type_q != NCL_double) {
    tmp_q = (double *)calloc(size_cape,sizeof(double));
    if(tmp_q == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_2d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  if(type_z != NCL_double) {
    tmp_z = (double *)calloc(size_cape,sizeof(double));
    if(tmp_z == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_2d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  if(type_zsfc != NCL_double) {
    tmp_zsfc = (double *)calloc(size_zsfc,sizeof(double));
    if(tmp_zsfc == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_2d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  if(type_psfc != NCL_double) {
    tmp_psfc = (double *)calloc(size_zsfc,sizeof(double));
    if(tmp_psfc == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"rip_cape_2d: Unable to allocate memory for coercing input arrays to double");
      return(NhlFATAL);
    }
  }

  /* Allocate space for errmsg*/
  errmsg = (char *) calloc(ERRLEN, sizeof(char))

/*
 * Call the Fortran routine.
 */ 
  index_cape        = index_zsfc = 0;
  index_output_cape = 0;
  index_output_cin  = size_left_zsfc;
  index_output_lcl  = 2 * size_left_zsfc;
  index_output_lfc  = 3 * size_left_zsfc;

  for(i = 0; i < ntime; i++) {
/*
 * Coerce subset of input arrays to double if necessary.
 */
    if(type_p != NCL_double) {
      coerce_subset_input_double(p,tmp_p,index_cape,type_p,size_cape,0,NULL,NULL);
    }
    else {
/*
 * Point tmp_p to appropriate location in p.
 */
      tmp_p = &((double*)p)[index_cape];
    }
    if(type_t != NCL_double) {
      coerce_subset_input_double(t,tmp_t,index_cape,type_t,size_cape,0,NULL,NULL);
    }
    else {
/*
 * Point tmp_t to appropriate location in t.
 */
      tmp_t = &((double*)t)[index_cape];
    }
    if(type_q != NCL_double) {
      coerce_subset_input_double(q,tmp_q,index_cape,type_q,size_cape,0,NULL,NULL);
    }
    else {
/*
 * Point tmp_q to appropriate location in q.
 */
      tmp_q = &((double*)q)[index_cape];
    }
    if(type_z != NCL_double) {
      coerce_subset_input_double(z,tmp_z,index_cape,type_z,size_cape,0,NULL,NULL);
    }
    else {
/*
 * Point tmp_z to appropriate location in z.
 */
      tmp_z = &((double*)z)[index_cape];
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
 * Call Fortran routine.
 */
    NGCALLF(dcapecalc3d,DCAPECALC3D)(tmp_p, tmp_t, tmp_q, tmp_z, tmp_zsfc,
                                     tmp_psfc, tmp_cape, tmp_cin, &cmsg, 
                                     &imiy, &imjx, &imkzh, &i3dflag, &iter,
                                     psa_file,errstat,errmsg,strlen(psa_file),
									 ERRLEN);

/* Terminate if there was an error */
	if (errstat != 0) {
		fprintf(stderr, errmsg);
		exit(errstat);
	}

/*
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
  if(type_p != NCL_double) NclFree(tmp_p);
  if(type_t != NCL_double) NclFree(tmp_t);
  if(type_q != NCL_double) NclFree(tmp_q);
  if(type_z != NCL_double) NclFree(tmp_z);
  if(type_zsfc != NCL_double) NclFree(tmp_zsfc);
  if(type_psfc != NCL_double) NclFree(tmp_psfc);
  NclFree(tmp_cape);
  NclFree(tmp_cin);
  NclFree(errmsg);
/*
 * Set up variable to return.
 */
  ret = NclReturnValue(cape,ndims_cape,dsizes_cape,&missing_cape,type_cape,0);
  NclFree(dsizes_cape);
  return(ret);
}

