#include <stdio.h>
#include "wrapper.h"

extern void NGCALLF(wrf_vintrp,WRF_VINTRP)(double *, double *, double *, 
                                           double *, double *, double *, 
                                           double *, double *, double *, 
                                           double *, double *, int *, 
                                           int *, int *, int *, int *, 
                                           int *, int *, int *, double *);

extern void NGCALLF(wrf_monotonic,WRF_MONOTONIC)(double *, double *, double *, 
                                                 double *, int *, double *, 
                                                 int *, int *, int *, int *);


NhlErrorTypes wrf_vintrp_W( void )
{

/*
 * Input variables
 */
/*
 * Argument # 0
 */
  void *field;
  double *tmp_field;
  int       ndims_field;
  ng_size_t dsizes_field[NCL_MAX_DIMENSIONS];
  int has_missing_field;
  NclScalar missing_field, missing_flt_field, missing_dbl_field;
  NclBasicDataTypes type_field;

/*
 * Argument # 1
 */
  void *pres;
  double *tmp_pres;
  int       ndims_pres;
  ng_size_t dsizes_pres[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_pres;

/*
 * Argument # 2
 */
  void *tk;
  double *tmp_tk;
  int       ndims_tk;
  ng_size_t dsizes_tk[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_tk;

/*
 * Argument # 3
 */
  void *qvp;
  double *tmp_qvp;
  int       ndims_qvp;
  ng_size_t dsizes_qvp[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_qvp;

/*
 * Argument # 4
 */
  void *ght;
  double *tmp_ght;
  int       ndims_ght;
  ng_size_t dsizes_ght[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_ght;

/*
 * Argument # 5
 */
  void *ter;
  double *tmp_ter;
  ng_size_t dsizes_ter[2];
  NclBasicDataTypes type_ter;

/*
 * Argument # 6
 */
  void *sfp;
  double *tmp_sfp;
  int       ndims_sfp;
  ng_size_t dsizes_sfp[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_sfp;

/*
 * Argument # 7
 */
  void *smsfp;
  double *tmp_smsfp;
  int       ndims_smsfp;
  ng_size_t dsizes_smsfp[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_smsfp;

/*
 * Argument # 8
 */
  void *vcarray;
  double *tmp_vcarray;
  int       ndims_vcarray;
  ng_size_t dsizes_vcarray[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_vcarray;

/*
 * Argument # 9
 */
  void *intrp_levels;
  double *tmp_intrp_levels;
  ng_size_t dsizes_intrp_levels[1];
  NclBasicDataTypes type_intrp_levels;

/*
 * Arguments # 10-13
 */
  int *icase, *extrap, *vcor, *logp;

/*
 * Return variable
 */
  void *field_out;
  double *tmp_field_out;
  ng_size_t *dsizes_field_out;
  NclScalar missing_field_out;
  NclBasicDataTypes type_field_out;
  
/*
 * Various
 */
  ng_size_t ninlev, nlat, nlon, ninlevlatlon, nlatlon;
  ng_size_t noutlev, noutlevlatlon;
  ng_size_t index_field, index_sfp, index_field_out;
  ng_size_t i, size_leftmost, size_output;
  int ininlev, inlat, inlon, inoutlev, ret;

/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
/*
 * Get argument # 0
 */
  field = (void*)NclGetArgValue(
           0,
           14,
           &ndims_field,
           dsizes_field,
           &missing_field,
           &has_missing_field,
           &type_field,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_field < 3 || ndims_field > 4) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: The field array must be 3D or 4D");
    return(NhlFATAL);
  }

/*
 * Coerce missing value to double if necessary.
 */
  coerce_missing(type_field,has_missing_field,&missing_field,
                 &missing_dbl_field,&missing_flt_field);

  ninlev = dsizes_field[ndims_field-3];
  nlat   = dsizes_field[ndims_field-2];
  nlon   = dsizes_field[ndims_field-1];

/*
 * Test dimension sizes.
 */
  if(ninlev > INT_MAX || nlat > INT_MAX || nlon > INT_MAX) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: one of bottom_top, south_north, or west_east is greater than INT_MAX");
    return(NhlFATAL);
  }
  ininlev = (int) ninlev;
  inlat   = (int) nlat;
  inlon   = (int) nlon;

/*
 * Get argument # 1
 */
  pres = (void*)NclGetArgValue(
           1,
           14,
           &ndims_pres,
           dsizes_pres,
           NULL,
           NULL,
           &type_pres,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_pres != ndims_field) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: The pres and field arrays must have the same dimensionality");
    return(NhlFATAL);
  }
  else {
    for(i = 0; i < ndims_field; i++) {
      if(dsizes_pres[i] != dsizes_field[i]) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: The pres and field arrays must have the same dimensionality");
        return(NhlFATAL);
      }
    }
  }

/*
 * Get argument # 2
 */
  tk = (void*)NclGetArgValue(
           2,
           14,
           &ndims_tk,
           dsizes_tk,
           NULL,
           NULL,
           &type_tk,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_tk != ndims_field) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: The tk and field arrays must have the same dimensionality");
    return(NhlFATAL);
  }
  else {
    for(i = 0; i < ndims_field; i++) {
      if(dsizes_tk[i] != dsizes_field[i]) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: The tk and field arrays must have the same dimensionality");
        return(NhlFATAL);
      }
    }
  }


/*
 * Get argument # 3
 */
  qvp = (void*)NclGetArgValue(
           3,
           14,
           &ndims_qvp,
           dsizes_qvp,
           NULL,
           NULL,
           &type_qvp,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_qvp != ndims_field) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: The qvp and field arrays must have the same dimensionality");
    return(NhlFATAL);
  }
  else {
    for(i = 0; i < ndims_field; i++) {
      if(dsizes_qvp[i] != dsizes_field[i]) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: The qvp and field arrays must have the same dimensionality");
        return(NhlFATAL);
      }
    }
  }

/*
 * Get argument # 4
 */
  ght = (void*)NclGetArgValue(
           4,
           14,
           &ndims_ght,
           dsizes_ght,
           NULL,
           NULL,
           &type_ght,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_ght != ndims_field) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: The ght and field arrays must have the same dimensionality");
    return(NhlFATAL);
  }
  else {
    for(i = 0; i < ndims_field; i++) {
      if(dsizes_ght[i] != dsizes_field[i]) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: The ght and field arrays must have the same dimensionality");
        return(NhlFATAL);
      }
    }
  }

/*
 * Get argument # 5
 */
  ter = (void*)NclGetArgValue(
           5,
           14,
           NULL,
           dsizes_ter,
           NULL,
           NULL,
           &type_ter,
           DONT_CARE);
/*
 * Check dimension sizes for ter.  It can either be 2D, or one fewer
 * dimensions than field.
 */
  if(dsizes_ter[0] != nlat || dsizes_ter[1] != nlon) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: The dimensions of ter must be south_north x west_east");
    return(NhlFATAL);
  }

/*
 * Get argument # 6
 */
  sfp = (void*)NclGetArgValue(
           6,
           14,
           &ndims_sfp,
           dsizes_sfp,
           NULL,
           NULL,
           &type_sfp,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_sfp != (ndims_field-1)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: The sfp array must have the same dimensionality as the field array, minus the level dimension");
    return(NhlFATAL);
  }
  else {
    for(i = 0; i < ndims_field-3; i++) {
      if(dsizes_sfp[i] != dsizes_field[i]) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: The sfp array must have the same dimensionality as the field array, minus the level dimension");
        return(NhlFATAL);
      }
    }
    if(dsizes_sfp[ndims_sfp-2] != nlat || dsizes_sfp[ndims_sfp-1] != nlon) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: The sfp array must have the same dimensionality as the field array, minus the level dimension");
      return(NhlFATAL);
    }
  }


/*
 * Get argument # 7
 */
  smsfp = (void*)NclGetArgValue(
           7,
           14,
           &ndims_smsfp,
           dsizes_smsfp,
           NULL,
           NULL,
           &type_smsfp,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_smsfp != (ndims_field-1)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: The smsfp array must have the same dimensionality as the field array, minus the level dimension");
    return(NhlFATAL);
  }
  else {
    for(i = 0; i < ndims_field-3; i++) {
      if(dsizes_smsfp[i] != dsizes_field[i]) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: The smsfp array must have the same dimensionality as the field array, minus the level dimension");
        return(NhlFATAL);
      }
    }
    if(dsizes_smsfp[ndims_smsfp-2] != nlat || 
       dsizes_smsfp[ndims_smsfp-1] != nlon) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: The smsfp array must have the same dimensionality as the field array, minus the level dimension");
      return(NhlFATAL);
    }
  }


/*
 * Get argument # 8
 */
  vcarray = (void*)NclGetArgValue(
           8,
           14,
           &ndims_vcarray,
           dsizes_vcarray,
           NULL,
           NULL,
           &type_vcarray,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_vcarray != ndims_field) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: The vcarray and field arrays must have the same dimensionality");
    return(NhlFATAL);
  }
  else {
    for(i = 0; i < ndims_field; i++) {
      if(dsizes_vcarray[i] != dsizes_field[i]) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: The vcarray and field arrays must have the same dimensionality");
        return(NhlFATAL);
      }
    }
  }

/*
 * Get argument # 9
 */
  intrp_levels = (void*)NclGetArgValue(
           9,
           14,
           NULL,
           dsizes_intrp_levels,
           NULL,
           NULL,
           &type_intrp_levels,
           DONT_CARE);

  noutlev = dsizes_intrp_levels[0];

/*
 * Test dimension sizes.
 */
  if(noutlev > INT_MAX) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: the # of output levels is greater than INT_MAX");
    return(NhlFATAL);
  }
  inoutlev = (int) noutlev;


/*
 * Get argument # 10
 */
  icase = (int*)NclGetArgValue(
           10,
           14,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);
/*
 * Get argument # 11
 */
  extrap = (int*)NclGetArgValue(
           11,
           14,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);
/*
 * Get argument # 12
 */
  vcor = (int*)NclGetArgValue(
           12,
           14,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);
/*
 * Get argument # 13
 */
  logp = (int*)NclGetArgValue(
           13,
           14,
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
  for(i = 0; i < ndims_field-3; i++) size_leftmost *= dsizes_field[i];

/* 
 * Allocate space for coercing input arrays.  If any of the input
 * is already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 */
/*
 * Allocate space for tmp_field.
 */
  nlatlon       = nlat * nlon;
  ninlevlatlon  = ninlev * nlatlon;
  noutlevlatlon = noutlev * nlatlon;

  if(type_field != NCL_double) {
    tmp_field = (double *)calloc(ninlevlatlon,sizeof(double));
    if(tmp_field == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: Unable to allocate memory for coercing field array to double");
      return(NhlFATAL);
    }
  }

/*
 * Allocate space for tmp_pres.
 */
  if(type_pres != NCL_double) {
    tmp_pres = (double *)calloc(ninlevlatlon,sizeof(double));
    if(tmp_pres == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: Unable to allocate memory for coercing pressure array to double");
      return(NhlFATAL);
    }
  }

/*
 * Allocate space for tmp_tk.
 */
  if(type_tk != NCL_double) {
    tmp_tk = (double *)calloc(ninlevlatlon,sizeof(double));
    if(tmp_tk == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: Unable to allocate memory for coercing tk array to double");
      return(NhlFATAL);
    }
  }

/*
 * Allocate space for tmp_qvp.
 */
  if(type_qvp != NCL_double) {
    tmp_qvp = (double *)calloc(ninlevlatlon,sizeof(double));
    if(tmp_qvp == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: Unable to allocate memory for coercing qvp array to double");
      return(NhlFATAL);
    }
  }

/*
 * Allocate space for tmp_ght.
 */
  if(type_ght != NCL_double) {
    tmp_ght = (double *)calloc(ninlevlatlon,sizeof(double));
    if(tmp_ght == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: Unable to allocate memory for coercing ght array to double");
      return(NhlFATAL);
    }
  }

/*
 * Coerce ter to double, if necessary.
 */
  tmp_ter = coerce_input_double(ter,type_ter,nlatlon,0,NULL,NULL);
  if(tmp_ter == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: Unable to coerce ter to double precision");
    return(NhlFATAL);
  }

/*
 * Allocate space for tmp_sfp.
 */
  if(type_sfp != NCL_double) {
    tmp_sfp = (double *)calloc(nlatlon,sizeof(double));
    if(tmp_sfp == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: Unable to allocate memory for coercing sfp array to double");
      return(NhlFATAL);
    }
  }

/*
 * Allocate space for tmp_smsfp.
 */
  if(type_smsfp != NCL_double) {
    tmp_smsfp = (double *)calloc(nlatlon,sizeof(double));
    if(tmp_smsfp == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: Unable to allocate memory for coercing smsfp array to double");
      return(NhlFATAL);
    }
  }

/*
 * Allocate space for tmp_vcarray.
 */
  if(type_vcarray != NCL_double) {
    tmp_vcarray = (double *)calloc(ninlevlatlon,sizeof(double));
    if(tmp_vcarray == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: Unable to allocate memory for coercing vcarray to double");
      return(NhlFATAL);
    }
  }
/*
 * Coerce intrp_levels to double, if necessary.
 */
  tmp_intrp_levels = coerce_input_double(intrp_levels,type_intrp_levels,
                                         noutlev,0,NULL,NULL);
  if(tmp_intrp_levels == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: Unable to allocate memory for coercing interp_levels array to double");
    return(NhlFATAL);
  }

/*
 * The output type defaults to float, unless one or more input 
 * arrays are double.
 */
  if(type_field   == NCL_double || type_pres  == NCL_double || 
     type_tk      == NCL_double || type_qvp   == NCL_double || 
     type_ght     == NCL_double || type_ter   == NCL_double || 
     type_sfp     == NCL_double || type_smsfp == NCL_double || 
     type_vcarray == NCL_double) {
    type_field_out = NCL_double;
  }
  else {
    type_field_out = NCL_float;
  }

/* 
 * Allocate space for output array and set a missing value.
 * Note: a missing value is returned even if the input doesn't
 * have a missing value, because the output may contain missing 
 * values after interpolation is done.
 */
  size_output = size_leftmost * noutlevlatlon;
  if(type_field_out != NCL_double) {
    field_out = (void *)calloc(size_output, sizeof(float));
    tmp_field_out = (double *)calloc(noutlevlatlon,sizeof(double));
    if(field_out == NULL || tmp_field_out == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: Unable to allocate memory for temporary output array");
      return(NhlFATAL);
    }
    missing_field_out = missing_flt_field;
  }
  else {
    field_out = (void *)calloc(size_output, sizeof(double));
    if(field_out == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
    missing_field_out = missing_dbl_field;
  }

/* 
 * Allocate space for output dimension sizes and set them.
 */
  dsizes_field_out = (ng_size_t*)calloc(ndims_field,sizeof(ng_size_t));  
  if( dsizes_field_out == NULL ) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: Unable to allocate memory for holding dimension sizes");
    return(NhlFATAL);
  }
  for(i = 0; i < ndims_field-3; i++) dsizes_field_out[i] = dsizes_field[i];
  dsizes_field_out[ndims_field-3] = noutlev;
  dsizes_field_out[ndims_field-2] = nlat;
  dsizes_field_out[ndims_field-1] = nlon;

/*
 * Loop across leftmost dimensions and call the Fortran routine for each
 * subsection of the input arrays.
 */
  index_field = index_sfp = index_field_out = 0;

  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of field (tmp_field) to double if necessary.
 */
    if(type_field != NCL_double) {
      coerce_subset_input_double(field,tmp_field,index_field,type_field,
                                 ninlevlatlon,0,NULL,NULL);
    }
    else {
      tmp_field = &((double*)field)[index_field];
    }

/*
 * Coerce subsection of pres (tmp_pres) to double if necessary.
 */
    if(type_pres != NCL_double) {
      coerce_subset_input_double(pres,tmp_pres,index_field,type_pres,ninlevlatlon,0,NULL,NULL);
    }
    else {
      tmp_pres = &((double*)pres)[index_field];
    }

/*
 * Coerce subsection of tk (tmp_tk) to double if necessary.
 */
    if(type_tk != NCL_double) {
      coerce_subset_input_double(tk,tmp_tk,index_field,type_tk,
                                 ninlevlatlon,0,NULL,NULL);
    }
    else {
      tmp_tk = &((double*)tk)[index_field];
    }

/*
 * Coerce subsection of qvp (tmp_qvp) to double if necessary.
 */
    if(type_qvp != NCL_double) {
      coerce_subset_input_double(qvp,tmp_qvp,index_field,type_qvp,
                                 ninlevlatlon,0,NULL,NULL);
    }
    else {
      tmp_qvp = &((double*)qvp)[index_field];
    }

/*
 * Coerce subsection of ght (tmp_ght) to double if necessary.
 */
    if(type_ght != NCL_double) {
      coerce_subset_input_double(ght,tmp_ght,index_field,type_ght,
                                 ninlevlatlon,0,NULL,NULL);
    }
    else {
      tmp_ght = &((double*)ght)[index_field];
    }

/*
 * Coerce subsection of sfp (tmp_sfp) to double if necessary.
 */
    if(type_sfp != NCL_double) {
      coerce_subset_input_double(sfp,tmp_sfp,index_sfp,type_sfp,
                                 nlatlon,0,NULL,NULL);
    }
    else {
      tmp_sfp = &((double*)sfp)[index_sfp];
    }

/*
 * Coerce subsection of smsfp (tmp_smsfp) to double if necessary.
 */
    if(type_smsfp != NCL_double) {
      coerce_subset_input_double(smsfp,tmp_smsfp,index_sfp,type_smsfp,
                                 nlatlon,0,NULL,NULL);
    }
    else {
      tmp_smsfp = &((double*)smsfp)[index_sfp];
    }

/*
 * Coerce subsection of vcarray (tmp_vcarray) to double if necessary.
 */
    if(type_vcarray != NCL_double) {
      coerce_subset_input_double(vcarray,tmp_vcarray,index_field,type_vcarray,
                                 ninlevlatlon,0,NULL,NULL);
    }
    else {
      tmp_vcarray = &((double*)vcarray)[index_field];
    }


/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_field_out == NCL_double) {
      tmp_field_out = &((double*)field_out)[index_field_out];
    }
    
/*
 * Call the Fortran routine.
 */
    NGCALLF(wrf_vintrp,WRF_VINTRP)(tmp_field, tmp_field_out, tmp_pres,
                                   tmp_tk, tmp_qvp, tmp_ght, tmp_ter, 
                                   tmp_sfp, tmp_smsfp, tmp_vcarray, 
                                   tmp_intrp_levels, &inoutlev, icase, 
                                   &inlon, &inlat, &ininlev, extrap, vcor, 
                                   logp, &missing_dbl_field.doubleval);

/*
 * Coerce output back to float if necessary.
 */
    if(type_field_out == NCL_float) {
      coerce_output_float_only(field_out,tmp_field_out,noutlevlatlon,
                               index_field_out);
    }
    index_field     += ninlevlatlon;
    index_field_out += noutlevlatlon;
    index_sfp       += nlatlon;
  }

/*
 * Free unneeded memory.
 */

  if(type_field        != NCL_double) NclFree(tmp_field);
  if(type_pres         != NCL_double) NclFree(tmp_pres);
  if(type_tk           != NCL_double) NclFree(tmp_tk);
  if(type_qvp          != NCL_double) NclFree(tmp_qvp);
  if(type_ght          != NCL_double) NclFree(tmp_ght);
  if(type_ter          != NCL_double) NclFree(tmp_ter);
  if(type_sfp          != NCL_double) NclFree(tmp_sfp);
  if(type_smsfp        != NCL_double) NclFree(tmp_smsfp);
  if(type_vcarray      != NCL_double) NclFree(tmp_vcarray);
  if(type_intrp_levels != NCL_double) NclFree(tmp_intrp_levels);
  if(type_field_out    != NCL_double) NclFree(tmp_field_out);

/*
 * Return value back to NCL script.
 */
  ret = NclReturnValue(field_out,ndims_field,dsizes_field_out,
                       &missing_field_out,type_field_out,0);

  NclFree(dsizes_field_out);
  return(ret);
}

NhlErrorTypes wrf_monotonic_W( void )
{

/*
 * Input variables
 */
/*
 * Argument # 0
 */
  void *x;
  double *tmp_x;
  int       ndims_x;
  ng_size_t dsizes_x[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_x;

/*
 * Argument # 1
 */
  void *pres;
  double *tmp_pres;
  int       ndims_pres;
  ng_size_t dsizes_pres[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_pres;

/*
 * Argument # 2
 */
  void *cor;
  double *tmp_cor;
  ng_size_t dsizes_cor[2];
  NclBasicDataTypes type_cor;

/*
 * Argument # 3
 */
  int *idir;
/*
 * Argument # 4
 */
  void *delta;
  double *tmp_delta;
  NclBasicDataTypes type_delta;

/*
 * Argument # 5
 */
  int *icorsw;
/*
 * Return variable
 */
  void *xout;
  double *tmp_xout;
  NclBasicDataTypes type_xout;


/*
 * Various
 */
  ng_size_t nlev, nlat, nlon, nlevlatlon, nlatlon;
  ng_size_t i, index_x, size_leftmost, size_output;
  int inlev, inlat, inlon, ret;
/*
 * Retrieve parameters.
 *
 * Note any of the pointer parameters can be set to NULL, which
 * implies you don't care about its value.
 */
/*
 * Get argument # 0
 */
  x = (void*)NclGetArgValue(
           0,
           6,
           &ndims_x,
           dsizes_x,
           NULL,
           NULL,
           &type_x,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_x < 3) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_monotonic: The x array must have at least three dimensions");
    return(NhlFATAL);
  }
  nlev = dsizes_x[ndims_x-3];
  nlat = dsizes_x[ndims_x-2];
  nlon = dsizes_x[ndims_x-1];

/*
 * Test dimension sizes.
 */
  if(nlev > INT_MAX || nlat > INT_MAX || nlon > INT_MAX) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_monotonic: one of bottom_top, south_north, or west_east is greater than INT_MAX");
    return(NhlFATAL);
  }
  inlev = (int) nlev;
  inlat = (int) nlat;
  inlon = (int) nlon;


/*
 * Get argument # 1
 */
  pres = (void*)NclGetArgValue(
           1,
           6,
           &ndims_pres,
           dsizes_pres,
           NULL,
           NULL,
           &type_pres,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_pres != ndims_x) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_monotonic: The pres and x arrays must have the same dimensionality");
    return(NhlFATAL);
  }
  else {
    for(i = 0; i < ndims_x; i++) {
      if(dsizes_pres[i] != dsizes_x[i]) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_monotonic: The pres and x arrays must have the same dimensionality");
        return(NhlFATAL);
      }
    }
  }

/*
 * Get argument # 2
 */
  cor = (void*)NclGetArgValue(
           2,
           6,
           NULL,
           dsizes_cor,
           NULL,
           NULL,
           &type_cor,
           DONT_CARE);

  if(dsizes_cor[0] != nlat || dsizes_cor[1] != nlon) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_vintrp: The dimensions of cor must be south_north x west_east");
    return(NhlFATAL);
  }

/*
 * Get argument # 3
 */
  idir = (int*)NclGetArgValue(
           3,
           6,
           NULL,
           NULL,
           NULL,
           NULL,
           NULL,
           DONT_CARE);
/*
 * Get argument # 4
 */
  delta = (void*)NclGetArgValue(
           4,
           6,
           NULL,
           NULL,
           NULL,
           NULL,
           &type_delta,
           DONT_CARE);
/*
 * Get argument # 5
 */
  icorsw = (int*)NclGetArgValue(
           5,
           6,
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
  for(i = 0; i < ndims_x-3; i++) size_leftmost *= dsizes_x[i];

/* 
 * Allocate space for coercing input arrays.  If any of the input
 * is already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 */

/*
 * Allocate space for tmp_x.
 */
  nlatlon     = nlat * nlon;
  nlevlatlon = nlev * nlatlon;

  if(type_x != NCL_double) {
    type_xout = NCL_float;
    tmp_x = (double *)calloc(nlevlatlon,sizeof(double));
    if(tmp_x == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_monotonic: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
  else {
    type_xout = NCL_double;
  }
/*
 * Allocate space for tmp_pres.
 */
  if(type_pres != NCL_double) {
    tmp_pres = (double *)calloc(nlevlatlon,sizeof(double));
    if(tmp_pres == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_monotonic: Unable to allocate memory for coercing input array to double");
      return(NhlFATAL);
    }
  }
/*
 * Coerce cor to double, if necessary.
 */
  tmp_cor = coerce_input_double(cor,type_cor,nlatlon,0,NULL,NULL);
  if(tmp_cor == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_monotonic: Unable to coerce cor to double precision");
    return(NhlFATAL);
  }

/*
 * Allocate space for tmp_delta.
 */
  tmp_delta = coerce_input_double(delta,type_delta,1,0,NULL,NULL);
  if(tmp_delta == NULL) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_monotonic: Unable to allocate memory for coercing input array to double");
    return(NhlFATAL);
  }

/*
 * Calculate size of output array, which is same as x.
 */
  size_output = size_leftmost * nlevlatlon;

/* 
 * Allocate space for output array.
 */
  if(type_xout != NCL_double) {
    xout = (void *)calloc(size_output, sizeof(float));
    tmp_xout = (double *)calloc(nlevlatlon,sizeof(double));
    if(xout == NULL || tmp_xout == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_monotonic: Unable to allocate memory for temporary output array");
      return(NhlFATAL);
    }
  }
  else {
    xout = (void *)calloc(size_output, sizeof(double));
    if(xout == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_monotonic: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }

/*
 * Loop across leftmost dimensions and call the Fortran routine for each
 * subsection of the input arrays.
 */
  index_x = 0;
  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of x (tmp_x) to double if necessary.
 */
    if(type_x != NCL_double) {
      coerce_subset_input_double(x,tmp_x,index_x,type_x,nlevlatlon,0,NULL,NULL);
    }
    else {
      tmp_x = &((double*)x)[index_x];
    }

/*
 * Coerce subsection of pres (tmp_pres) to double if necessary.
 */
    if(type_pres != NCL_double) {
      coerce_subset_input_double(pres,tmp_pres,index_x,type_pres,nlevlatlon,0,NULL,NULL);
    }
    else {
      tmp_pres = &((double*)pres)[index_x];
    }

/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_xout == NCL_double) tmp_xout = &((double*)xout)[index_x];

/*
 * Call the Fortran routine.
 */
    NGCALLF(wrf_monotonic,WRF_MONOTONIC)(tmp_xout, tmp_x, tmp_pres, 
                                         tmp_cor, idir, tmp_delta, &inlon,
                                         &inlat, &inlev, icorsw);
/*
 * Coerce output back to float if necessary.
 */
    if(type_xout == NCL_float) {
      coerce_output_float_only(xout,tmp_xout,nlevlatlon,index_x);
    }
    index_x += nlevlatlon;
  }

/*
 * Free unneeded memory.
 */
  if(type_x     != NCL_double) NclFree(tmp_x);
  if(type_pres  != NCL_double) NclFree(tmp_pres);
  if(type_cor   != NCL_double) NclFree(tmp_cor);
  if(type_delta != NCL_double) NclFree(tmp_delta);
  if(type_xout  != NCL_double) NclFree(tmp_xout);

/*
 * Return value back to NCL script.
 */
  ret = NclReturnValue(xout,ndims_x,dsizes_x,NULL,type_xout,0);
  return(ret);
}
