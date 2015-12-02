#include <stdio.h>
#include "wrapper.h"

extern void NGCALLF(wrfcttcalc,WRFCTTCALC)(double *, double *, double *, 
                                           double *, double *, double *, 
                                           double *, double *, int *, 
                                           int *, int *, int *);

extern NclDimRec *get_wrf_dim_info(int,int,int,ng_size_t*);


NhlErrorTypes wrf_ctt_W( void )
{

/*
 * Input variables
 */
/*
 * Argument # 0
 */
  void *pres;
  double *tmp_pres;
  int       ndims_pres;
  ng_size_t dsizes_pres[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_pres;

/*
 * Argument # 1
 */
  void *tk;
  double *tmp_tk;
  int       ndims_tk;
  ng_size_t dsizes_tk[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_tk;

/*
 * Argument # 2
 */
  void *qci;
  double *tmp_qci;
  int       ndims_qci;
  ng_size_t dsizes_qci[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_qci;

/*
 * Argument # 3
 */
  void *qcw;
  double *tmp_qcw;
  int       ndims_qcw;
  ng_size_t dsizes_qcw[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_qcw;

/*
 * Argument # 4
 */
  void *qvp;
  double *tmp_qvp;
  int       ndims_qvp;
  ng_size_t dsizes_qvp[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_qvp;

/*
 * Argument # 5
 */
  void *ght;
  double *tmp_ght;
  int       ndims_ght;
  ng_size_t dsizes_ght[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_ght;

/*
 * Argument # 6
 */
  void *ter;
  double *tmp_ter;
  int       ndims_ter;
  ng_size_t dsizes_ter[NCL_MAX_DIMENSIONS];
  NclBasicDataTypes type_ter;

/*
 * Arguments # 7
 */
  int *haveqci;

/*
 * Variable for getting/setting dimension name info.
 */
  NclDimRec *dim_info      = NULL;
  NclDimRec *dim_info_ght = NULL;

/*
 * Return variable and attributes
 */
  void *ctt;
  NclQuark *description, *units;
  char *cdescription, *cunits;
  double *tmp_ctt;
  int       ndims_ctt;
  ng_size_t *dsizes_ctt;
  NclBasicDataTypes type_ctt;
  NclObjClass type_obj_ctt;
  
/*
 * Various
 */
  ng_size_t nlev, nlat, nlon, nlevlatlon, nlatlon;
  ng_size_t index_pres, index_ter, index_ctt;
  ng_size_t i, size_leftmost, size_output;
  int inlev, inlat, inlon;

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

/*
 * Get argument # 1
 */
  pres = (void*)NclGetArgValue(
           0,
           8,
           &ndims_pres,
           dsizes_pres,
           NULL,
           NULL,
           &type_pres,
           DONT_CARE);

  if(ndims_pres < 3 || ndims_pres > 4) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: The pres array must be 3D or 4D");
    return(NhlFATAL);
  }

  nlev = dsizes_pres[ndims_pres-3];
  nlat = dsizes_pres[ndims_pres-2];
  nlon = dsizes_pres[ndims_pres-1];

/*
 * Test dimension sizes.
 */
  if(nlev > INT_MAX || nlat > INT_MAX || nlon > INT_MAX) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: one of bottom_top, south_north, or west_east is greater than INT_MAX");
    return(NhlFATAL);
  }
  inlev = (int) nlev;
  inlat = (int) nlat;
  inlon = (int) nlon;

/*
 * Get argument # 1
 */
  tk = (void*)NclGetArgValue(
           1,
           8,
           &ndims_tk,
           dsizes_tk,
           NULL,
           NULL,
           &type_tk,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_tk != ndims_pres) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: The tk and pres arrays must have the same dimensionality");
    return(NhlFATAL);
  }
  else {
    for(i = 0; i < ndims_pres; i++) {
      if(dsizes_tk[i] != dsizes_pres[i]) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: The tk and pres arrays must have the same dimensionality");
        return(NhlFATAL);
      }
    }
  }


/*
 * Get argument # 2
 */
  qci = (void*)NclGetArgValue(
           2,
           8,
           &ndims_qci,
           dsizes_qci,
           NULL,
           NULL,
           &type_qci,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_qci != ndims_pres) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: The qci and pres arrays must have the same dimensionality");
    return(NhlFATAL);
  }
  else {
    for(i = 0; i < ndims_pres; i++) {
      if(dsizes_qci[i] != dsizes_pres[i]) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: The qci and pres arrays must have the same dimensionality");
        return(NhlFATAL);
      }
    }
  }

/*
 * Get argument # 3
 */
  qcw = (void*)NclGetArgValue(
           3,
           8,
           &ndims_qcw,
           dsizes_qcw,
           NULL,
           NULL,
           &type_qcw,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_qcw != ndims_pres) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: The qcw and pres arrays must have the same dimensionality");
    return(NhlFATAL);
  }
  else {
    for(i = 0; i < ndims_pres; i++) {
      if(dsizes_qcw[i] != dsizes_pres[i]) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: The qcw and pres arrays must have the same dimensionality");
        return(NhlFATAL);
      }
    }
  }

/*
 * Get argument # 4
 */
  qvp = (void*)NclGetArgValue(
           4,
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
  if(ndims_qvp != ndims_pres) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: The qvp and pres arrays must have the same dimensionality");
    return(NhlFATAL);
  }
  else {
    for(i = 0; i < ndims_pres; i++) {
      if(dsizes_qvp[i] != dsizes_pres[i]) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: The qvp and pres arrays must have the same dimensionality");
        return(NhlFATAL);
      }
    }
  }

/*
 * Get argument # 5
 */
  ght = (void*)NclGetArgValue(
           5,
           8,
           &ndims_ght,
           dsizes_ght,
           NULL,
           NULL,
           &type_ght,
           DONT_CARE);

/*
 * Check dimension sizes.
 */
  if(ndims_ght != ndims_pres) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: The ght and pres arrays must have the same dimensionality");
    return(NhlFATAL);
  }
  else {
    for(i = 0; i < ndims_pres; i++) {
      if(dsizes_ght[i] != dsizes_pres[i]) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: The ght and pres arrays must have the same dimensionality");
        return(NhlFATAL);
      }
    }
  }

/*
 * Get argument # 6
 */
  ter = (void*)NclGetArgValue(
           6,
           8,
           &ndims_ter,
           dsizes_ter,
           NULL,
           NULL,
           &type_ter,
           DONT_CARE);

/*
 * Check dimension sizes for ter.  It can either be 2D, or one fewer
 * dimensions than pres.
 */
  if(ndims_ter != 2 && ndims_ter != (ndims_pres-1)) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: ter must either be a 2D array dimensioned south_north x west_east or it must have the same dimensionality as the pres array, minus the level dimension");
    return(NhlFATAL);
  }

  if(ndims_ter == 2) {
    if(dsizes_ter[0] != nlat || dsizes_ter[1] != nlon) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: The dimensions of ter must be south_north x west_east");
      return(NhlFATAL);
    }
  }
  else {
    for(i = 0; i < ndims_pres-3; i++) {
      if(dsizes_ter[i] != dsizes_pres[i]) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: ter must either be a 2D array dimensioned south_north x west_east or it must have the same dimensionality as the pres array, minus the level dimension");
        return(NhlFATAL);
      }
    }
  }


/*
 * Get argument # 7
 */
  haveqci = (int*)NclGetArgValue(
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
  for(i = 0; i < ndims_pres-3; i++) size_leftmost *= dsizes_pres[i];

/* 
 * Allocate space for coercing input arrays.  If any of the input
 * is already double, then we don't need to allocate space for
 * temporary arrays, because we'll just change the pointer into
 * the void array appropriately.
 */
/*
 * Allocate space for tmp_pres.
 */
  nlatlon    = nlat * nlon;
  nlevlatlon = nlev * nlatlon;

  if(type_pres != NCL_double) {
    tmp_pres = (double *)calloc(nlevlatlon,sizeof(double));
    if(tmp_pres == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: Unable to allocate memory for coercing pressure array to double");
      return(NhlFATAL);
    }
  }

/*
 * Allocate space for tmp_tk.
 */
  if(type_tk != NCL_double) {
    tmp_tk = (double *)calloc(nlevlatlon,sizeof(double));
    if(tmp_tk == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: Unable to allocate memory for coercing tk array to double");
      return(NhlFATAL);
    }
  }

/*
 * Allocate space for tmp_qci.
 */
  if(type_qci != NCL_double) {
    tmp_qci = (double *)calloc(nlevlatlon,sizeof(double));
    if(tmp_qci == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: Unable to allocate memory for coercing qci array to double");
      return(NhlFATAL);
    }
  }

/*
 * Allocate space for tmp_qcw.
 */
  if(type_qcw != NCL_double) {
    tmp_qcw = (double *)calloc(nlevlatlon,sizeof(double));
    if(tmp_qcw == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: Unable to allocate memory for coercing qcw array to double");
      return(NhlFATAL);
    }
  }

/*
 * Allocate space for tmp_qvp.
 */
  if(type_qvp != NCL_double) {
    tmp_qvp = (double *)calloc(nlevlatlon,sizeof(double));
    if(tmp_qvp == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: Unable to allocate memory for coercing qvp array to double");
      return(NhlFATAL);
    }
  }

/*
 * Allocate space for tmp_ght.
 */
  if(type_ght != NCL_double) {
    tmp_ght = (double *)calloc(nlevlatlon,sizeof(double));
    if(tmp_ght == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: Unable to allocate memory for coercing ght array to double");
      return(NhlFATAL);
    }
  }

/*
 * Coerce ter to double, if necessary.
 */
  if(ndims_ter == 2) {
    tmp_ter = coerce_input_double(ter,type_ter,nlatlon,0,NULL,NULL);
    if(tmp_ter == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: Unable to allocate memory for coercing ter array to double");
      return(NhlFATAL);
    }
  }
  else {
/*
 * Allocate space for tmp_ter.
 */
    if(type_ter != NCL_double) {
      tmp_ter = (double *)calloc(nlatlon,sizeof(double));
      if(tmp_ter == NULL) {
        NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: Unable to allocate memory for coercing ter array to double");
        return(NhlFATAL);
      }
    }
  }


/*
 * The output type defaults to float, unless one or more input 
 * arrays are double.
 */
  if(type_pres == NCL_double || type_tk  == NCL_double || 
     type_qci  == NCL_double || type_qcw == NCL_double || 
     type_qvp  == NCL_double || type_ght == NCL_double || 
     type_ter  == NCL_double) {
    type_ctt     = NCL_double;
    type_obj_ctt = nclTypedoubleClass;
  }
  else {
    type_ctt     = NCL_float;
    type_obj_ctt = nclTypefloatClass;
  }

/* 
 * Allocate space for output array.
 */
  size_output = size_leftmost * nlatlon;
  if(type_ctt != NCL_double) {
    ctt = (void *)calloc(size_output, sizeof(float));
    tmp_ctt = (double *)calloc(nlatlon,sizeof(double));
    if(ctt == NULL || tmp_ctt == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: Unable to allocate memory for temporary output array");
      return(NhlFATAL);
    }
  }
  else {
    ctt = (void *)calloc(size_output, sizeof(double));
    if(ctt == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: Unable to allocate memory for output array");
      return(NhlFATAL);
    }
  }

/* 
 * Allocate space for output dimension sizes and set them.
 */
  ndims_ctt  = ndims_pres-1;
  dsizes_ctt = (ng_size_t*)calloc(ndims_ctt,sizeof(ng_size_t));  
  if( dsizes_ctt == NULL ) {
    NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: Unable to allocate memory for holding dimension sizes");
    return(NhlFATAL);
  }
  for(i = 0; i < ndims_ctt-2; i++) dsizes_ctt[i] = dsizes_pres[i];
  dsizes_ctt[ndims_ctt-2] = nlat;
  dsizes_ctt[ndims_ctt-1] = nlon;

/*
 * Get dimension info to see if we have named dimensions.
 * Using "ght" here, because it is more likely than "pres"
 * to have metadata attached to it. 
 * 
 * This will be used for return variable.
 */
  dim_info_ght = get_wrf_dim_info(5,8,ndims_ght,dsizes_ght);
  if(dim_info_ght != NULL) {
    dim_info = malloc(sizeof(NclDimRec)*ndims_ctt);
    if(dim_info == NULL) {
      NhlPError(NhlFATAL,NhlEUNKNOWN,"wrf_ctt: Unable to allocate memory for holding dimension information");
      return(NhlFATAL);
    }
    for(i = 0; i < ndims_ght-3; i++) {
      dim_info[i] = dim_info_ght[i];
    }
    dim_info[ndims_ctt-1] = dim_info_ght[ndims_ght-1];
    dim_info[ndims_ctt-2] = dim_info_ght[ndims_ght-2];
  }

/*
 * Loop across leftmost dimensions and call the Fortran routine for each
 * subsection of the input arrays.
 */
  index_pres = index_ter = index_ctt = 0;

  for(i = 0; i < size_leftmost; i++) {
/*
 * Coerce subsection of pres (tmp_pres) to double if necessary.
 */
    if(type_pres != NCL_double) {
      coerce_subset_input_double(pres,tmp_pres,index_pres,
                                 type_pres,nlevlatlon,0,NULL,NULL);
    }
    else {
      tmp_pres = &((double*)pres)[index_pres];
    }

/*
 * Coerce subsection of tk (tmp_tk) to double if necessary.
 */
    if(type_tk != NCL_double) {
      coerce_subset_input_double(tk,tmp_tk,index_pres,type_tk,
                                 nlevlatlon,0,NULL,NULL);
    }
    else {
      tmp_tk = &((double*)tk)[index_pres];
    }

/*
 * Coerce subsection of qci (tmp_qci) to double if necessary.
 */
    if(type_qci != NCL_double) {
      coerce_subset_input_double(qci,tmp_qci,index_pres,type_qci,
                                 nlevlatlon,0,NULL,NULL);
    }
    else {
      tmp_qci = &((double*)qci)[index_pres];
    }

/*
 * Coerce subsection of qcw (tmp_qcw) to double if necessary.
 */
    if(type_qcw != NCL_double) {
      coerce_subset_input_double(qcw,tmp_qcw,index_pres,type_qcw,
                                 nlevlatlon,0,NULL,NULL);
    }
    else {
      tmp_qcw = &((double*)qcw)[index_pres];
    }

/*
 * Coerce subsection of qvp (tmp_qvp) to double if necessary.
 */
    if(type_qvp != NCL_double) {
      coerce_subset_input_double(qvp,tmp_qvp,index_pres,type_qvp,
                                 nlevlatlon,0,NULL,NULL);
    }
    else {
      tmp_qvp = &((double*)qvp)[index_pres];
    }

/*
 * Coerce subsection of ght (tmp_ght) to double if necessary.
 */
    if(type_ght != NCL_double) {
      coerce_subset_input_double(ght,tmp_ght,index_pres,type_ght,
                                 nlevlatlon,0,NULL,NULL);
    }
    else {
      tmp_ght = &((double*)ght)[index_pres];
    }

/*
 * Coerce subsection of ter (tmp_ter) to double if necessary.
 */
    if(ndims_ter != 2) {
      if(type_ter != NCL_double) {
        coerce_subset_input_double(ter,tmp_ter,index_ter,type_ter,
                                   nlatlon,0,NULL,NULL);
      }
      else {
        tmp_ter = &((double*)ter)[index_ter];
      }
    }

/*
 * Point temporary output array to void output array if appropriate.
 */
    if(type_ctt == NCL_double) {
      tmp_ctt = &((double*)ctt)[index_ctt];
    }
    
/*
 * Call the Fortran routine.
 */
    NGCALLF(wrfcttcalc,WRFCTTCALC)(tmp_pres, tmp_tk, tmp_qci, tmp_qcw,
                                   tmp_qvp, tmp_ght, tmp_ter, tmp_ctt,
                                   haveqci,&inlev, &inlat, &inlon);

/*
 * Coerce output back to float if necessary.
 */
    if(type_ctt == NCL_float) {
      coerce_output_float_only(ctt,tmp_ctt,nlatlon,
                               index_ctt);
    }
    index_pres  += nlevlatlon;
    index_ctt  += nlatlon;
    if(ndims_ter != 2) { 
      index_ter += nlatlon;
    }
  }

/*
 * Free unneeded memory.
 */
  if(type_pres != NCL_double) NclFree(tmp_pres);
  if(type_tk   != NCL_double) NclFree(tmp_tk);
  if(type_qci  != NCL_double) NclFree(tmp_qci);
  if(type_qcw  != NCL_double) NclFree(tmp_qcw);
  if(type_qvp  != NCL_double) NclFree(tmp_qvp);
  if(type_ght  != NCL_double) NclFree(tmp_ght);
  if(type_ter  != NCL_double) NclFree(tmp_ter);
  if(type_ctt  != NCL_double) NclFree(tmp_ctt);

/*
 * Set up some attributes ("description" and "units") to return.
 */
  cdescription = (char *)calloc(22,sizeof(char));
  cunits       = (char *)calloc(2,sizeof(char));
  strcpy(cdescription,"Cloud Top Temperature");
  strcpy(cunits,"K");
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
                            (void*)ctt,
                            NULL,
                            ndims_ctt,
                            dsizes_ctt,
                            TEMPORARY,
                            NULL,
                            type_obj_ctt
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

  if(dim_info != NULL) NclFree(dim_info);
  NclFree(dim_info_ght);

/*
 * Return output grid and attributes to NCL.
 */
  return_data.kind = NclStk_VAR;
  return_data.u.data_var = tmp_var;
  _NclPlaceReturn(return_data);
  return(NhlNOERROR);
}
