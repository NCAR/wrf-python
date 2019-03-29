.. _internals:

WRF-Python Internals
========================================

WRF-Python is a collection of diagnostic and interpolation routines for 
WRF-ARW data. The API is kept to a minimal set of functions, since we've found
this to be the easiest to teach to new programmers, students, and scientists. 
Future plans include adopting the Pangeo xarray/dask model, along with an 
object oriented API, but this is not currently supported as of this user 
guide.

A typical use case for a WRF-Python user is to:

1) Open a WRF data file (or sequence of files) using NetCDF4-python or PyNIO.
2) Compute a WRF diagnostic using :meth:`wrf.getvar`.
3) Perform any additional computations using methods outside of WRF-Python.
4) Create a plot of the output using matplotlib (basemap or cartopy) or 
   PyNGL.
   
The purpose of this guide is to explain the internals of item (2) so that 
users can help contribute or support the computational diagnostic 
routines.


Overview of a :meth:`wrf.getvar` Diagnostic Computation
---------------------------------------------------------------

A diagnostic computed using the :meth:`wrf.getvar` function consists of the 
following steps:

1) Call the appropriate 'getter' function based on the specified diagnostic 
   label. This step occurs in the :meth:`wrf.getvar` routine in routines.py. 
2) Extract the required variables from the NetCDF data file (or files).
3) Compute the diagnostic using a wrapped Fortran, C, or Python routine.
4) Convert to the desired units (if applicable).
5) Set the metadata (if desired) and return the result as an 
   :class:`xarray.DataArray`, or return a :class:`numpy.ndarray` if no 
   metadata is required.
   
In the source directory, the :meth:`wrf.getvar` 'getter' routines have a 
"\g\_" prefix for the naming convention (the "g" is for "get"). 

The unit conversion is handled by a :mod:`wrapt` decorator that can be found 
in *decorators.py*. The setting of the metadata is handled using a :mod:`wrapt` 
decorator, which can be found in the *metadecorators.py* file.


Overview of Compiled Computational Routines
---------------------------------------------------------

Currently, the compiled computational routines are written in Fortran 
90 and exposed the Python using f2py. The routines have been aquired over 
decades, originated from NCL's Fortran77 codebase, the WRF model itself, 
or other tools like RIP (Read Interpolate Plot), and do not necessarily 
conform to a common programming mindset (e.g. 1D arrays, 2D arrays, etc).

The raw Fortran routines are compiled in to the :mod:`wrf._wrffortran` 
extension module, but are not particularly useful for applications in their 
raw form. These routines are imported in the *extension.py* module, where 
additional functionality is added to make the routines more user friendly.

The typical behavior for a fully exported Fortran routine in *extension.py* 
is:

1) Verify that the supplied arguments are valid in shape. Although f2py does 
   this as well, the errors thrown by f2py are confusing to users, so this 
   step helps provide better error messages.

2) Allocate an ouput array based on the output shape of the algorithm, 
   number of "leftmost"[1]_ dimensions, and size of the data.
   
3) Iterate over the leftmost [1]_ dimensions and compute output for argument 
   data slices that are of the same dimensionality as the compiled algorithm. 
   
4) Cast the argument arrays (or array slices) in to the dtype used in the 
   compiled routine. For WRF data, the conversion is usually from a 4-byte 
   float to an 8-byte double.
   
5) Extract the argument arrays out of xarray in to numpy arrays 
   (if applicable) and transpose them in to Fortran ordering. Note that this 
   does not actually do any copying of the data, it simply reorders the shape 
   tuple for the data and sets the Fortran ordering flag. This allows data
   pointers from the output array slices to be passed directly to the 
   Fortran routine, which eliminates the need to copy the result to the output 
   array.
   
The steps described above are handled in :mod:`wrapt` decorators that can be 
found in *decorators.py*. For some routines that produce multiple outputs or 
have atypical behavior, the special case decorators are located in 
*specialdec.py*. 

.. [1] If the Fortran algorithm is written for a 2-dimensional array, 
       and a users passes in a 5-dimensional array, there are 3 "leftmost" 
       dimensions.


Example
----------------------------

The above overviews are better explained by an example. Although there are a 
few exceptions (e.g. ll_to_xy), many of the routines in WRF-Python behave this 
way. 

For this example, let's make a routine that adds a variable's base state  
to its perturbation. This is the kind of thing that you'd normally use numpy 
for (e.g. Ptot = PB + P), but you could do this if you wanted concurrency 
for this operation via OpenMP rather than using dask (in a future release of 
WRF-Python, both OpenMP and dask will be available). 


Fortran Code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Below is the Fortran 90 code, which will be written to a file called 
example.f90. 

.. code:: fortran

   SUBROUTINE pert_add(base, pert, total, nx, ny)

   !f2py threadsafe
   !f2py intent(in,out) :: result
   
   REAL(KIND=8), INTENT(IN), DIMENSION(nx, ny) :: base, pert
   REAL(KIND=8), INTENT(OUT), DIMENSION(nx, ny) :: total
   INTEGER, INTENT(IN) :: nx, ny

   INTEGER :: i

   !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(runtime)
   DO j=1, ny
       DO i=1, nx
           total(i, j) = base(i, j) + pert(i, j)
       END DO
   END DO
   !$OMP END PARALLEL DO


   END SUBROUTINE pert_add

This code adds the 2D base and perturbation variables and stores the result in 
a 2D output array. (For this example, we're using a 2D array to help 
illustrate leftmost indexing below, but it could have been written using 
a 1D or 3D array). 

At the top, there are these two f2py directives:

.. code::

   !f2py threadsafe
   !f2py intent(in,out) :: total
   
The *threadsafe* directive tells f2py to release Python's Global Interpreter 
Lock (GIL) before calling the Fortran routine. The Fortran code no longer 
uses Python variables, so you should relese the GIL before running the 
computation. This way, Python threads will contine to run, which may be 
important if you are using this in a webserver or in some other 
threaded environment like dask's threaded scheduler. 

The *intent(in,out)* f2py directive is used because we will 
be supplying a slice of the output array directly to this routine and don't 
want to have to copy the result from Fortran back in to the result array. By 
specifying intent(in,out), we're telling f2py to use the pointer to our 
output array directly.

Finally, for the OpenMP directive, the scheduler is set to use runtime 
scheduling via *SCHEDULE(runtime)*. By using runtime scheduling, users 
can set the scheduling type within Python, but for most users the default is 
sufficient.


Building the Fortran Code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To build the Fortran code, the *example.f90* source code should be placed in 
the *fortran* directory of the source tree. 

Next, we need to update the numpy.distutils.core.Extension section of 
*setup.py* in the root directory of the source tree.

.. code:: python

   ext1 = numpy.distutils.core.Extension(
   name="wrf._wrffortran",
   sources=["fortran/wrf_constants.f90",
            "fortran/wrf_testfunc.f90",
            "fortran/wrf_user.f90",
            "fortran/rip_cape.f90",
            "fortran/wrf_cloud_fracf.f90",
            "fortran/wrf_fctt.f90",
            "fortran/wrf_user_dbz.f90",
            "fortran/wrf_relhl.f90",
            "fortran/calc_uh.f90",
            "fortran/wrf_user_latlon_routines.f90",
            "fortran/wrf_pvo.f90",
            "fortran/eqthecalc.f90",
            "fortran/wrf_rip_phys_routines.f90",
            "fortran/wrf_pw.f90",
            "fortran/wrf_vinterp.f90",
            "fortran/wrf_wind.f90",
            "fortran/omp.f90",
            "fortran/example.f90 # New file added here
            ]
    )
    
The easiest way to build your code is to use one of the build scripts located 
in the *build_scripts* directory of the source tree. These scripts contain 
variants for compiling with or without OpenMP support. Unless you are 
debugging a problem, building with OpenMP is recommended. 

For this example, we're going to assume you already followed how to 
:ref:`dev_setup`. Below are the build instructions for compiling with 
OpenMP enabled on GCC (Linux or Mac):

.. code::

   pip uninstall wrf-python
   cd build_scripts
   sh ./gnu_omp.sh
   
The above command will build and install the new routine, along with the 
other Fortran routines. If you recieve errors, then your code failed to 
build sucessfully. Otherwise, your new routine can be called as 
:meth:`wrf._wrffortran.pert_add`. 


Creating a Thin Python Wrapper
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The new Fortran pert_add routine will work well for a 2D slice of data. 
However, if you want to extend the functionality
to work with any dimensional array, you'll need to add a thin wrapper 
with some extra functionality that make use of :mod:`wrapt` decorators.

First, let's start by creating a very thin wrapper in Python in *extension.py*.

.. code:: python
    
   from wrf._wrffortran import pert_add
   
   .
   .
   .
   
   def _pert_add(base, pert, outview=None):
       """Wrapper for pert_add.

       Located in example.f90.

       """
       if outview is None:
           outview = np.empty(base.shape[0:2], base.dtype, order="F")

       result = pert_add(base,
                         pert,
                         outview)

       return result

Despite being only a few lines of code, there is quite a bit going on in the 
wrapper. The first thing to note is the arguments to the wrapper function. The
only arguments that we need for the wrapper are the inputs to the function 
and an "outview" keyword argument. At this point in the call chain, the 
arguments are assumed to be Fortran-ordered, in that the Fortran ordering flag 
is set and the shape is transposed from a usual C-ordered numpy array 
(the data itself remains in the same order that it was created). By passing 
numpy arrays with the Fortran order flag set, f2py will pass the pointer 
directly through to the Fortran routine.

The *outview* keyword argument is used during leftmost dimension indexing to 
send slices of the output array to the Fortran routine to be filled. If there 
are no leftmost dimensions (e.g. this routine is called with 2D data), then the 
outview argument will be None and an outview variable will be created with the 
same number of dimensions as the *base* argument. It should be created with 
Fortran ordering so that the pointer is directly passed to the Fortran routine.

When the actual :meth:`wrf._wrffortran.pert_add` Fortran routine is called, 
the nx and ny arguments are ommitted because f2py will supply this for us 
based on the shape of the numpy arrays we are supplying as input arguments. 
F2py also likes to return an array as a result, so even though we supplied 
outview as an array to be filled by the Fortran routine, we will still get a 
result from the function call that is pointing to the same thing as outview. 
(We could have chosen to ignore the result and return outview instead).


Extract and Transpose
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   
The arrays that are being passed to the _pert_add thin wrapper need to be 
numpy arrays in Fortran ordering, but they won't come this way from 
users. They will come in as either :class:`numpy.ndarray` 
or :class:`xarray.DataArray` and will be C-ordered. So, we need to make 
sure that a Fortran-ordered :class:`numpy.ndarray` is what is passed to 
the thin wrapper.

Since this type of operation is repeated for many diagnostic functions, a 
decorator has been written in *decorators.py* for this purpose. Let's decorate 
our thin wrapper with this function.


.. code:: python
    
   @extract_and_transpose()
   def _pert_add(base, pert, outview=None):
       """Wrapper for pert_add.

       Located in example.f90.

       """
       if outview is None:
           outview = np.empty(base.shape[0:2], base.dtype, order="F")

       result = pert_add(base,
                         pert,
                         outview)

       return result


The :meth:`extract_and_transpose` decorator converts any argument to _pert_add
that are of type :class:`xarray.DataArray` to :class:`numpy.ndarray`, and then 
gets the :attr:`numpy.ndarray.T` attribute, and passes this on to the 
_pert_add wrapper.

Following the computation, we want the result to be returned back as the 
same C-ordered array types that went in as arguments, so this decorator takes 
the result of the computation and returns the :attr:`numpy.ndarray.T` from the 
Fortran-ordered result. This result gets passed back up the decorator chain.


Cast to Fortran Array Types
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Fortran routine expects a specific data type for the arrays (usually 
REAL(KIND=8)). WRF files typically store their data as 4-byte floating point 
numbers to save space. The arrays being passed to the 
:meth:`wrf.decorators.extract_and_transpose` decorator need to be converted 
to the type used in the Fortran routine (e.g. double), then converted back to 
the original type (e.g. float) after the computation is finished. This is 
handled by the :meth:`wrf.decorators.cast_type` decorator function in 
*decorators.py*.

.. code:: python
   
   @cast_type(arg_idxs=(0, 1))
   @extract_and_transpose()
   def _pert_add(base, pert, outview=None):
       """Wrapper for pert_add.

       Located in example.f90.

       """
       if outview is None:
           outview = np.empty(base.shape[0:2], base.dtype, order="F")

       result = pert_add(base,
                         pert,
                         outview)

       return result
       
The :meth:`wrf.decorators.cast_type` decorator function takes an 
*arg_idxs* argument to specify which positional arguments need to be cast to 
the Fortran algorithm type, in this case arguments 0 and 1 (base and pert). 

Following the computation, the result will be cast back to the original type 
for the input arguments (usually float), and passed back up the decorator 
chain.


Leftmost Dimension Indexing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The WRF-Python algorithms written in Fortran are usually written for fixed 
size arrays of 1, 2, or 3 dimensions. If your input arrays have more than 
the number of dimensions specified for the Fortran algorithm, then we need to 
do the following:

1. Determine how many leftmost dimensions are used.

2. Create an output array that has a shape that contains the leftmost 
   dimensions concatenated with the shape of the result from the Fortran 
   algorithm.
   
3. Iterate over the leftmost dimensions and send slices of the input arrays 
   to the Fortran algorithm.
   
4. Along with the input arrays above, send a slice of the output array to be 
   filled by the Fortran algorithm.
   
5. Return the fully calculated output array.
   
The :meth:`wrf.decorators.left_iteration` is general purpose decorator 
contained in *decorators.py* to handle most leftmost index iteration cases. 
(Note: Some products, like cape_2d, return multiple products in the output 
and don't fall in to this generic category, so those decorators can be found 
in *specialdec.py*)

Let's look at how this is used below.

.. code:: python

   @left_iteration(2, 2, ref_var_idx=0)
   @cast_type(arg_idxs=(0, 1))
   @extract_and_transpose()
   def _pert_add(base, pert, outview=None):
       """Wrapper for pert_add.

       Located in example.f90.

       """
       if outview is None:
           outview = np.empty(base.shape[0:2], base.dtype, order="F")

       result = pert_add(base,
                         pert,
                         outview)

       return result
   

The :meth:`wrf.decorators.left_iteration` decorator handles many different 
use cases with its arguments, but this example is one of the more common cases. 
The 0th positional argument tells the decorator that the "reference" input 
variable should provide at least two dimensions. This should be set to 
the same number of dimensions as in the Fortran algorithm, which is two in this 
case. Dimensions to the left of these two dimensions are considered "leftmost" 
dimensions. 

The next positional argument (value of 2) tells the decorator that the 
newly created output variable should retain the shape of the reference 
variable's right two dimensions. This only applies when your output has less 
dimensions than the reference variable (e.g. sea level pressure uses 
geopotential height for the reference but produces 2D output). Since we are 
not reducing the output dimensions, it should be set to the same value as the 
previous argument. 

The final keyword argument of *ref_ver_idx* tells the decorator to use 
positional argument 0 (for the _pert_add function) as the reference 
variable. 

The result of this decorator will be the fully computed output array, which 
gets passed back up the chain.


Checking Argument Shapes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before any computations can be performed, the argument shapes are checked to 
verify their sizes. Although f2py will catch problems at the 
entry point to the Fortran routine, the error thrown is confusing to 
users. 

The :meth:`wrf.decorators.check_args` decorator is used to verify that the 
arguments are the correct size before proceeding. 

Here is how it is used:


.. code:: python

   @check_args(0, 2, (2, 2))
   @left_iteration(2, 2, ref_var_idx=0)
   @cast_type(arg_idxs=(0, 1))
   @extract_and_transpose()
   def _pert_add(base, pert, outview=None):
       """Wrapper for pert_add.

       Located in example.f90.

       """
       if outview is None:
           outview = np.empty(base.shape[0:2], base.dtype, order="F")

       result = pert_add(base,
                         pert,
                         outview)

       return result

The 0th positional argument (value of 0), tells 
:meth:`wrf.decorators.check_args` that the 0th positional argument of 
_pert_add is the reference variable. 

The next postional argument (value of 2) tells 
:meth:`wrf.decorators.check_args` that it should expect at least 2 dimensions 
for the reference variable. This should be set to the number of dimensions 
used in the Fortran algorithm, which is two in this case.

The final positional argument is a tuple with the number of dimensions that 
are expected for each array argument. Again, this should be set to the same 
number of dimensions expected in the Fortran routine for each positional 
argument. If an argument to your wrapped function is not an array type, you 
can use None in the tuple to ignore it, but that is not applicable for this 
example.


Putting It All Together
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The previous sections showed how the decorator chain was built up from the 
_pert_add function. However, when you actually make a call to _pert_add, the 
decorators are called from top to bottom. This means check_args is called 
first, then left_iteration, then cast_type, then extract_and_transpose, 
and finally _pert_add. After _pert_add is finished, the result is passed 
back up the chain and back to the user.

Now that we have a fully wrapped compiled routine, how might we use this?

Let's make a new :meth:`wrf.getvar` diagnostic called 'total_pressure'. A  
similar diagnostic already exists in WRF-Python, but this is just for 
illustration of how to use our newly wrapped Fortran routine.


Make a 'getter' Function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, we need a 'getter' routine that extracts the required input variables 
from the WRF NetCDF file(s) to perform the computation. In this case, the 
variables are P and PB.

The current naming convention in WRF-Python is to prefix the 'getter' 
functions with a '\g\_', so let's call this file g_totalpres.py and make a 
function get_total_pressure inside of it. 

The contents of this file will be:

.. code:: python

   # g_totalpres.py
   
   from .extension import _pert_add
   from .util import extract_vars

   @copy_and_set_metadata(copy_varname="P", name="total_pressure",
                          description="total pressure",
                          units="Pa")
   def get_total_pressure(wrfin, timeidx=0, method="cat", squeeze=True, 
                          cache=None, meta=True, _key=None):  
       """Return total pressure.

        This functions extracts the necessary variables from the NetCDF file
        object in order to perform the calculation.
    
        Args:
    
            wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
                iterable): WRF-ARW NetCDF
                data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
                or an iterable sequence of the aforementioned types.
    
            timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The
                desired time index. This value can be a positive integer,
                negative integer, or
                :data:`wrf.ALL_TIMES` (an alias for None) to return
                all times in the file or sequence. The default is 0.
    
            method (:obj:`str`, optional): The aggregation method to use for
                sequences.  Must be either 'cat' or 'join'.
                'cat' combines the data along the Time dimension.
                'join' creates a new dimension for the file index.
                The default is 'cat'.
    
            squeeze (:obj:`bool`, optional): Set to False to prevent dimensions
                with a size of 1 from being automatically removed from the 
                shape of the output. Default is True.
    
            cache (:obj:`dict`, optional): A dictionary of (varname, ndarray)
                that can be used to supply pre-extracted NetCDF variables to 
                the computational routines.  It is primarily used for internal
                purposes, but can also be used to improve performance by
                eliminating the need to repeatedly extract the same variables
                used in multiple diagnostics calculations, particularly when 
                using large sequences of files.
                Default is None.
    
            meta (:obj:`bool`, optional): Set to False to disable metadata and
                return :class:`numpy.ndarray` instead of
                :class:`xarray.DataArray`.  Default is True.
    
            _key (:obj:`int`, optional): A caching key. This is used for 
                internal purposes only.  Default is None.
    
        Returns:
            :class:`xarray.DataArray` or :class:`numpy.ndarray`: Omega.
            If xarray is
            enabled and the *meta* parameter is True, then the result will be a
            :class:`xarray.DataArray` object.  Otherwise, the result will be a
            :class:`numpy.ndarray` object with no metadata.
    
       """ 
       
       # Get the base and perturbation pressures
       varnames = ("PB", "P")
       ncvars = extract_vars(wrfin, timeidx, varnames, method, squeeze, cache,
                             meta=False, _key=_key)

       pb = ncvars["PB"]
       p = ncvars["P"]

       total_pres = _pert_add(pb, p)

       return total_pres


This getter function extracts the PB and P (base and pertrubation pressure) 
variables and calls the _pert_add function and returns the result. The 
arguments *wrfin*, *timeidx*, *method*, *squeeze*, *cache*, *meta*, and 
*_key* are used for every getter function and you can read what they do in 
the docstring. 

The getter function is also decorated with a  
:meth:`wrf.decorators.copy_and_set_metadata` decorator. This is a general 
purpose decorator that is used for copying metadata from an input variable 
to the result. In this case, the variable to copy is 'P'. The *name* parameter 
specifies that the :attr:`xarray.DataArray.name` attribute for the variable 
(the name that will be written to a NetCDF variable). The *description* is a 
brief description for variable that will be placed in the 
:attr:`xarray.DataArray.attrs` dictionary along with the *units* parameter.


Make Your New Diagnostic Available in :meth:`wrf.getvar`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The final step is to make the new 'total_pressure' diagnostic available from 
:meth:`wrf.getvar`.  To do this, modifications need to be made to 
*routines.py*.

First, import your new getter routine at the top of routines.py.

.. code:: python
   
   from __future__ import (absolute_import, division, print_function)

   from .util import (get_iterable, is_standard_wrf_var, extract_vars, 
                      viewkeys, get_id)
   from .g_cape import (get_2dcape, get_3dcape, get_cape2d_only,
                        get_cin2d_only, get_lcl, get_lfc, get_3dcape_only,
                        get_3dcin_only)
   .
   .
   .
   from .g_cloudfrac import (get_cloudfrac, get_low_cloudfrac, 
                             get_mid_cloudfrac, get_high_cloudfrac)
   from .g_totalpres import get_total_pressure


Next, update _FUNC_MAP to map your diagnostic label ('total_pressure') 
to the getter routine (get_total_pres).

.. code:: python

   _FUNC_MAP = {"cape2d": get_2dcape,
                "cape3d": get_3dcape,
                .
                .
                .
                "high_cloudfrac": get_high_cloudfrac,
                "total_pressure": get_total_pressure
                }
                

Finally, update _VALID_KARGS to inform :meth:`wrf.getvar` of any additional 
keyword argument names that this routine might use. The :meth:`wrf.getvar` 
routine will check keyword arguments and throws an error when it gets any that 
are not declared in this map.
 
In this case, there aren't any addtional keyword arguments, so we'll just 
supply an empty list.

.. code:: python

   _VALID_KARGS = {"cape2d": ["missing"],
                   "cape3d": ["missing"],
                   "dbz": ["do_variant", "do_liqskin"],
                   "maxdbz": ["do_variant", "do_liqskin"],
                   .
                   .
                   .
                   "high_cloudfrac": ["vert_type", "low_thresh",
                                      "mid_thresh", "high_thresh"],
                   "total_pressure": []
                   }
                   
After this is complete, your new routine is now available for use from 
:meth:`wrf.getvar`.




       