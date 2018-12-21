.. _diagnostic-table:

Table of Available Diagnostics
=================================

.. include:: _templates/product_table.txt


.. _subdiagnostic-table:

Table of Subproduct Diagnostics
----------------------------------

Some diagnostics (e.g. cape_2d) include multiple products in its 
output. These products have been broken out in to individual diagnostics 
to help those utilities that are unable to work with multiple outputs. 
These individual diagnostics can be requested like any other diagnostic 
using :meth:`wrf.getvar`. These are summarized in the table below.


.. include:: _templates/subproducts.txt
