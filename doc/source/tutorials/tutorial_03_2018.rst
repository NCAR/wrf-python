WRF-Python Tutorial 2018
=========================

NCAR will be providing a four hour tutorial for wrf-python on Wednesday, March
7, 2018 from 8 AM to 12 PM. The tutorial is free, but seating is limited to 
only 16 students, so registration is required. 

The tutorial will take place at NCAR's corporate training center in Boulder, 
Colorado. 

`Corporate Technical Training Center <https://www2.fin.ucar.edu/it/about-cttc>`_
3085 Center Green Drive, Building CG-2, Room #3024
Boulder, Colorado

Overview
--------------

WRF-Python is a collection of diagnostic and interpolation routines for use 
with output from the Weather Research and Forecasting (WRF-ARW) Model. The 
package provides over 30 diagnostic calculations, 
several interpolation routines, and utilities to help with plotting 
via cartopy, basemap, or PyNGL. The functionality is similar to what is 
provided by the NCL WRF package. 

.. note:: 

   WRF-Python is NOT a tool for running the WRF-ARW model using Python.

This tutorial provides an introduction to wrf-python. The tutorial is beginner 
friendly for new users of wrf-python, but this is NOT an introduction to the 
Python programming language (see :ref:`prereq`). Due to limited seating, if you 
do not have any previous experience with Python, please do not register 
for this tutorial.

.. note::

   For online training that provides an introduction to the Python 
   programming language itself, please see the 
   `Unidata Python Training Page <https://unidata.github.io/online-python-training/>`_.

Computers will be provided, but feel free to use your own laptop if you prefer. 
We will be covering how to install wrf-python via conda as part of the 
tutorial.

Students are encouraged to bring their own data sets, but data will be provided
if this is not an option. Students will be provided a jupyter notebook workbook
which can be modified to accommodate their data. 

Topics include:

- How to install wrf-python via conda
- A brief introduction to jupyter notebook
- Overview of WRF data files
- WRF-Python basics
- Plotting with cartopy
- Overview of OpenMP features and other performance tips
- Open lab for students 


Registration
---------------

Please register prior to February 21, 2018. The registration form is here:

`Registration Form <https://goo.gl/forms/is5VExf3w4bFGXUb2>`_

Registration consists of a brief survey, which will help give the instructor
a brief overview of your background and will help tailor the tutorial to 
your expectations.

.. _prereq:

Prerequisites
---------------

This tutorial assumes that you have basic knowledge of how to type commands 
in to a terminal window using your preferred operating system.  You 
should know some basic directory commands like *cd*, *mkdir*, *cp*, *mv*.

This tutorial assumes that you have prior experience programming in Python.
Below is a list of some Python concepts that you will see in the examples, 
but don't worry if you aren't familiar with everything.  

- Opening a Python interpreter and entering commands.
- Importing packages via the import statement.
- Familiarity with some of the basic Python types: str, list, tuple, dict, bool, float, int, None.
- Creating a list, tuple, or dict with "[ ]", "( )", "{ }" syntax (e.g. my_list = [1,2,3,4,5]).
- Accessing dict/list/tuple items with the "x[ ]" syntax (e.g. my_list_item = my_list[0]).
- Slicing str/list/tuple with the ":" syntax (e.g. my_slice = my_list[1:3]).
- Using object methods and attributes with the "x.y" syntax (e.g. my_list.append(6)).
- Calling functions (e.g. result = some_function(x, y))
- Familiarity with numpy would be helpful, as only a very brief introduction
  is provided.
- Familiarity with matplotlib would be helpful, as only a very brief 
  introduction is provided.
  


