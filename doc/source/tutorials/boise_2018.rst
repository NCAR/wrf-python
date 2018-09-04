WRF-Python and VAPOR Workshop 2018 (Boise State University)
=============================================================

The Department of Geosciences at Boise State University is partnering with 
staff from the National Center for Atmospheric Research (NCAR) to host a free, 
2-day workshop in the Environmental Research Building (ERB) lab 2104 at 
Boise State University on September 26-27, 2018. The tutorial will be centered 
on the WRF-Python and VAPOR tools for analyzing and visualizing data from the 
Weather Research and Forecasting (WRF) regional weather and climate model. 

Users must be registered to attend this tutorial (see :ref:`registration`).

Location
---------------------

September 26-27, 2018 9:00 AM - 4:00 PM

Boise State University, Environmental Research Building (ERB) lab #2104.


WRF-Python Overview
---------------------

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
Python programming language (see :ref:`prereq_boise`). Due to limited seating, 
if you do not have any previous experience with Python, please do not register 
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

.. _registration:

Registration
---------------

Please register prior to September 19, 2018. The registration form is here:

`Registration Form <https://goo.gl/forms/ASb8bP7Bz2Boxye23>`_

Registration consists of a brief survey, which will help give the instructor
a brief overview of your background and will help tailor the tutorial to 
your expectations.

.. _prereq_boise:

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
  

-------------------------------------------------

Instructions for Computer Lab Installation
-------------------------------------------------

Step 1: Download Miniconda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For this tutorial, you will need to download and install Miniconda.  We are 
going to use Python 3.6+. 

Please use the appropriate link below to download Miniconda for your operating 
system. 

.. note:: 

   64-bit OS recommended  

`Win64 <https://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86_64.exe>`_

`Mac <https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh>`_

`Linux <https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh>`_

For more information, see: https://conda.io/miniconda.html


Step 2: Install Miniconda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Windows:

    1. Browse to the directory where you downloaded Miniconda3-latest-Windows-x86_64.exe.
    
    2. Double click on Miniconda3-latest-Windows-x86_64.exe. 
     
    3. Follow the instructions.
    
    4. For Windows 10, use the Anaconda command prompt found under the Anaconda2
       menu (Start Menu -> Anaconda2 -> Anaconda Prompt). Otherwise, open a 
       regular command prompt.
    
Mac and Linux:

    For Mac and Linux, the installer is a bash script. 
    
    1. Using a terminal, you need to execute the bash shell script that you downloaded by
       doing::
    
            bash /path/to/Miniconda3-latest-MacOSX-x86_64.sh [Mac]
            
            bash /path/to/Miniconda3-latest-Linux-x86_64.sh [Linux]
    
    2. Follow the instructions.  
    
    3. At the end of the installation, it will ask if you want to add the 
       miniconda3 path to your bash environment.  If you are unsure what to do,
       you should say "yes".  If you say "no", we're going to assume you know
       what you are doing.
       
       If you said "yes", then once you restart your shell, the miniconda3 Python 
       will be found instead of the system Python when you type the "python" 
       command.  If you want to undo this later, then you can edit 
       either ~/.bash_profile or ~/.bashrc (depending on OS used) and 
       comment out the line that looks similar to::
    
            # added by Miniconda3 x.x.x installer
            export PATH="/path/to/miniconda3/bin:$PATH"
            
    4. Restart your command terminal.
    
    5. [Linux and Mac Users Only] Miniconda only works with bash.  If bash is 
       not your default shell, then you need to activate the bash shell by typing 
       the following in to your command terminal::
       
           bash
           
    6. Verify that your system is using the correct Python interpreter by typing
       the following in to your command terminal::
       
           which python
           
       You should see the path to your miniconda installation.  If not, see the 
       note below. 
       
       .. note::

           If you have already installed another Python distribution, like Enthought 
           Canopy, you will need to comment out any PATH entries for that distribution
           in your .bashrc or .bash_profile.  Otherwise, your shell environment may 
           pick to wrong Python installation.
           
           If bash is not your default shell type, and the PATH variable has been 
           set in .bash_profile by the miniconda installer, try executing 
           "bash -l" instead of the "bash" command in step 5.  
           
   
Step 3: Set Up the Conda Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are new to the conda package manager, one of the nice features of conda 
is that you can create isolated Python environments that prevent package 
incompatibilities. This is similar to the *virtualenv* package that some 
Python users may be familiar with.  However, conda is not compatible with 
virtualenv, so only use conda environments when working with conda.

The name of our conda environment for this tutorial is: **tutorial_backup**.

Follow the instructions below to create the tutorial_backup environment.

   1. Open a command terminal if you haven't done so.
   
   2. [Linux and Mac Users Only] The conda package manager only works with bash, 
      so if bash is not your current shell, type::
      
          bash
      
   3. Add the conda-forge channel to your conda package manager. 
   
      Type or copy the command below in to your command terminal. You should 
      run this command even if you have already done it in the past.  
      This will ensure that conda-forge is set as the highest priority channel.
      
      :: 
   
          conda config --add channels conda-forge
          
      .. note:: 
         
         Conda-forge is a community driven collection of packages that are 
         continually tested to ensure compatibility.  We highly recommend using
         conda-forge when working with conda.  See https://conda-forge.github.io/
         for more details on this excellent project.
        
   4. Create the backup conda environment for the tutorial.
   
      Students will create a conda environment during the tutorial, but if 
      they run in to problems, we're going to create a backup environment.
   
      Type or copy this command in to your command terminal::
      
          conda create -n tutorial_backup python=3.6 matplotlib cartopy netcdf4 jupyter git ffmpeg wrf-python
          
      Type "y" when prompted.  It will take several minutes to install everything.
          
      This command creates an isolated Python environment named *tutorial_backup*, and installs 
      the python interpreter, matplotlib, cartopy, netcdf4, jupyter, git, ffmpeg, and wrf-python 
      packages.  
         
     .. note::
     
         When the installation completes, your command terminal might post a message similar to:
         
         .. code-block:: none
         
             If this is your first install of dbus, automatically load on login with:
             
             mkdir -p ~/Library/LaunchAgents
             cp /path/to/miniconda3/envs/tutorial_test/org.freedesktop.dbus-session.plist ~/Library/LaunchAgents/
             launchctl load -w ~/Library/LaunchAgents/org.freedesktop.dbus-session.plist
             
         This is indicating that the dbus package can be set up to automatically load on login.  You 
         can either ignore this message or type in the commands as indicated on your command terminal.  
         The tutorial should work fine in either case.
      
   5. Activate the conda environment.
   
      To activate the tutorial_backup Python environment, type the following 
      in to the command terminal:
      
      For Linux and Mac (using bash)::
          
          source activate tutorial_backup
          
      For Windows::
      
          activate tutorial_backup
          
      You should see (tutorial_backup) on your command prompt.
      
      To deactivate your conda environment, type the following in to the 
      command terminal:
      
      For Linux and Mac::
      
          source deactivate
          
      For Windows::
      
          deactivate tutorial_backup
      

Step 4: Download the Student Workbook
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The student workbook for the tutorial is available on GitHub.  The tutorial_backup
conda environment includes the git application needed to download the repository.

These instructions download the tutorial in to your home directory.  If you want 
to place the tutorial in to another directory, we're going to assume you know 
how to do this yourself.

To download the student workbook, follow these instructions:

    1. Activate the tutorial_backup conda environment following the instructions 
       in the previous step (*source activate tutorial_backup* or 
       *activate tutorial_backup*).
    
    2. Change your working directory to the home directory by typing the 
       following command in to the command terminal:
    
       For Linux and Mac:: 
       
           cd ~
           
       For Windows:: 
       
           cd %HOMEPATH%
           
    3. Download the git repository for the tutorial by typing the following 
       in to the command terminal::
       
           git clone https://github.com/NCAR/wrf_python_tutorial.git
           
    4. There may be additional changes to the tutorial after you have downloaded 
       it. To pull down the latest changes, type the following in to the 
       command terminal:
       
       For Linux and Mac::
       
           source activate tutorial_backup
           
           cd ~/wrf_python_tutorial/boise_workshop_2018
           
           git pull
           
       For Windows::
       
           activate tutorial_2018
           
           cd %HOMEPATH%\wrf_python_tutorial\boise_workshop_2018
           
           git pull
       
       .. note::
       
           If you try the "git pull" command and it returns an error indicating
           that you have made changes to the workbook, this is probably because 
           you ran the workbook and it contains the cell output.  To fix this, 
           first do a checkout of the workbook, then do the pull.  
           
           .. code-block:: none
           
               git checkout -- .
               git pull
               

Step 5:  Verify Your Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Verifying that your environment is correct involves importing a few 
packages and checking for errors (you may see some warnings for matplotlib 
or xarray, but you can safely ignore these). 

    1. Activate the tutorial_backup conda environment if it isn't already active 
       (see instructions above).
       
    2. Open a python terminal by typing the following in to the command 
       terminal::
       
           python
       
    3. Now type the following in to the Python interpreter::
    
           >>> import netCDF4
           >>> import matplotlib
           >>> import xarray
           >>> import wrf
       
   4. You can exit the Python interpreter using **CTRL + D**
    

Step 6: Install WRF Output Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A link will be provided in an email prior to the tutorial for the WRF-ARW 
data files used for the examples. 


    1. The link in the email should take you to a location on an Amazon cloud 
       drive.
       
    2. If you hover your mouse over the wrf_tutorial_data.zip file, you'll see 
       an empty check box appear next to the file name.  Click this check 
       box.
       
    3. At the bottom of the screen, you'll see a Download button next to a 
       cloud icon.  Click this button to start the download.
       
    4. The download was most likely placed in to your ~/Downloads folder 
       [%HOMEPATH%\\Downloads for Windows]. Using your preferred method of choice 
       for unzipping files, unzip this file in to your home directory.  Your data 
       should now be in ~/wrf_tutorial_data 
       [%HOMEPATH%\\wrf_tutorial_data for Windows].
       
    5. Verify that you have three WRF output files in that directory. 
  


