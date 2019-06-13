WRF Users' Workshop 2019
=========================

Welcome WRF-Python Tutorial attendees!

If you wish to actively participate in the tutorial, please bring your own 
laptop. **Due to limited time constraints, the instructions below should be 
completed prior to arriving at the tutorial**.

I will be executing the same cells as the student workbook, so if you prefer 
to sit and watch, that's OK too. Following the tutorial, I will upload 
the instructor slides in to the same GitHub location as the student workbook if 
you want to try it out later.

Prerequisites
---------------

This tutorial assumes that you have basic knowledge of how to type commands 
in to a command terminal using your preferred operating system.  You 
should know some basic directory commands like *cd*, *mkdir*, *cp*, *mv*.

Regarding Python, to understand the examples in this tutorial, you
should have some experience with Python basics.  Below is a list of some 
Python concepts that you will see in the examples, but don't worry if you aren't 
familiar with everything.  

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
  
If you are completely new to Python, that shouldn't be a problem, since 
most of the examples consist of basic container types and function calls.  It 
would be helpful to look at some introductory material before arriving at the 
tutorial.  If you've programmed before, picking up Python isn't too difficult.  

Here are some links:

https://www.learnpython.org/

https://developers.google.com/edu/python/


Step 1: Open a Command Terminal
--------------------------------

To begin, you will first need to know how to open a command line terminal for 
your operating system.   

For Windows:

.. code-block:: none

    WINDOWS + r
    type cmd in the run window
    
For Mac:

.. code-block:: none

     Finder -> Applications -> Utilities -> Terminal
     
For Linux:

.. code-block:: none

    Try one of the following:
    
    CTRL + ALT + T
    CTRL + ALT + F2


Step 2: Download Miniconda
----------------------------

For this tutorial, you will need to download and install Miniconda.  We are 
going to use Python 3.7, but it will also work with Python 2.7. 

Please use the appropriate link below to download Miniconda for your operating 
system. 

.. note:: 

   64-bit OS recommended  

`Win64 <https://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86_64.exe>`_

`Mac <https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh>`_

`Linux <https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh>`_

For more information, see: https://conda.io/miniconda.html

.. note::

    **What is Miniconda?**

    If you have used the Anaconda distribution for Python before, then you will 
    be familiar with Miniconda.  The Anaconda Python distribution includes numerous 
    scientific packages out of the box, which can be difficult for users to build and 
    install. More importantly, Anaconda includes the conda package manager. 
    
    The conda package manager is a utility (similar to yum or apt-get) that installs 
    packages from a repository of pre-compiled Python packages.  These repositories 
    are called channels.  Conda makes it easy for Python users to install and 
    uninstall packages, and also can be used to create isolated Python environments 
    (more on that later).
    
    Miniconda is a bare bones implementation of Anaconda and only includes the 
    conda package manager. Since we are going to use the conda-forge channel to 
    install our scientific packages, Miniconda avoids any complications between 
    packages provided by Anaconda and conda-forge. 


Step 3: Install Miniconda
----------------------------

Windows:

    1. Browse to the directory where you downloaded Miniconda3-latest-Windows-x86_64.exe.
    
    2. Double click on Miniconda3-latest-Windows-x86_64.exe. 
     
    3. Follow the instructions.
    
    4. Restart your command terminal.
    
Mac and Linux:

    For Mac and Linux, the installer is a bash script. 
    
    1. Using a terminal, you need to execute the bash shell script that you downloaded by
       doing::
    
            bash /path/to/Miniconda3-latest-MacOSX-x86_64.sh [Mac]
            
            bash /path/to/Miniconda3-latest-Linux-x86_64.sh [Linux]
    
    2. Follow the instructions.  
    
    3. At the end of the installation, it will ask if you "wish the installer to
       initialize Miniconda3".  If you are unsure what to do, you should say "yes".
       If you say "no", we're going to assume you know what you are doing.

       If you said "yes", then once you restart your shell, the "base" conda environment
       will be activated by default and the miniconda3 Python will be found instead
       of the system Python when you type the "python" command.  If you want to undo
       this later, then you can run the following "conda" command::

            conda config --set auto_activate_base false

       Note that this will not affect your ability to run "conda" commands, it will just
       prevent the miniconda3 Python from overriding any existing Python environments you
       may already have on your machine.

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
           pick the wrong Python installation.
           
           If bash is not your default shell type, and the PATH variable has been 
           set in .bash_profile by the miniconda installer, try executing 
           "bash -l" instead of the "bash" command in step 5.  
           
   
Step 4: Set Up the Conda Environment
--------------------------------------

If you are new to the conda package manager, one of the nice features of conda 
is that you can create isolated Python environments that prevent package 
incompatibilities. This is similar to the *virtualenv* package that some 
Python users may be familiar with.  However, conda is not compatible with 
virtualenv, so only use conda environments when working with conda.

The name of our conda environment for this tutorial is: **tutorial_2019**.

Follow the instructions below to create the tutorial_2019 environment.

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
        
   4. Create the conda environment for the tutorial.
   
      Type or copy this command in to your command terminal::
      
          conda create -n tutorial_2019 python=3.7 matplotlib cartopy netcdf4 jupyter git ffmpeg wrf-python
          
      Type "y" when prompted.  It will take several minutes to install everything.
          
      This command creates an isolated Python environment named *tutorial_2019*, and installs 
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
   
      To activate the tutorial_2019 Python environment, type the following 
      in to the command terminal:
      
      For Linux and Mac (using bash)::
          
          conda activate tutorial_2019
          
      For Windows::
      
          activate tutorial_2019
          
      You should see (tutorial_2019) on your command prompt.
      
      To deactivate your conda environment, type the following in to the 
      command terminal:
      
      For Linux and Mac::
      
          conda deactivate
          
      For Windows::
      
          deactivate tutorial_2019
      

Step 5: Download the Student Workbook
---------------------------------------

The student workbook for the tutorial is available on GitHub.  The tutorial_2019 
conda environment includes the git application needed to download the repository.

These instructions download the tutorial in to your home directory.  If you want 
to place the tutorial in to another directory, we're going to assume you know 
how to do this yourself.

To download the student workbook, follow these instructions:

    1. Activate the tutorial_2019 conda environment following the instructions 
       in the previous step (*conda activate tutorial_2019* or 
       *activate tutorial_2019*).
    
    2. Change your working directory to the home directory by typing the 
       following command in to the command terminal:
    
       For Linux and Mac:: 
       
           cd ~
           
       For Windows:: 
       
           cd %HOMEPATH%
           
    3. Download the git repository for the tutorial by typing the following 
       in to the command terminal::
       
           git clone --recursive https://github.com/NCAR/wrf_python_tutorial.git
           
    4. There may be additional changes to the tutorial after you have downloaded 
       it. To pull down the latest changes, type the following in to the 
       command terminal:
       
       For Linux and Mac::
       
           conda activate tutorial_2019
           
           cd ~/wrf_python_tutorial/wrf_workshop_2019
           
           git pull
           
       For Windows::
       
           activate tutorial_2019
           
           cd %HOMEPATH%\wrf_python_tutorial\wrf_workshop_2019
           
           git pull
       
       .. note::
       
           If you try the "git pull" command and it returns an error indicating
           that you have made changes to the workbook, this is probably because 
           you ran the workbook and it contains the cell output.  To fix this, 
           first do a checkout of the workbook, then do the pull.  
           
           .. code-block:: none
           
               git checkout -- .
               git pull
               

Step 6:  Verify Your Environment
----------------------------------

Verifying that your environment is correct involves importing a few 
packages and checking for errors (you may see some warnings for matplotlib 
or xarray, but you can safely ignore these). 

    1. Activate the tutorial_2019 conda environment if it isn't already active 
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
    

Step 7: Obtain WRF Output Files
----------------------------------

The wrf_python_tutorial git repository linked to in Step 5 includes a directory
containing several WRF-ARW data files which will be used for examples during
the tutorial.

You also have the option of using your own data files for the tutorial by 
modifying the first Jupyter Notebook cell to point to your data set. 
However, there is no guarantee that every cell in your workbook will work 
without some modifications (e.g. cross section lines will be drawn outside of 
your domain).

    1. If you have recently cloned the wrf_python_tutorial git repository, then
       you should have a "wrf_tutorial_data" directory at the root level of the
       "wrf_python_tutorial" directory.

    2. If this directory does not exist, try running the following commands
       from within the "wrf_python_tutorial" directory to update your local
       copy of the git repository::

           git checkout -- .
           git pull
           git submodule init
           git submodule update
           cd wrf_tutorial_data
           git checkout -- .

    3. Verify that you have three WRF output files in the "wrf_tutorial_data"
       directory::

           $ ls wrf_tutorial_data
           wrfout_d01_2005-08-28_00_00_00
           wrfout_d01_2005-08-28_12_00_00
           wrfout_d01_2005-08-29_00_00_00


Getting Help
----------------

If you experience problems during this installation, please send a question 
to the :ref:`google-group` support mailing list.  


We look forward to seeing you at the tutorial!
