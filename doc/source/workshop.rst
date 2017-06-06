WRF Workshop 2017
=====================

Welcome wrf-python tutorial attendees!

The instructions below should be completed prior to arriving at the tutorial.  
There will not be enough time to do this during the tutorial.

Prerequisites
---------------

This tutorial assumes that you have basic knowledge of how to type commands 
in to a command terminal using your preferred operating system.  You 
should know some basic directory commands like *cd*, *mkdir*, *cp*, *mv*.

Regarding Python, to understand the examples in this tutorial, you
should have some experience with Python basics.  This includes:

- Opening a Python interpreter and entering commands.
- Importing packages via the import statement.
- Using the basic Python types: str, list, tuple, dict, boolean, numeric, None.
- Creating a list, tuple, or dict with "[ ]", "( )", "{ }" syntax.
- Accessing dict/list/tuple items with the "x[ ]" syntax.
- Slicing str/list/tuple with the ":" syntax.
- Using object methods and attributes with the "x.y" syntax.
- Familiarity with numpy would be helpful, as only a very brief introduction
  is provided.
- Familiarity with matplotlib would be helpful, as only a very brief 
  introduction is provided.
  
If you are completely new to Python, that shouldn't be a problem, since 
most of the examples consist of basic container types and function calls.  It 
would be helpful to look at some introductory material before arriving at the 
tutorial.  If you've programmed before, picking up Python is fairly easy.  

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
going to use Python 2.7, but it should also work with Python 3.5+.  However, 
due to limitations with open source compilers on conda-forge, only Python 2.7 
is available for Windows.

Please use the appropriate link below to download Miniconda for your operating 
system. 

.. note:: 

   64-bit OS recommended  

`Win64 <https://repo.continuum.io/miniconda/Miniconda2-latest-Windows-x86_64.exe>`_

`Mac <https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh>`_

`Linux <https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh>`_

For more information, see: https://conda.io/miniconda.html

.. note::

    **What is Miniconda?**

    If you have used the Anaconda distribution for Python before, then you will be 
    familiar with Miniconda.  The Anaconda Python distribution includes numerous 
    scientific packages out of box, which can be difficult for users to build and 
    install. More importantly, Anaconda includes the conda package manager. 
    
    The conda package manager is a utility (similar to yum or apt-get) that installs 
    packages from a repository of pre-compiled Python packages.  These repositories 
    are called channels.  Conda makes it easy for Python users to install and 
    uninstall packages, and also can be used to create isolated Python environments 
    (more on that later).
    
    Miniconda is a bare bones implementation of Anaconda and only includes the 
    conda package manager.  Since we are going to use the conda-forge channel to 
    install our scientific packages, Miniconda avoids any complications between 
    packages provided by Anaconda and conda-forge. 


Step 3: Install Miniconda
----------------------------

Windows:

    1. Browse to the directory where you downloaded Miniconda2-latest-Windows-x86_64.exe.
    2. Double click on Miniconda2-latest-Windows-x86_64.exe.  
    3. Follow the instructions.
    4. Restart your command terminal.
    
Mac and Linux:

    For Mac and Linux, the installer is a bash script.  
    
    1. Using a terminal, you need to execute the shell script that you downloaded by
       doing::
    
            bash /path/to/Miniconda2-latest-MacOSX-x86_64.sh [Mac]
            
            bash /path/to/Miniconda2-latest-Linux-x86_64.sh [Linux]
    
    2. Follow the instructions.  
    
    3. At the end of the installation, it will ask if you want to add the 
       miniconda2 path to your bash environment.  If you are unsure what to do,
       you should say "yes".  If you say "no", we're going to assume you know
       what you are doing.
       
       If you said "yes", then once you restart your shell, the miniconda2 Python 
       will be found instead of the system Python when you type the "python" 
       command.  If you want to undo this later, then you can edit 
       either ~/.bash_profile or ~/.bashrc (depending on OS used) and 
       comment out the line that looks similar to::
    
            # added by Miniconda2 4.1.11 installer
            export PATH="/path/to/miniconda2/bin:$PATH"
            
    4. Restart your command terminal.
    
.. note:: 

   If you have already installed another Python distribution, like Enthought 
   Canopy, you will need to comment out any PATH entries for that distribution
   in your .bashrc or .bash_profile.  Otherwise, your shell environment may 
   pick to wrong Python installation.  You can use the command "which python" 
   to verify the Python installation being used.
   
    
    
Step 4: Set Up the Conda Environment
--------------------------------------

If you are new to the conda package manager, one of the nice features of conda 
is that you can create isolated Python environments that prevent package 
incompatibilities. This is similar to the virtualenv package that some 
Python users may be familiar with.  However, conda is not compatible with 
virtualenv, so only use conda environments when working with conda.

The name of our conda environment for this tutorial is: **tutorial_2017**.

Follow the instructions below to create the tutorial_2017 environment.

   1. Open a command terminal if you haven't done so.
   
   2. [Linux and Mac Users Only] The conda package manager only works with bash, 
      so if bash is not your current shell, type::
      
          bash
      
   3. Add the conda-forge channel to your conda package manager. 
   
      Type or copy this command in to your command terminal::
   
          conda config --add channels conda-forge
        
      .. note:: 
         
         Conda-forge is a community driven collection of packages that are 
         continually tested to ensure compatibility.  We highly recommend using
         conda-forge when working with conda.  See https://conda-forge.github.io/
         for more details on this excellent project.
        
   4. Create the conda environment for the tutorial.
   
      Type or copy this command in to your command terminal::
      
          conda create -n tutorial_2017 python=2.7 matplotlib=1.5.3 cartopy netcdf4 jupyter git ffmpeg wrf-python
          
      Type "y" when prompted.  It will take several minutes to install everything.
          
      This command creates an isolated Python environment named *tutorial_2017*, and installs 
      the python interpreter, matplotlib, cartopy, netcdf4, jupyter, git, ffmpeg, and wrf-python 
      packages.  
      
      .. note:: 
         
         In this tutorial, we need to use matplotlib v1.5.3 due to some issues with cartopy, which 
         should be fixed in a later version of cartopy.  Be sure to supply the version number as 
         indicated in the command above.
      
   5. Activate the conda environment.
   
      To activate the tutorial_2017 Python environment, type the following 
      in to the command terminal:
      
      For Linux and Mac (using bash)::
          
          source activate tutorial_2017
          
      For Windows::
      
          activate tutorial_2017
          
      You should see (tutorial_2017) on your command prompt.
      
      To deactivate your conda environment, type the following in to the 
      command terminal:
      
      For Linux and Mac::
      
          source deactivate
          
      For Windows::
      
          deactivate tutorial_2017
      

Step 5: Download the Student Workbook
---------------------------------------

The student workbook for the tutorial is available on GitHub.  The tutorial_2017 
conda environment includes the git application needed to download the repository.

These instructions download the tutorial in to your home directory.  If you want 
to place the tutorial in to another directory, we're going to assume you know 
how to do this yourself.

To download the student workbook, follow these instructions:

    1. Activate the tutorial_2017 conda environment following the instructions 
       in the previous step (*source activate tutorial_2017* or 
       *activate tutorial_2017*).
    
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
       
           source activate tutorial_2017
           
           cd ~/wrf_python_tutorial
           
           git pull
           
       For Windows::
       
           activate tutorial_2017
           
           cd %HOMEPATH%\wrf_python_tutorial
           
           git pull
    

Step 6: Obtain WRF Output Files
----------------------------------

For this tutorial, we strongly recommend that you use your own WRF output files.  
The tutorial includes an easy way to point to your own data files.  The WRF 
output files should all be from the same WRF run and use the same domain.  
If your files are located on another system (e.g. yellowstone), then copy 2 or 
3 of these files to your local computer prior to the tutorial.

If you do not have any of your own WRF output files, then you can download the 
instructor data files from a link that should have been provided to you in an 
email prior to the tutorial.


Getting Help
----------------

If you experience problems during this installation, please send a question 
to the :ref:`google-group` support mailing list.  


We look forward to seeing you at the tutorial!
