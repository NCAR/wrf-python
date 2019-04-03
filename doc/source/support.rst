Support
==========  

Github
----------------

For crashes and other issues related to the software, 
you can submit an issue to the 
`Github repository <https://github.com/NCAR/wrf-python>`_.

This should be used strictly for crashes and bugs.  For general usage
questions, please use the :ref:`google-group`.

.. _submitting_files:

Submitting Files
-------------------

When issues are encountered, we often need to see the file that is problematic.
As long as the file is not so large that it fills up the FTP server, users 
can submit files via FTP using the following instructions:

.. code-block:: none

   ftp ftp.cgd.ucar.edu
   anonymous
   <use your email address for the password>
   cd incoming
   put <your file>
   put <your file>
   .
   .
   .
   quit
   
.. note::
   
   For security reasons, you cannot list the contents of this directory, 
   and neither can we.  You must tell us the *exact* names of the 
   files that you uploaded in order for us to retrieve them.
   

.. _google-group:

Google Group
---------------

The best way to receive support is to use the 
`[wrfpython-talk] 
<https://groups.google.com/a/ucar.edu/d/forum/wrfpython-talk>`_ group on 
Google. Users can ask usage questions, submit bug reports, make feature 
requests, etc.  All users of wrf-python are welcome to join.

Below, you can view the most recent discussions:

-------------------------

.. raw:: html

   <iframe id="forum_embed"
      src="javascript:void(0)"
      scrolling="no"
      frameborder="0"
      width="700"
      height="700">
    </iframe>
    <script type="text/javascript">
      document.getElementById('forum_embed').src =
         'https://groups.google.com/a/ucar.edu/forum/embed/?place=forum/wrfpython-talk'
         + '&showsearch=true&showpopout=true&showtabs=false'
         + '&parenturl=' + encodeURIComponent(window.location.href);
    </script>
    