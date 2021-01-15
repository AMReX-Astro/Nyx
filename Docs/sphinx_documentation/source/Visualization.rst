
 .. role:: cpp(code)
    :language: c++

.. _Visualization:

Visualization
=============

Nyx generates plotfile in the native AMReX format as well as HDF5.

There are several visualization tools that can be used for AMReX plotfiles, 
specifically VisIt, ParaView, and yt.

See the `AMReX documentation <https://amrex-codes.github.io/amrex/docs_html/Visualization_Chapter.html>`_ for 
guidance on how to use each of these tools.

The following inputs control when plotfiles are written

+-------------------------+-----------------------------------------------------------------------+-------------+-----------+
|                         | Description                                                           |   Type      | Default   |
+=========================+=======================================================================+=============+===========+
| amr.plot_int            | Frequency of plotfile output;                                         |    Int      | -1        |
|                         | if -1 then no plotfiles will be written                               |             |           |
+-------------------------+-----------------------------------------------------------------------+-------------+-----------+
| amr.plotfile_on_restart | Should we write a plotfile when we restart (only used if plot_int>0)  |   Bool      | False     |
+-------------------------+-----------------------------------------------------------------------+-------------+-----------+
| amr.plot_file           | Prefix to use for plotfile output                                     |  String     | plt       |
+-------------------------+-----------------------------------------------------------------------+-------------+-----------+

and the following control which variables appear in the plotfile

+----------------------------+---------------------------------------------------+------------+-----------+
|                            | Description                                       |   Type     | Default   |
+============================+===================================================+============+===========+
| amr.plot_vars              | Name of state variables to be in plotfiles        |   Strings  | All       |
+----------------------------+---------------------------------------------------+------------+-----------+
| amr.plot_dervive_plot_vars | Name of derived variables to be in plotfiles      |   Strings  | All       |
+----------------------------+---------------------------------------------------+------------+-----------+
