Running the Code
================

The Nyx executable reads run-time information from an “inputs” file which you designate on the command line.
Values can be specified either in the inputs file or on the command line.
If a value is specified on the command line, that value will override a value specified in the inputs file.

See the Inputs section for more detail about what parameters can be specified at run-time.

#. From within the build directory type:

   ::

      <executable-name> inputs.32

   where ``<executable-name>`` is  ``Nyx3d.Linux.gnu.ex1`` if you built Nyx with GNU Make, or
   ``nyx_<name-of-build-directory>`` if you built Nyx with CMake.

#. You will notice that running the code generates directories that look like
   plt00000, plt00020, etc.,
   and chk00000, chk00020, etc. These are “plotfiles” and
   “checkpoint” files. The plotfiles are used for visualization,
   and the checkpoint files for restarting the code.

See the Visualization chapter for how to visualize these plotfiles.

.. include:: NyxSundials.rst
