Running the Code
================

The Nyx executable reads run-time information from an “inputs” file which you designate on the command line.
Values can be specified either in the inputs file or on the command line.
If a value is specified on the command line, that value will override a value specified in the inputs file.

See the Inputs section for more detail about what parameters can be specified at run-time.

#. From within the build directory type:

   ::

      <executable-name> <inputs-name>

   ``<executable-name>`` is  ``Nyx3d.Linux.gnu.ex1`` if you built Nyx with GNU Make, or
   ``nyx_<name-of-build-directory>`` if you built Nyx with CMake.

   ``<inputs-name>`` for a small test problem is  ``inputs.32`` for the MiniSB example,
   and ``inputs.rt`` for the LyA example. Most executable directories have an ``inputs``
   for a larger problem, and an ``inputs.rt`` or ``inputs.regtest`` for regression-test
   sized problems.

   .. note::
      For certain HPC systems, you may want to have a run directory separate from your compile / build directory.

      In that case, copy the executable and inputs file from the build directory to your run directory on scratch.
      Runs starting from a ``binary_particle_file`` need an absolute path in <inputs-name> or a symlink in the run directory.
      Runs with heating-cooling must have access to the TREECOOL_middle file in the run directory, and ascent in-situ runs
      need access to ascent_actions.yaml.
      

#. You will notice that running the code generates directories that look like
   plt00000, plt00020, etc.,
   and chk00000, chk00020, etc. These are “plotfiles” and
   “checkpoint” files. The plotfiles are used for visualization,
   and the checkpoint files for restarting the code.

See the Visualization chapter for how to visualize these plotfiles.
