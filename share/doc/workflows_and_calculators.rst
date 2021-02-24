Workflow and Calculator Structure
#################################

The main functions of XTP are distributed over three command line tools ``xtp_run`` and ``xtp_parallel`` which execute so called *calculators* and ``xtp_tools`` which execute *tools*. The main difference between the two is that calculators perform tasks on a topology (a state file) such as, QM/MM calculations, performing KMC simulations, computing transfer integrals between pairs of molecules etc. The tools on the other hand perform stand-alone tasks such as, calculations on a single molecule, preparing and converting input files and post processing of data.

There is another command line tool called ``xtp_map`` that is used to perform the specific task of converting a GROMACS or LAMMPS MD topology to an XTP state file. The tool is capable of partitioning a topology into conjugated segments (hopping sites) and rigid fragments. If desired the tool can also map groundstate geometries to the distorted MD geometries and it can map multipole and polarization data to atomic expansion sites.

Basic Usage of the Tools
~~~~~~~~~~~~~~~~~~~~~~~~

Tools are used in the following way

.. code-block:: bash

   xtp_tools -e <toolName> -n <basename> -o <options>.xml

the ``-e`` flag is used for what to execute e.g. ``dftgwbse``, ``orb2mol`` etc. The ``-n`` flag is used for the "basic" file on which the tool operates, it only makes sense in relation to the tool chosen, as an example, the basename for a DFT-GWBSE calculation will be the name of the .xyz file, without the file extension, while in the case of a file converter it will be the name of the file to convert etc. Options are passed to the tool with the help of an xml file. 

Basic Usage of the Calculators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To use one of the calculators you will first need to create an XTP state with the ``xtp_map`` tool. This will create a state in an hdf5 format, which can be inspected with an hdf5 viewer such as Silx. After the creation of the state a calculator is run (depending on the task with either ``xtp_run`` or ``xtp_parallel``) in the following way

.. code-block:: bash
   
   xtp_run -e <calculatorName> -f <state>.hdf5 -o <options>.xml

Further flags and options might be necessary and will be discussed in the documentation of the specific calculators and in the tutorials.

An Overview of the Tools and Calculators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   xtp_map_overview
   xtp_parallel_overview
   xtp_run_overview
   xtp_tools_overview
