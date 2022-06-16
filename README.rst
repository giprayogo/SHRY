README
==========

.. image:: logo/logo.jpg
    :alt: logo

SHRY (\ **S**\ uite for \ **H**\ igh-th\ **r**\ oughput generation of models
with atomic substitutions implemented by p\ **y**\ thon)
is a tool for generating unique ordered structures
corresponding to a given disordered structure.

Installation
------------

SHRY can be obtained from PyPI

.. code-block:: console

    pip install shry

For Windows Users
^^^^^^^^^^^^^^^^^

For Windows users,
if you don't have Python already,
you can try, for example,
installing Python from the Microsoft store
following instructions on
`this page`_.

.. _`this page`: https://docs.microsoft.com/en-us/windows/python/beginners.

Then install SHRY just like above
within PowerShell or your favourite terminal.

.. code-block:: console

    pip install shry

Development
^^^^^^^^^^^

If you prefer to install from source,
instead follow the procedures below.

.. code-block:: console

    git clone https://github.com/giprayogo/SHRY.git
    cd SHRY
    pip install -e .

Quick use
---------

Preparation of an input file (a CIF file)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can prepare a CIF file with partial occupations.

.. code-block::

    # label element x y z occupation
    Sm1 Sm 0.000 0.00 0.00 1.000
    Fe1 Fe 0.250 0.25 0.25 0.400
    Nb1 Nb 0.250 0.25 0.25 0.600
    Fe2 Fe 0.278 0.50 0.00 1.000

``SHRY`` will automatically stop if the total occupancy of a site is
either less or more than 1.0. To simulate vacancies, create a pseudo
atom with species ``X``.

Check total symmetry-inequivalent structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can readily check the number of total symmetry-inquivalent
structures using the following command.

.. code-block:: console

    shry --count-only STRUCTURE_CIF

This operation is based on Polya enumeration and takes much less time
than a proper generation.

Creating supercell
^^^^^^^^^^^^^^^^^^

Sometimes a supercell is required to fit in finer concentrations.
``SHRY`` accepts either 3-digit (diagonal) or 9-digit (non-diagonal)
format to specify the supercell's scaling matrix. For example a 2x2x1
supercell can be specified by either

.. code-block:: console

    shry -s 2 2 1 --count-only STRUCTURE_CIF

or

.. code-block:: console

    shry -s 2 0 0 0 2 0 0 0 1 --count-only STRUCTURE_CIF

Generating unique structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, you can generate symmetry-inequivalent structures using the
following command:

.. code-block:: console

    shry -s 2 2 1 STRUCTURE_CIF

The generated symmetry-inequivalent structures are saved in sliceXX
directories.

Additional information
^^^^^^^^^^^^^^^^^^^^^^

For additional information, you can use the help command:

.. code-block:: console

    shry -h

or you can refer to the documentation.

Documentation
-------------

The documentation is available `here <https://shry.readthedocs.io/en/latest/>`_.

How to cite
-------------

When using our code please cite

- |DOI|
- |Paper|

.. |DOI| image:: https://zenodo.org/badge/425687455.svg
   :target: https://zenodo.org/badge/latestdoi/425687455

.. |Paper| image:: https://img.shields.io/static/v1?label=arXiV&message=2111.13409&color=b31b1b
   :target: https://arxiv.org/abs/2111.13409

Contributing to the project
---------------------------

If you want to contribute to the project, report a bug, or ask for
a new feature, please `raise an issue <https://github.com/giprayogo/SHRY/issues>`_.
