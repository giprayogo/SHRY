Welcome to SHRY's documentation
===============================

``SHRY`` (**S**\ uite for **H**\ igh-th\ **r**\ oughput generation of models
with atomic substitutions implemented by p\ **y**\ thon)
is a tool for generating unique ordered structures
corresponding to a given disordered structure.

============
Installation
============

SHRY is available via PyPI (pending)

.. code-block:: console

    pip install shry

-----------
Development
-----------

Installing from source

.. code-block:: console

    git clone https://github.com/giprayogo/SHRY.git
    cd SHRY
    pip install .

==========
How to use
==========

If you have a disorder CIF file (i.e. containing partial occupations.),
simply use it as the first argument

.. code-block:: console

    shry STRUCTURE_CIF

there can be millions of order structure,
so it is advisable to start with a dry run

.. code-block:: console

    shry --count-only STRUCTURE_CIF

If you don't have a disorder structure,
and only order structure is there,
then choose one of the below
way to make the disorder CIF

---------------------------------------
Editing CIF file
---------------------------------------

Crystallographic Information File (CIF) file
is typically supplied in below format

.. code-block::

    loop_
    _atom_site_label
    _atom_site_type_symbol
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    _atom_site_occupancy
    _atom_site_U_iso_or_equiv
    Sm1 Sm 0.000 0.00 0.00 1.000 0.0
    Fe1 Fe 0.250 0.25 0.25 1.000 0.0
    Fe2 Fe 0.278 0.50 0.00 1.000 0.0

Suppose that here we want to replace Fe1 with
40/60 mix of Nb. Then do as below

.. code-block::

    loop_
    _atom_site_label
    _atom_site_type_symbol
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    _atom_site_occupancy
    _atom_site_U_iso_or_equiv
    Sm1 Sm 0.000 0.00 0.00 1.000 0.0
    Fe1 Fe 0.250 0.25 0.25 0.400 0.0
    Nb1 Nb 0.250 0.25 0.25 0.600 0.0
    Fe2 Fe 0.278 0.50 0.00 1.000 0.0

SHRY specifically read ``_atom_site_label``
to separate sites, so be careful.

----------------------------------------------
Using command line argument to create disorder
----------------------------------------------

Using the same example, you can do the same in command line

.. code-block::

    shry [-f/--from-species] Fe1 [-t/--to-species] [composition B] STRUCTURE_CIF

Which is equivalent to the above example.
Note that if instead ``Fe`` is supplied, then
all ``Fe`` sites are going to be replaced with single
``Fe0.4Nb0.6``, effectively merging ``Fe1`` sites and ``Fe2``
together, which may or may not be what you want.

-------------------------
More command line options
-------------------------

^^^^^^^^^^
Count only
^^^^^^^^^^

Before performing the whole run, it is a good
idea to check how many structures are going to be written.
By using

.. code-block::

    shry --count-only

SHRY will perform Polya enumeration and display how many there are.

^^^^^^^^^^^^^^^^^^
Creating supercell
^^^^^^^^^^^^^^^^^^

Sometimes the specified concentration needs a larger
supercell to fit.
You can enlarge the unit cell in command line by using

.. code-block::

    shry [-s/--scaling-matrix] [scaling_matrix]

In either 3-digit (diagonal) or 9-digit (non-diagonal) format.
For example a 2x2x1 supercell can be specified by either

.. code-block::

    shry -s 2 2 1 ...

or

.. code-block::

    shry -s 2 0 0 0 2 0 0 0 1 ...

^^^^^^^^^^^^^
Disorder only
^^^^^^^^^^^^^

If you just want to make the disorder structure,
you can add

.. code-block::

    shry --mod-only ...

in which case SHRY will output the disorder CIF without the order equivalents.

^^^^^^^^^^^^^
Other options
^^^^^^^^^^^^^

Other options can be found from the help menu

.. code-block::

    shry -h