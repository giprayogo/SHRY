How to use
==========

Quick use

.. code-block:: console

    shry STRUCTURE_CIF

You can prepare CIFs with partial occupations
by one of the ways below

---------------------------------------
Editing CIF file
---------------------------------------

The important part is the ``_atom_site_occupancy``
and ``_atom_site_label``, which are typically
grouped together in a loop

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

Suppose that here we want to replace ``Fe1``
to a 40/60 mix together with Nb.
Copy and edit the ``Fe1`` line, adjusting
the labels and occupations.

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

``SHRY`` will automatically stop if the total
occupancy of a site is either less or more than 1.0.
To simulate vacancies, create a pseudatom with species ``X``.

----------------------------------------------
Using SHRY option to create partial occupation
----------------------------------------------

The below example achieves the same modification.
You can also choose to above by the below options.

.. code-block::

    shry [-f/--from-species] Fe1 [-t/--to-species] Fe0.4Nb0.6 STRUCTURE_CIF

Note that ``SHRY`` targets either ``_atom_site_label`` or ``_atom_site_label``.
If instead ``Fe`` is used in the first argument,
all iron sites including ``Fe2`` will be replaced by ```Fe0.4Nb0.6``.

--------------------------------------------
Check total symmetry-inequivalent structures
--------------------------------------------

.. code-block:: console

    shry --count-only STRUCTURE_CIF

This operation is based on Polya enumeration
and takes much less time than a proper generation.

-------------------------
More command line options
-------------------------

^^^^^^^^^^^^^^^^^^
Creating supercell
^^^^^^^^^^^^^^^^^^

Sometimes a supercell is required to fit in finer concentrations.
``SHRY`` accepts either 3-digit (diagonal) or 9-digit (non-diagonal) format
to specify the supercell's scaling matrix.
For example a 2x2x1 supercell can be specified by either

.. code-block::

    shry -s 2 2 1 ...

or

.. code-block::

    shry -s 2 0 0 0 2 0 0 0 1 ...

^^^^^^^^^^^^^
Disorder only
^^^^^^^^^^^^^

If you just want to modify the CIF, without making the unique structures,
you can add

.. code-block::

    shry --mod-only ...

^^^^^^^^^^^^^
Other options
^^^^^^^^^^^^^

Other options can be found in the help menu

.. code-block::

    shry -h