# Using SHRY in Python

This directory contains examples on
how to use SHRY as a Python module.

## Basic function

(See `example1.py` for this section).

`Substitutor` is the main interface for using
functions implemented in SHRY.
It uses Pymatgen's `Structure` as the
structure representation.

```python
from pymatgen.core import Structure
from shry import Substitutor

...

structure = Structure.from_file(file_name)
substitutor = Substitutor(structure)
```

We recommend enumerating the unique structures
using `Substitutor.count()` before generating them.
It is an implementation of Polya enumeration,
which is almost instant for most cases.

The structures (either as CIF or `Structure`), weights,
configuration letters, etc. can then be obtained from
`Substitutor.quantities(string_tuple)`.
The `string_tuple` may contain any of these keywords:

- `cifwriter`. Pymatgen's `CifWriter` instances.
- `structure`. Pymatgen's `Structure` instances.
  Use this for writing CIF files.
- `weight`. How many configurations
- `letter`. Configuration letter ('aaa', 'bab', etc.)
  corresponding to the substitution.
- `ewald`. Ewald energy for the given structure.

`Substitutor.quantities(string_tuple)` is a generator
of a dictionary with the previous keywords as keys.
For example, if you want to get the CIFs and weights, do

```python
for i, packet in enumerate(substitutor.quantities(("cifwriter", "weight"))):
    cifwriter = packet["cifwriter"]
    weight = packet["weight"]

    filename=f"cif_i{i}w{weight}.cif"
    cifwriter.write_file(filename=os.path.join(output_dir, filename))
```

There are also individual generators for each of the quantities:

- `cifwriters()`
- `structure_writers()`
- `weights()`
- `letters()`
- `ewalds()`

however these will invoke one full run for every call,
so if more than one quantities is required,
it will be slower.

## Enumlib equivalent

See `example2.py` for comparison with equivalent
`enumlib` functions through Pymatgen's `EnumlibAdaptor`.

## Advanced use

### LabeledStructure

`LabeledStructure` is a modified Pymatgen's `Structure`,
but it tracks the CIF's `_atom_site_label`.
This is useful if you want to group sites
together regardless of supercell use,
which can sometimes split the sites.

It also implements `replace_species` for substituting sites
with a slightly more convenient syntax (see `example3.py`).

### Saving substitutor instance

`Substitutor` can automatically "remap" pattern
generated in other structure if the structures are
symmetrically similar.
This can save a lot of time if you are dealing with
multiple concentrations or a set of symmetrically similar sytems.

See `example4a.py` for how to enable caching and pickle
the `Substitutor` instance, and `example4b.py` for the later reloading.
