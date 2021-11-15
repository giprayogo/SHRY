# README.md

SHRY (**S**uite for **H**igh-th**r**oughput generation of models
with atomic substitutions implemented by p**y**thon)
is a tool for generating unique ordered structures
from a given disordered structure,
for use in computational materials simulations.
SHRY uses `pymatgen.Structure` extensively
for easy interfacing with Pymatgen.

## Installation and requirements

Clone, and

```console
pip install .
```

Requirements:

- `numpy` >= 1.17
- `pymatgen` >= 2019.1.1
- `scipy` >= 1.4.1
- `sympy` >= 1.5.1
- `tqdm`

## How to use

### Basic

Quick use: `shry STRUCTURE_CIF`.
See `shry -h` for more options.

### Examples

See `SHRY_INSTALLDIR/examples` for examples
on how to use SHRY as a Python module.
All examples use `SmFe12.cif` as the base CIF.

#### Site substitutions

The sample CIF has 4 Wyckoff positions.
Suppose that we would like to replace
the Fe sites, Fe1, Fe2, and Fe3, into Fe3Ti each.
With the command line interface:

```console
shry SmFe12.cif -f Fe1 Fe2 Fe3 -t Fe3Ti Fe3Ti Fe3Ti
```

Or as a python script (`example1.py`)

```python
import os

from shry import LabeledStructure, Substitutor

cif_file = 'SmFe12.cif'
structure = LabeledStructure.from_file(cif_file)
structure.replace_species({'Fe1': 'Fe3Ti'})
structure.replace_species({'Fe2': 'Fe3Ti', 'Fe3': 'Fe3Ti'})

substitutor = Substitutor(structure)
substitutor.make_patterns()
os.makedirs("output", exist_ok=True)
# A generator for Pymatgen's CifWriters
for i, cifwriter in enumerate(substitutor.cifwriters()):
    # Some naming logic.
    output_filename = f"output/{i}.cif"

    cifwriter.write_file(filename=output_filename)

```

`LabeledStructure` is basically Pymatgen's `Structure`,
but adds `_atom_site_label` as `site_properties`
with `replace_species` modified accordingly.
This allows for grouping the sites based on their label
when using supercell that is not consistent with some symmetry.
Otherwise, it is equivalent to `pymatgen.Structure` or `pymatgen.Molecule`.
See the last example for non-label-based grouping.

`Substitutor.make_patterns()` is the one that generates
the disordered structures.

#### Merging two or more positions

Here we would like to instead substitute all Fe sites into Fe3Ti.
In addition to unique structures equivalent to the one
given by the first example, this also includes cases
where the individual positions can be FeTi, FeTi3, etc.,
as long as the overall concentration adds up to Fe3Ti.
By command line:

```console
shry SmFe12.cif -f Fe -t Fe3Ti
```

As a python script (`examples/example2.py`)

```python
...

cif_file = 'SmFe12.cif'
structure = LabeledStructure.from_file(cif_file)
structure.replace_species({'Fe': 'Fe3Ti'})
...
```

#### Random sampling

Randomly choose 50 ordered structures from the output:

```console
shry SmFe12.cif -f Fe -t Fe3Ti --sample 50
```

or simply add `sample=50` (`examples/example3.py`)

```python
...
substitutor = Substitutor(structure, sample=50)
...
```

Note that each run will yield a different set of structures.

#### Using supercells

(Both are equivalent)

```console
shry SmFe12.cif -f Fe -t Fe7Ti -s 1 2 1
shry SmFe12.cif -f Fe -t Fe7Ti -s 1 0 0 0 2 0 0 0 1
```

The interface for making supercells of `LabeledStructure`
is the same as `pymatgen.Structure` (`examples/example4.py`)

```python
...
structure *= [1, 2, 1]
# 9-digits case
# structure *= [[1, 0, 0], [0, 2, 0], [0, 0, 1]]
...
```

#### Multiple target concentrations

It is preferable to reuse the same `Substitutor`
instance for all pattern generations.
This is because the algorithm works by recursively generating
substitution patterns for higher concentration from the lower ones.
Moreover, if two distinct systems have identical set of
permutations within the targeted sites, `Subsitutor`
can identify this and automatically re-arrange the site indices
to create the new substitutions from the previously generated patterns.

For example, consider substitution case of one of the Fe
sites on SmFe12, into three concentrations: FeTi, Fe3Ti, and FeTi3.
(`examples/example5.py`)

```python
...
structure1 = structure.copy()
structure2 = structure.copy()
structure3 = structure.copy()
structure1.replace_species({'Fe1': 'Fe7Ti'})
# Higher Ti concentration
structure2.replace_species({'Fe2': 'Fe3Ti'})
# Reverse of above
structure3.replace_species({'Fe2': 'FeTi3'})

os.makedirs("output1", exist_ok=True)
os.makedirs("output2", exist_ok=True)
os.makedirs("output3", exist_ok=True)

substitutor = Substitutor(structure1)
substitutor.make_patterns()
for i, cifwriter in enumerate(substitutor.cifwriters()):
    output_filename = f"output1/{i}.cif"
    cifwriter.write_file(filename=output_filename)

substitutor.structure = structure2
substitutor.make_patterns()
for i, cifwriter in enumerate(substitutor.cifwriters()):
    output_filename = f"output2/{i}.cif"
    cifwriter.write_file(filename=output_filename)

substitutor.structure = structure3
substitutor.make_patterns()
for i, cifwriter in enumerate(substitutor.cifwriters()):
    ...
```

For larger projects, it is perhaps useful to pickle
`PatternMaker` instances for later use.
See `examples/example6a.py` (save) and `examples/example6b.py` (load).

#### Ewald summation

Uses Pymatgen's `EwaldSummation` (see `examples/example7.py`).
Note that all sites must have their oxidation states defined.

#### Custom grouping

By default, SHRY groups sites by their labels within the CIF files
(typically corresponding to their crystallographic orbit).
Suppose that instead we would like to group them by their species label,
we can define a custom `groupby` as in `examples/example8.py`.

```python
# Some other processing
...
structure = ... # Structure, Molecule, or SiteCollection

# Loops over Structure.sites
def groupby(site):
   return site.species

substitutor = Substitutor(structure, groupby=groupby)

...
```

`groupby` loops over all sites within the `Structure`-like object,
meaning that `Site`/`PeriodicSite` are the input,
with group identifier (i.e. constant within the group but unique
to the group) as the output.
Anything hashable is valid as the identifier.
