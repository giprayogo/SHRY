# README.md

SHRY (**S**uite for **H**igh-th**r**oughput generation of models
with atomic substitutions implemented by p**y**thon)
is a tool for generating symmetrically unique order structures,
from a given disordered structure, ready for computational
materials simulations.
Uses `pymatgen.Structure` extensively for easy interface with Pymatgen.

## Installation and requirements

Available on PyPI (Soon! Use `pip install -e .` for now).

```console
pip install shry
```

SHRY requires:

- `numpy` >= 1.17
- `pymatgen` >= 2019.1.1
- `scipy` >= 1.4.1
- `sympy` >= 1.5.1
- `tqdm`

## How to use

### Command line

```console
shry input.cif [-f FROM_SPECIES] [-t TO_SPECIES] [-s SUPERCELL] [--sample SAMPLE] [--symprec SYMPREC] [--dir-size DIR_SIZE] [--count-only] [--no-write-symm]
```

Arguments:

- `input`. CIF structure file, or an `*.ini` file containing
  the rest of command line arguments. If using an `*.ini` file,
  write all arguments under `DEFAULT` section (See `./examples`).

Optional arguments:

- `--from-species`, `-f`. Replace species/label from the input
  structure into `--to` species.
  Matches either `_atom_site_label` or `_atom_site_type_symbol`
  as defined within the `input` CIF file.
- `--to-species`, `-t`. Final chemical formula of the `--from` sites
  after substitution. Recognizes decimal as well as oxidations states,
  thus one may write, for example, Sm2+0.5Sm3+Fe10.5, etc.
  Different oxidation states are treated as symmetrically distinct sites.
- `--supercell`, `-s`. Rescales the unit cell. Accepts comma or space separated
  3 or 9 digits values.
- `--sample`. Sample `SAMPLE` from the generated structures (no replacement).
- `--symprec`. Precision used in symmetry searches.
- `--dir-size`. Divides output structures into several directories
  of `DIR_SIZE` files each. Default to 10000.
- `--count-only`. Skip structures generation, instead only prints
  the final number of structure (from Polya enumeration).
- `--write-symm`. Write symmetries for the generated structures (slower).

Also available on the help menu:

```console
shry -h
```

Module-based usages are given in the examples.

### Examples

See `examples` for the input files.

#### Fe substitution on SmFe12

Using `examples/SmFe12.cif` as the base structure,
suppose that we would like to substituted Fe1, Fe2, and Fe3,
separately into Fe3Ti each. Using the script:

```console
shry SmFe12.cif -f Fe1,Fe2,Fe3 -t Fe7Ti,Fe7Ti,Fe7Ti
```

which is almost equivalent to (`examples/example1.py`)

```python
from shry import LabeledStructure, Substitutor

cif_file = 'SmFe12.cif'
structure = LabeledStructure.from_file(cif_file)
structure.replace_species({'Fe1': 'Fe7Ti'})
structure.replace_species({'Fe2': 'Fe7Ti'})
structure.replace_species({'Fe3': 'Fe7Ti'})

substitutor = Substitutor(structure)
os.makedirs("output", exist_ok=True)
# A generator of Structure instances.
for i, _structure in enumerate(substitutor.order_structure()):
    # Some naming logic.
    output_filename = f"output/{i}.cif"

    _structure.to(output_filename)
```

`LabeledStructure` inherits Pymatgen's `Structure`,
but adds `_atom_site_label` into the `site_properties`.
This allows us to keep track of the relationship between
the `LabeledStructure`'s sites and those defined within the CIF file.
The overwritten `replace_species()` reads this property or the sites' element.

While `Subsitutor` accepts `pymatgen.Structure`
or `pymatgen.Molecule` instances just fine,
each distinct orbits of the disorder sites will be treated
as separate groups.
Beware that using a non-conventional cell may result in unexpected orbit.
If similar grouping is desired, either manually add the `_atom_site_label`,
or write a custom `groupby` for the `Substitutor` (see the last example).
Note that `Subsitutor` only treats disorder sites.

#### Grouping distinct orbits

Here we would like to instead substitute _all_ the Fe sites into Fe3Ti.
In addition to structures given by the first example,
this also includes cases where the individual orbits are not Fe3Ti,
but the overall sites adds up to Fe3Ti.
Therefore, more structures are generated.
By command line:

```console
shry SmFe12.cif -f Fe -t Fe3Ti
```

Or with a script (`examples/example2.py`)

```python
...

cif_file = 'SmFeTi.cif'
structure = LabeledStructure.from_file(cif_file)
structure.replace_species({'Fe': 'Fe7Ti'})

substitutor = Substitutor(structure)
os.makedirs("output", exist_ok=True)
for i, _structure in enumerate(substitutor.order_structure()):
    ...
```

`replace_species()` automatically re-labeled the Fe sites into a single `Fe`,
so that `Substitutor` recognizes them as a single group.

#### Random sampling

To take 50 random samples from the output:

```console
shry SmFe12.cif -f Fe -t Fe7Ti --sample 50
```

or simply (`examples/example3.py`)

```python
...

substitutor.sample = 50
for _structure in substitutor.order_structure():
    ...
```

#### Using supercells

(Both are equivalent)

```console
shry SmFe12.cif -f Fe -t Fe7Ti -s 1 2 1
shry SmFe12.cif -f Fe -t Fe7Ti -s 1 0 0 0 2 0 0 0 1
```

The interface for `LabeledStructure` is the same as `pymatgen.Structure`
(`examples/example4.py`)

```python
import numpy as np

...

structure *= [1, 2, 1]
# 9-digits case.
structure *= np.array([[1, 0, 0], [0, 2, 0], [0, 0, 1]])

...
```

#### Multiple target concentrations

`Substitutor` can re-use previously generated patterns
if the structure is symmetrically identical
to a previously used structures.
As an example, let us consider substituting the three Fe
sites on SmFe12, separately into FeTi, Fe2Ti, and Fe3Ti.
(`examples/example5.py`)

```python
...

cif_file = 'SmFe12.cif'
structure = LabeledStructure.from_file(cif_file)
structure *= [1, 2, 1]
structure1 = structure.copy()
structure2 = structure.copy()
structure3 = structure.copy()
structure4 = structure.copy()
structure1.replace_species({'Fe1': 'Fe7Ti1'})
structure2.replace_species({'Fe2': 'Fe6Ti2'})
structure3.replace_species({'Fe3': 'FeTi'})
structure4.replace_species({'Fe2': 'FeTi3'})

substitutor = Substitutor(structure1)
for _structure in substitutor.order_structure():
    ...

# Replace structure with structure2
subsitutor.structure = structure2
for _structure in substitutor.order_structure():
    ...

subsitutor.structure = structure3
for _structure in substitutor.order_structure():
    ...

# Equivalent to structure2 but with reversed Fe-Ti concentration
subsitutor.structure = structure4
for _structure in substitutor.order_structure():
    ...
```

The `structure2`, `structure3`, and `structure4` part re-used
the same instance of `PatternGenerator`.
Generation for Fe6Ti2 (Fe=75% Ti=25%)
continues from the generated patterns of Fe7Ti1 (Fe=86% Ti=14%).
Patterns of `structure4` are automatically re-labeled from those of
`structure2`, as they are symmetrically identical.

Even when working on different structures,
as long as they have the same symmetry,
`Substitutor` will automatically re-use or build upon previously
made patterns.
This works internally by keeping a collection of `PatternGenerator`
instances saved in a dictionary with their `signature()` as keys.
In larger projects, it is perhaps advantageous to pickle
`PatternGenerator` instances for later use (`examples/example6a.py`):

```python
import pickle

...

# Save pattern generators for later use.
with open('pg.pkl', 'wb') as f:
    pickle.dump(substitutor.pattern_generators, f)
```

Reload the `pattern_generators` before substituting

```python
...

# Load previous patterns
with open('pg.pkl', 'rb') as f:
    pattern_generators = pickle.load(f)
substitutor.pattern_generators = pattern_generators

# Substitute other structures
for _structure in substitutor.order_structure():
    ...
```

#### Ewald summation

Uses Pymatgen's `EwaldSummation` (see `examples/example7.py`).
Note that all sites must have their oxidation states defined.

#### Custom grouping

Suppose that we have processed `pymatgen.Structure` or `pymatgen.Molecule`
from some other place in our pipeline,
and we would like to replace sites by their atomic species (`examples/example8.py`).

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

Since `groupby` basically loops over the sites,
any custom grouping utilizing Pymatgen's `Site` / `PeriodicSite`
properties is possible.
