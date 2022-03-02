# SmFeTiMo, benchmark test
SmFeTiMo, 1525023.cif (unsubstituted), $1\times1\times1\times$ supercell.

First you should set enviroment variables in ``bench.py``.

```python:bench.py
TIME_LOOP=3 # how many time you are going to measure
BABEL=False # True for supercell ver.1.2, False for that >= ver.2.0
SUPERCELL_PATH="your path to supercell BINARY"
SHRY_PATH = "your path to shry DIRECTORY"
SUPERCELL_FLAG=True # if false, supercell jobs are skipped.
SHRY_FLAG=False # if false, shry jobs are skipped.
```

Before the benchmark test, you should establish your python enviroment in which all the depenedent python packages are installed.

```bash:
conda activate shry
```

Finally, you can run the PbTe benchmark test

```bash:
python ../bench.py 1525023_partial.cif 1x1x1
```
