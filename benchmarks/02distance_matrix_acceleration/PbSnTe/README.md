# PbSnTe, benchmark test
PbSnTe, 9011358.cif (unsubstituted), $2\times2\times2\times$ supercell.

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
python ../bench.py 9011358_partial.cif 2x2x2
```
