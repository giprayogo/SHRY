# check cif files
python bench_symmetry_check.py

# split the xsl files
python bench_split_xsl.py

# run the benchmark for each SG
./bench_run.sh

# combine the results
python bench_combine_xsl.py