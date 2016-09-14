source activate py27
# PYTHONPATH=../. py.test -m "current" --color=yes -svv -f tests/triform/test_merge_peaks.py
PYTHONPATH=../. py.test -m "current" --color=yes -svv -f tests/test_matrix/test_create_matrix.py
