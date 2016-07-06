source activate py27
PYTHONPATH=../. py.test -m "current" --color=yes -svv -f tests/triform/test_find_peaks.py
