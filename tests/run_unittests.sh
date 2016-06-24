source activate py27
PYTHONPATH=../. py.test -m "unit" --color=yes -svv -f tests
