source activate py27
PYTHONPATH=../. py.test -m "unit" --color=yes -svv -n 10 -f tests
