source activate py27
PYTHONPATH=../. py.test -m "current" --color=yes -svv -f tests/preprocess/test_make_chromosome_cover_files.py
