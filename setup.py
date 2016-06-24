import sys
from setuptools import setup, find_packages

from triform.version import __version__
install_requires = ["scipy", "pandas", "numpy", "natsort", "joblib", "pyfaidx"]

setup(
    name="triform",
    packages=find_packages(),
    package_dir={'triform': 'triform'},
    package_data={'triform': ['*.R', "scripts/chromsizes/*.chromsizes"]},
    # scripts=["bin/epic", "bin/epic-effective"],
    version=__version__,
    description=
    "Improved sensitivity, specificity and control of false discovery rates in ChIP-Seq peak finding.",
    author="Endre Bakken Stovner",
    author_email="endrebak85@gmail.com",
    url="http://github.com/endrebak/triform",
    keywords=["ChIP-Seq"],
    license=["MIT"],
    install_requires=install_requires,
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Development Status :: 2 - Pre-Alpha",
        "Environment :: Other Environment", "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Topic :: Scientific/Engineering"
    ],
    long_description=
    "Improved sensitivity, specificity and control of false discovery rates in ChIP-Seq peak finding.")
