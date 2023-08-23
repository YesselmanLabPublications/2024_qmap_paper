#!/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()

readme = open("README.md").read()

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="rna_secstruct",
    version="0.1.0",
    description="analysis for q_dms_ttr_paper",
    long_description=readme,
    long_description_content_type="test/markdown",
    author="Joe Yesselman",
    author_email="jyesselm@unl.edu",
    url="https://github.com/jyesselm/q_dms_ttr_paper",
    packages=[
        "q_dms_ttr_paper",
    ],
    package_dir={"q_dms_ttr_paper": "q_dms_ttr_paper"},
    py_modules=[
        "q_dms_ttr_paper/cli",
        "q_dms_ttr_paper/paths",
        "q_dms_ttr_paper/data_processing",
        "q_dms_ttr_paper/logger",
    ],
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords="rna_secstruct",
    classifiers=[
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: Implementation :: PyPy",
    ],
    entry_points={"console_scripts": []},
)
