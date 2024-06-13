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
    name="qmap_paper",
    version="1.0.0",
    description="analysis for qmap_paper https://www.biorxiv.org/content/10.1101/2024.03.11.584472v1",
    long_description=readme,
    long_description_content_type="test/markdown",
    author="Joe Yesselman",
    author_email="jyesselm@unl.edu",
    url="https://github.com/YesselmanLabPublications/2024_qmap_paper",
    packages=[
        "qmap_paper",
    ],
    package_dir={"qmap_paper": "qmap_paper"},
    py_modules=[
        "qmap_paper/cli",
        "qmap_paper/construct_design",
        "qmap_paper/data_processing",
        "qmap_paper/farfar_modeling",
        "qmap_paper/logger",
        "qmap_paper/mutation_characterize",
        "qmap_paper/paths",
        "qmap_paper/plotting",
        "qmap_paper/sasa",
        "qmap_paper/titration",
    ],
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords="qmap_paper",
    classifiers=[
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: Implementation :: PyPy",
    ],
    entry_points={"console_scripts": []},
)
