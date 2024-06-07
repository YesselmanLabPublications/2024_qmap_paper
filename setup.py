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
    version="0.1.0",
    description="analysis for qmap_paper https://www.biorxiv.org/content/10.1101/2024.03.11.584472v1",
    long_description=readme,
    long_description_content_type="test/markdown",
    author="Joe Yesselman",
    author_email="jyesselm@unl.edu",
    url="https://github.com/jyesselm/qmap_paper",
    packages=[
        "qmap_paper",
    ],
    package_dir={"qmap_paper": "qmap_paper"},
    py_modules=[
        "qmap_paper/cli",
        "qmap_paper/data_processing",
        "qmap_paper/logger",
        "qmap_paper/paths",
        "qmap_paper/plotting",
        "qmap_paper/titration",
    ],
    include_package_data=True,
    # TODO uncomment when ready to publish
    # install_requires=requirements,
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
