#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from os import path
from setuptools import setup, find_packages

# Read the version from file
with open(
    path.join(path.dirname(__file__), "fmri_preproc", "resources", "version.txt")
) as fid:
    version = fid.read().strip()

# Read the contents of the README file
with open(
    path.join(path.abspath(path.dirname(__file__)), "README.md"), encoding="utf-8"
) as fid:
    readme = fid.read()

if __name__ == "__main__":
    setup(
        name="fmri_preproc",
        version=version,
        packages=find_packages(exclude=["tests"]),
        classifiers=[
            "Intended Audience :: Science/Research",
            "Natural Language :: English",
            "Operating System :: MacOS",
            "Operating System :: Linux",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3 :: Only",
            "Topic :: System :: Software Distribution",
            "Topic :: Utilities",
        ],
        python_requires=">=3.7",
        setup_requires=["pytest-runner"],
        tests_require=["pytest", "pytest-cov", "coverage"],
        scripts=["fmri_preproc.py"],
        keywords="Neonatal rs-fMRI preprocessing neuroimaging",
        description="Performs neonatal neuroimaging preprocessing of rs-fMRI data.",
        long_description=readme,
        long_description_content_type="text/markdown",
        author="Adebayo Braimah",
        author_email="adebayo.braimah@gmail.com",
        include_package_data=True,
        url="https://github.com/AdebayoBraimah/fmri_preproc",
    )
