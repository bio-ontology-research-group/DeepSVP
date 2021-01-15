#!/usr/bin/env python3
import os
import sys
import setuptools.command.egg_info as egg_info_cmd
from setuptools import setup

SETUP_DIR = os.path.dirname(__file__)
README = os.path.join(SETUP_DIR, "README.md")

try:
    import gittaggers
    tagger = gittaggers.EggInfoFromGit
except ImportError:
    tagger = egg_info_cmd.egg_info

install_requires = [
    "click < 8", "scikit-learn", "pandas", "tensorflow==2.3.0",
    "numpy", "scipy" , 
    "gensim==3.8.3",
    "sklearn",
    "mygene",
    "h5py",
    "progress",
]
needs_pytest = {"pytest", "test", "ptr"}.intersection(sys.argv)
pytest_runner = ["pytest < 6", "pytest-runner < 5"] if needs_pytest else []

setup(
    name="deepsvp",
    version="1.0.0",
    description="DeepSVP: Integration of Genomics and Phenotypes forStructural Variant Prioritization using Deep Learning",
    long_description=open(README).read(),
    long_description_content_type="text/markdown",
    author="Azza Althagafi",
    author_email="azza.althagafi@kaust.edu.sa",
    download_url="https://github.com/bio-ontology-research-group/deepsvp/archive/v1.0.0.tar.gz",
    license="Apache 2.0",
    packages=["deepsvp",],
    package_data={"deepsvp": [],},
    install_requires=install_requires,
    extras_require={},
    setup_requires=[] + pytest_runner,
    tests_require=["pytest<5"],
    entry_points={
        "console_scripts": [
            "deepsvp=deepsvp.main:main",
        ]
    },
    zip_safe=True,
    cmdclass={"egg_info": tagger},
    python_requires=">3",
)
