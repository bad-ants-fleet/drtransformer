#!/usr/bin/env python
from setuptools import setup

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

setup(
    name = 'drtransformer',
    version = '0.10',
    description = 'Heuristic cotranscriptional folding using the nearest neighbor energy model.',
    long_description = LONG_DESCRIPTION,
    long_description_content_type = 'text/markdown',
    author = 'Stefan Badelt',
    author_email = 'bad-ants-fleet@posteo.eu',
    maintainer = 'Stefan Badelt',
    maintainer_email = 'bad-ants-fleet@posteo.eu',
    url = 'https://github.com/bad-ants-fleet/drtransformer',
    license = 'MIT',
    classifiers = [
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Intended Audience :: Science/Research',
        ],
    python_requires = '>=3.8',
    install_requires = [
        'numpy',
        'scipy',
        'matplotlib'],
    packages = ['drtransformer'],
    test_suite = 'tests',
    entry_points = {
        'console_scripts': [
            'DrTransformer=drtransformer.drtransformer:main',
            'DrPlotter=drtransformer.plotting:main',
            'DrFindpath=drtransformer.rnafolding:call_findpath_exe',
            'DrCoarseGrain=drtransformer.rnafolding:top_down_coarse_graining_exe',
            'DrSimulate=drtransformer.linalg:main',
            ],
        }
)

