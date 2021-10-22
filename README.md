# DrTransformer -- heuristic cotranscriptional folding.

DrTransformer (short for "DNA-to-RNA transformer") is a program for heuristic
and deterministic cotranscriptional folding simulations of RNA molecules.

## Installation
```sh
  ~$ python setup.py install
  ~$ python -m pytest tests/ -v -s
```

## Usage
Until further documentation is available, please use the *--help* options of the 
commandline executables:
```sh
  ~$ DrTransformer --help
  ~$ DrPlotter --help
```

## ViennaRNA dependencies
This package uses the [ViennaRNA] library. Generally, every official release of
DrTransformer should be compatible with the latest [ViennaRNA bioconda] version.
For features available only on the development branches, you may have to install
a more recent version from the [ViennaRNA source] on github.

## Cite
Badelt et al. (in preparation)
 
## Version
v0.9-dev -- standalone package
  * extraction from the [ribolands] package to a standalone Python package.
  * using scipy and numpy for matrix exponentials (instead of [treekin])

## License
MIT

[//]: References
[ViennaRNA]: <http://www.tbi.univie.ac.at/RNA>
[ViennaRNA source]: <https://github.com/ViennaRNA/ViennaRNA>
[ViennaRNA bioconda]: <https://anaconda.org/bioconda/viennarna>
[ribolands]: <https://github.com/bad-ants-fleet/ribolands>
[treekin]: <https://github.com/ViennaRNA/Treekin>

