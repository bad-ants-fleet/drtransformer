# DrTransformer -- heuristic cotranscriptional folding.

DrTransformer (short for "DNA-to-RNA transformer") is a program for heuristic
and deterministic cotranscriptional folding simulations of RNA molecules. The
software uses the [ViennaRNA] package which is available through the [ViennaRNA
license].

## Installation
```sh
  ~$ python setup.py install
  ~$ python -m pytest tests/ -v -s
```

### ViennaRNA dependencies
This package uses the [ViennaRNA] library. Generally, every official release of
DrTransformer should be compatible with the latest [ViennaRNA bioconda] version.
For features available only on the development branches, you may have to install
a more recent version from the [ViennaRNA source] on github.

## Usage
Until further documentation is available, please use the *--help* options of the 
commandline executables:
```sh
  ~$ DrTransformer --help
  ~$ DrPlotter --help
```

## Version
v0.9 -- standalone package
  * extraction from the [ribolands] package to a standalone Python package.
  * using scipy and numpy for matrix exponentials (instead of [treekin])
  * implemented lookahead to skip pruning of potentially relevant future structures

## License
Same as the [ViennaRNA license]. 

## Cite
Badelt et al. (in preparation)
 
[//]: References
[ViennaRNA]: <http://www.tbi.univie.ac.at/RNA>
[ViennaRNA source]: <https://github.com/ViennaRNA/ViennaRNA>
[ViennaRNA bioconda]: <https://anaconda.org/bioconda/viennarna>
[ViennaRNA license]: <https://github.com/ViennaRNA/ViennaRNA/blob/master/license.txt>
[ribolands]: <https://github.com/bad-ants-fleet/ribolands>
[treekin]: <https://github.com/ViennaRNA/Treekin>

