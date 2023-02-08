## Version
v1.0 -- publication release  
  * switched from setup.py to pyproject.toml
  * added citation 

v0.12 -- perparing for official release
  * changed --t-lin, --t-log defaults and fixed --t-lin=1, --t-log=1
  * fixed potential issues with --t-end = --t-ext
  * adapted README example to publication 

v0.11 -- using lonely base-pairs
  * removed the --noLP default (added parameter setting)
  * added profiling option for runtime optimization
  * using --cg-auto default paramter
  * using k0=1e5, t-ext=0.04 default parameter
  * added new visulization types and fixed motif file input
  * added epsilon to t-fast sanity check

v0.10 -- moved to beta status (first official release)
  * changes in parameter defaults 
  * bugfix in linalg
  * new DrPlotter simulation layout and motif plotting
  * repaired code to enable plotting including pause sites

v0.9 -- standalone package (no official release)
  * extraction from the [ribolands] package to a standalone Python package.
  * using scipy and numpy for matrix exponentials (instead of treekin)
  * implemented lookahead to skip pruning of potentially relevant future structures

[//]: References
[ribolands]: <https://github.com/bad-ants-fleet/ribolands>
