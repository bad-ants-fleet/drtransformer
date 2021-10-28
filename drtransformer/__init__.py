#
# /drtransformer/__init__.py
#
# Distributed under the same license as the ViennaRNA package:
# https://github.com/ViennaRNA/ViennaRNA
#
# See README for instructions.
#

__version__ = "0.9"

_MIN_VRNA_VERSION = "2.4.13"

import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

from packaging import version

import RNA
if version.parse(RNA.__version__) < version.parse(_MIN_VRNA_VERSION):
    raise ImportError(f'ViennaRNA v{_MIN_VRNA_VERSION} required (found ViennaRNA v{RNA.__version__})!')

