"""Init for mylibrary."""
from __future__ import absolute_import
from .fragments import Fragments
from . import utilities
from .correct import correction
from .segment import cnv
from .sequence import *
from .plot import plot_plt
from .io import read_sam
from .metrics import gene_activity

# This is extracted automatically by the top-level setup.py.
__version__ = '2.4.2'