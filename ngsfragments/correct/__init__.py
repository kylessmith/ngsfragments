"""Init for mylibrary."""
from .correction import gaussian_smooth, correct
from .correct_intervals import calculate_interval_bias
from .cylowess.cylowess import lowess