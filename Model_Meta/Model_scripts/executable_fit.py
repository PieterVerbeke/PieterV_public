import Fitting_function as Fit
import numpy as np
from   multiprocessing import Pool

"""
Before running remember to change dataset variable in Fitting_function file
"""

cp = 18

worker_pool = []

with Pool(cp) as pool:
    result = pool.map(Fit.Fitting_execution, np.arange(cp))
