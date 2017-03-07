import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

pandas2ri.activate()

r_base = None
r_sva = None


def combat(expr, batch):
    """Normalizes expression data for batch effects using ComBat.

    Parameters
    ----------
    expr : pandas.DataFrame
        Dataframe of samples-by-genes containing expression values.

    """

    assert expr.shape[1] == batch.shape[0]

    global r_base
    global r_sva

    # Perform imports if needed.
    if r_base is None:
        r_base = importr('base')

    if r_sva is None:
        r_sva = importr('sva')

    # Filter for minimum variance.
    expr_filt = expr.ix[expr.var(axis=1) > 1e-12]

    # Run ComBat.
    expr_combat = r_sva.ComBat(
        r_base.as_matrix(expr_filt), batch=r_base.as_character(batch))

    # Convert back to frame and return.
    return pd.DataFrame(
        np.array(expr_combat),
        index=expr_filt.index,
        columns=expr_filt.columns)
