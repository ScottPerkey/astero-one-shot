from __init__ import TESTDATA

import os
from nonad import plot_modes
import numpy as np
import pandas as pd

def test_plot_modes():
    plot_modes.main(profile=TESTDATA + 'profile3.data', plot=True)

    # require that the dnu_scal, dnu_obs, and numax found all match to within absolute(a - b) <= (atol + rtol * absolute(b))
    # there is only one profile that is anallyzed, so just pull that row using .iloc[0]
    out = pd.read_csv(TESTDATA + 'profile3.data.astero.in', sep='\s+').iloc[0].values
    out_true = pd.read_csv(TESTDATA + 'profile3.data.astero.in.truth', sep='\s+').iloc[0].values

    # also check that a plot was made --- might not be the same byte-by-byte depending on local plotting parameters, so not checking content of the plot, just that one was made
    plot_made = os.path.exists(TESTDATA + 'profile3.data_freqs.png')
    # now remove the files that were made just in case need to do a test again.

    os.remove(TESTDATA + 'profile3.data_freqs.png')
    os.remove(TESTDATA + 'profile3.data.astero.in')
    assert (np.allclose(out, out_true, rtol=1e-05, atol=1e-08))*(plot_made)
    
if __name__ == '__main__':
    test_plot_modes()

