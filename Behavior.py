import Atom
import numpy as np
from scipy import stats


class Behavior:

    def __init__(self, name, vals=None, num_bins=4):
        """
        a Behavior represents one variable being analyzed for Genomic discretization
        :param name: name of the behavior
        :param vals: numpy array of all values for this behavior in the data set being analyzed for
        """
        self.num_bins = num_bins
        self.name = name
        self.vals = vals[~np.isnan(vals)] if vals is not None else None
        self.atoms = None

    def remove_zscore(self, zscore=2):
        """

        :param zscore: zscore threshold of values to remove
        :return:
        """

        x = self.vals
        x = stats.zscore(x)

        for i, val in enumerate(self.vals):
            if np.abs(x[i]) > zscore:
                self.vals[i] = np.NaN
        self.vals = self.vals[~np.isnan(self.vals)]

    def create_atom_bins(self):
        """
        creates bins for atoms based on min and max value of the behavior

        """

        # create numbins atoms
        self.atoms = [Atom.Atom() for _ in range(0, self.num_bins)]

        # find interval of each bin based on min and max
        b_min = np.min(self.vals) - 0.0001
        b_max = np.max(self.vals) + 0.0001
        interval = (b_max - b_min)/self.num_bins

        # iterate over atoms and assign upper and lower bounds
        count = 0
        for atom in self.atoms:
            atom.set_lower_bound(b_min + interval * count)
            count += 1
            atom.set_upper_bound(b_min + interval * count)

    def set_vals(self, vals):
        self.vals = vals

