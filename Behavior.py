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
        print self.name
        print self.vals
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
        # behavior class no longer needs the values
        np.delete(self.vals, 1)

    def set_vals(self, vals):
        self.vals = vals

    def add_vals(self, vals):

        if self.vals is None:
            self.vals = np.array(vals)
        else:
            self.vals = np.append(self.vals, vals)

    def set_atom_codes(self, codes):
        for atom in self.atoms:
            atom.set_code(codes.pop()[0])

    def get_code_for_value(self, value):
        for atom in self.atoms:
            if atom.is_in_atom(value):
                return atom.get_code()
        return 'AAAAA'