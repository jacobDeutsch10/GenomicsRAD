import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import Atom


class RADFramework:

    def __init__(self):

        self.df = None
        self.subunit_categories = ['vv', 'thetaToR', 'zt', 'curvature']
        self.numBins = 4
        self.numAtoms = (len(self.subunit_categories)*self.numBins)
        self.atomList = []*self.numAtoms
        self.atom_dict = None
        self.rad_df = None

    def set_subunits(self, categories):

        self.subunit_categories = categories

    def read_table_from_csv(self, csv):
        self.df = pd.read_excel(csv, sheet_name=0, header=3)

    def drop_columns(self):

        columns_we_want = ['dt']

        # find column headers that we want
        for header in self.df.columns:
            # temporary fix should find better way to sort this out
            if self.convert_header_to_sub(header) in self.subunit_categories:
                columns_we_want.append(header)

        self.df.drop(self.df.columns.difference(columns_we_want), 1, inplace=True)

        print self.df.to_string()

    def create_histograms(self):
        for header in self.df.columns[1:]:
            x = self.df[header]
            x = x[~np.isnan(x)]
            plt.hist(x, bins=4)
            plt.ylabel(header)
            plt.show()

    @staticmethod
    def find_complement(code):
        complement = ''
        for i in code:
            if i == 'A':
                complement += 'T'
            elif i == 'T':
                complement += 'A'
            elif i == 'C':
                complement += 'G'
            elif i == 'G':
                complement += 'C'
        return complement

    def create_atom_codes(self):
        dna = ["A", "G", "C", "T"]

        atoms = []

        for i in dna:
            for j in dna:
                for k in dna:
                    atom = i+j+k
                    complement = self.find_complement(atom)
                    if atom not in atoms and complement not in atoms:
                          atoms.append(atom)
        # random.shuffle(atoms)
        complements = [self.find_complement(i) for i in atoms]

        self.atomList = [(atoms[i], complements[i]) for i in range(0, self.numAtoms)]

        print self.atomList

    def assign_atoms_bins(self):

        self.atom_dict = dict.fromkeys(self.subunit_categories)
        for sub in self.subunit_categories:
            self.atom_dict[sub] = [Atom.Atom(self.atomList.pop()[0]) for i in range(0, self.numBins)]
        for key in self.atom_dict:
            sub_min = self.get_subunit_min(key)-0.0001
            sub_max = self.get_subunit_max(key)
            interval = (sub_max - sub_min)/self.numBins
            count = 0
            for atom in self.atom_dict[key]:

                atom.set_lower_bound(sub_min + interval*count)
                count += 1
                atom.set_upper_bound(sub_min + interval*count)
                print "%s %s lower: %f upper: %f" % (key, atom.get_code(), atom.get_lower_bound(), atom.get_upper_bound())

    def get_subunit_min(self, sub):
        column_names = []
        for col in self.df.columns:
            name = self.convert_header_to_sub(col)
            if name == sub:
                column_names.append(col)
        print column_names
        return min(self.df[column_names].min())

    def get_subunit_max(self, sub):
        column_names = []
        for col in self.df.columns:
            name = self.convert_header_to_sub(col)
            if name == sub:
                column_names.append(col)
        print column_names
        return max(self.df[column_names].max())

    @staticmethod
    def convert_header_to_sub(name):
        return '_'.join(name.split("_")[:-1])

    def create_rad_matrix(self):

        columns = []

        for col in self.df.columns[1:]:
            for i in range(0, self.numBins):
                columns.append(col + '_' + str(i))
        print columns

        self.rad_df = pd.DataFrame()
        self.rad_df['dt'] = self.df['dt']
        for col in columns:
            eq_col = self.convert_header_to_sub(col)
            num = int(col.split("_")[-1])
            sub_unit = self.convert_header_to_sub(eq_col)
            self.rad_df[col] = [self.atom_dict[sub_unit][num].get_code() if (self.atom_dict[sub_unit][num]).is_in_atom(i)
                           else (self.atom_dict[sub_unit][num]).get_complement() for i in self.df[eq_col]]

        print self.rad_df.to_string()

    def average_df_over_time_step(self, step):

        avg_dict = dict.fromkeys(self.df.columns)
        num_rows = self.df.shape[0]-1
        print num_rows
        for key in avg_dict:
            avg_dict[key] = []
        count = 1
        first = 0
        for index, row in self.df.iterrows():
            if index == num_rows:

                for key in avg_dict.keys():
                    if key != 'dt':
                        avg_dict[key].append(np.mean(self.df[key][first:index + 1]))
                avg_dict['dt'].append(step * count)
            elif pd.isnull(self.df['dt'][index]):
                continue
            elif self.df['dt'][index] - self.df['dt'][first] < step:
                continue
            else:
                for key in avg_dict.keys():
                    if key != 'dt':
                        avg_dict[key].append(np.mean(self.df[key][first:index+1]))

                avg_dict['dt'].append(step*count)
                count += 1
                first = index
        temp = pd.DataFrame(columns=self.df.columns)
        for col in temp:
            temp[col] = avg_dict[col]
        self.df = temp
        # drop rows with >= 2 Nan values
        self.df = self.df.dropna(thresh=self.df.shape[1]-2)
        print self.df.to_string()