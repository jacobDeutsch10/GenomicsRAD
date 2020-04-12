import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import Atom
from scipy import stats


class RADFrame:

    def __init__(self, behaviors=None, frame=None, num_bins=4):

        self.df = frame
        self.subunit_categories = ['vv', 'thetaToR', 'zt', 'curvature'] if behaviors is None else behaviors
        self.numBins = num_bins
        self.numAtoms = (len(self.subunit_categories)*self.numBins)
        self.atomList = []*self.numAtoms
        self.atom_dict = None
        self.rad_df = None
        self.filename = None

    def set_subunits(self, categories):

        self.subunit_categories = categories

    def set_df(self, df):

        self.df = df

    def read_table_from_xl(self, xlname):
        self.df = pd.read_excel(xlname, sheet_name=0, header=3)
        self.df = self.df.dropna(how='all')
        self.filename = xlname

    def get_behavior_vals(self, behavior):
        vals = np.empty(0)
        for col in self.df.columns:
            if self.convert_header_to_sub(col) == behavior:
                vals = np.append(vals, self.df[col].values)
        return vals[~np.isnan(vals)]

    def drop_columns(self):

        columns_we_want = ['dt']

        # find column headers that we want
        for header in self.df.columns:
            # temporary fix should find better way to sort this out
            if self.convert_header_to_sub(header) in self.subunit_categories:
                columns_we_want.append(header)

        self.df.drop(self.df.columns.difference(columns_we_want), 1, inplace=True)

        print self.df.to_string()

    def create_histograms(self, color='b'):
        for header in self.df.columns[1:]:
            x = self.df[header]
            x = x[~np.isnan(x)]
            print header
            print stats.zscore(x)
            plt.hist(x, color=color, bins=10)
            plt.ylabel(header)
            plt.show()

    @classmethod
    def find_complement(cls, code):
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

    def create_atom_codes(self, num=3):
        dna = ["A", "G", "C", "T"]
        first = 'A'*num
        atoms = [first]
        print atoms
        if num == 3:
            for i in dna:
                for j in dna[-1:] + dna[:-1]:
                    for k in dna[-2:] + dna[:-2]:
                        atom = i+j+k
                        complement = self.find_complement(atom)
                        if atom not in atoms and complement not in atoms:
                            atoms.append(atom)
        elif num == 4:
            for i in dna:
                for j in dna[-1:] + dna[:-1]:
                    for k in dna[-2:] + dna[:-2]:
                        for l in dna[-3:] + dna[:-3]:
                            atom = i+j+k+l
                            complement = self.find_complement(atom)
                            if atom not in atoms and complement not in atoms:
                                atoms.append(atom)
        elif num == 5:
            for i in dna:
                for j in dna[-1:] + dna[:-1]:
                    for k in dna[-2:] + dna[:-2]:
                        for l in dna[-3:] + dna[:-3]:
                            for m in dna[-4:]+ dna[:-4]:
                                atom = i+j+k+l+m
                                complement = self.find_complement(atom)
                                if atom not in atoms and complement not in atoms:
                                    atoms.append(atom)
        print atoms
        complements = [self.find_complement(i) for i in atoms]

        self.atomList = [(atoms[i], complements[i]) for i in range(0, self.numAtoms+1)]

    def assign_atoms_bins(self):

        self.atom_dict = dict.fromkeys(self.subunit_categories)
        start = Atom.Atom(self.atomList.pop(0)[0])
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
        self.atom_dict["start"] = start

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

    def remove_outliers_zscore(self, z_thresh=2):

        for col in self.df.columns[1:]:
            x = self.df[col]
            x = x[~np.isnan(x)]
            x = stats.zscore(x)
            print col
            print x
            count = 0
            for i, val in enumerate(self.df[col]):
                if np.isnan(val):
                    continue
                elif np.abs(x[count]) > z_thresh:
                    count += 1
                    self.df[col][i] = np.NaN
                else:
                    count += 1
        print self.df.to_string()

    def add_start_stop(self):

        self.rad_df = self.rad_df.drop("dt", axis=1)
        num_objects = len(self.rad_df.columns)/self.numAtoms
        rad_by_object = np.array_split(self.rad_df, num_objects, axis = 1)
        for rad in rad_by_object:
            rad.insert(0, 'start', self.atom_dict['start'].get_code())
            rad.insert(len(rad.columns), 'stop', self.atom_dict['start'].get_complement())

        for i in rad_by_object:
            print self.convert_rad_matrix_to_string(i)

    @staticmethod
    def convert_rad_matrix_to_string(df):

        radString = df.to_string(header=False, index=False, index_names=False).split('\n')

        radString = ''.join(radString)
        return ''.join(radString.split(' '))

    def get_key_from_xl(self):

        plotting_sheet = pd.read_excel(self.filename, sheet_name=2, header=0)

        data_name = plotting_sheet.columns[0].split('_')
        data_name = '_'.join([data_name[1], data_name[3]])
        plotting_sheet = plotting_sheet[18:]
        frames = ''
        object_names = []
        for index, row in plotting_sheet.iterrows():
            if row[0] == '*':

                frames = str(row[10]) + '_' + str(row[11])
                object_names = [row[1], row[2]]

        keys = [str(o) + '_' + data_name + '_' + frames for o in object_names]
        print keys
        return keys

    def __str__(self):
        return self.df.to_string()
