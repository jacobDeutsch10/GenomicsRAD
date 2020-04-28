import RADFrame
import Atom
import Behavior
import pandas as pd
import numpy as np


class MultiFrame:
    """
    a MultiFrame holds a list of RADFrames for all animals in a dataset
    """

    def __init__(self, frames=None, behaviors=None, num_bins=4):
        """

        :param frames: a list of pandas frames representing the behavior of each animal
        :param behaviors: a list of variables to be used for Genomic conversion
        """

        default = ['vv', 'thetaToR', 'zt', 'curvature']
        self.frames = frames if frames is not None else []
        self.filename = None
        self.num_bins = num_bins
        if behaviors is None:
            self.behaviors = {key: Behavior.Behavior(name=key, num_bins=self.num_bins) for key in default}
        else:
            self.behaviors = {key: Behavior.Behavior(name=key, num_bins=self.num_bins) for key in behaviors}

        self.keys = []
        self.atomList = []

    def read_from_xl_FULL(self, filename):
        """
        Reads a MultiFrame in from an excel file with data from many animals
        :param filename: the name of the excel file
        """
        self.filename = filename
        # readfile and drop empty columns
        df = pd.read_excel(filename, index=False, sheet_name=0, header=2)
        df = df.dropna(how='all')
        index = 0
        frames = []
        print "DUDE plus df.iloc =" + df.iloc[0][0]
        # iterate over rows of excel file and separate into a list of dataframes
        while index < len(df):
            found = False
            if df.iloc[index][0] == 'dt':
                start = index
                found = True

                for i in range(start + 1, len(df)):
                    if df.iloc[i][0] == 'dt' or i == len(df)-1:
                        df2 = df[start + 1:i - 1]
                        df2.columns = df.iloc[start]
                        df2.reset_index(inplace=True, drop=True)
                        frames.append(df2)
                        index = i - 1
                        break
            if not found:
                index += 1

        # create a list of RADframes from list of read in dataframes
        self.frames = [RADFrame.RADFrame(behaviors=self.behaviors.keys(), frame=i, num_bins=self.num_bins) for i in frames]

        # drop unused columns
        for frame in self.frames:
            frame.drop_columns()
        print len(self.frames)

    def get_keys(self):
        """
        get keys for each object from the plotting sheet in order to match RAD strings
        to data
        """
        plotting_sheet = pd.read_excel(self.filename, sheet_name=2, header=0)
        data_name = plotting_sheet.columns[0].split('_')
        data_name = '_'.join([data_name[1], data_name[3]])
        plotting_sheet = plotting_sheet[18:].reset_index(drop=True)
        plotting_sheet.columns = plotting_sheet.iloc[0]
        frames = 'NA'
        object_names = []
        for index, row in plotting_sheet.iterrows():

            if plotting_sheet['ToUse'][index] == '*':
                dn = data_name
                if not np.isnan(plotting_sheet['Start frame'][index]):
                    frames = str(plotting_sheet['Start frame'][index]) + '_' + str(plotting_sheet['End frame'][index])
                else:
                    frames = str(plotting_sheet['group'][index]) + '_' + str(plotting_sheet['colors 1'][index])
                object_names = [row[1], row[2]] if not np.isnan(row[2]) else [row[1]]
                if isinstance(plotting_sheet['From another project'][index], str):
                    dn = plotting_sheet['From another project'][index].split('_')
                    dn = '_'.join([dn[1], dn[3]])
                self.keys.append([str(o) + '_' + dn + '_' + frames for o in object_names])
        print self.keys

    def __str__(self):
        return str([str(i) for i in self.frames])

    def create_behavior_bins(self):

        for key, behavior in zip(self.behaviors.keys(), self.behaviors.values()):
            for frame in self.frames:
                behavior.add_vals(frame.get_behavior_vals(key))
            behavior.remove_zscore(zscore=2)
            behavior.create_atom_bins()
            print key
            print behavior.vals

    def create_atom_codes(self, num=3):
        dna = ["A", "G", "C", "T"]
        first = 'A'*num
        atoms = [first]

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
        complements = [self.find_complement(i) for i in atoms]
        self.atomList = [(atoms[i], complements[i]) for i in range(0, self.num_bins*len(self.behaviors.keys())+1)]

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

    def avg_over_time_step(self, step):
        for frame in self.frames:
            frame.average_df_over_time_step(step)

    def assign_atom_codes(self):

        for behavior in self.behaviors.values():
            behavior.set_atom_codes(self.atomList)

    def print_multi(self):
        for frame in self.frames:
            print(frame.df.to_string())

    def convert_frames_to_rad(self):
        strings = []
        for frame in self.frames:
            frame.create_rad_matrix(self.behaviors)
            frame.add_start_stop(self.atomList[0][0])
            strings.append(frame.convert_rad_matrix_to_string(frame.rad_df))
        strings_df = pd.DataFrame(columns=['keys', 'Strings'])
        print(len(strings))
        print(self.keys)
        print(len(self.keys))
        strings_df['keys'] = self.keys
        strings_df['Strings'] = strings
        strings_df.to_csv('RAD_STRINGS.csv')