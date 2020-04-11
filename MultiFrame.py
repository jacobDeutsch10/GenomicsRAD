import RADFrame
import Atom
import Behavior
import pandas as pd




class MultiFrame:
    """
    a MultiFrame holds a list of RADFrames for all animals in a datset
    """

    def __init__(self, frames=None, behaviors=None):
        """

        :param frames: a list of pandas frames representing the behavior of each animal
        :param behaviors: a list of variables to be used for Genomic conversion
        """

        default = ['vv', 'thetaToR', 'zt', 'curvature']
        self.frames = frames if frames is not None else []
        self.filename = None
        if behaviors is None:
            self.behaviors = {key: Behavior.Behavior(name=key) for key in default}
        else:
            self.behaviors = {key: Behavior.Behavior(name=key) for key in behaviors}

        self.plotting_sheet = None

    def read_from_xl_FULL(self, filename):
        """
        Reads a MultiFrame in from an excel file with data from many animals
        :param filename: the name of the excel file
        """
        self.filename = filename
        # readfile and drop empty columns
        df = pd.read_excel(filename, index=False, sheet_name=0)
        df = df.dropna(how='all')
        index = 0
        frames = []

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
        self.frames = [RADFrame.RADFrame(behaviors=self.behaviors, frame=i) for i in frames]







"""
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
                        complement = RADFrame.find_complement(atom)
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
"""