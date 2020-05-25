
class Atom:

    def __init__(self, code=None):
        self.code = code if code is not None else None
        self.complement = None
        if code is not None:
            self.find_complement()
        self.lowerBound = 0
        self.upperBound = 100

    def set_code(self, code):
        self.code = code
        self.find_complement()

    def get_code(self):
        return self.code

    def get_complement(self):
        return self.complement

    def find_complement(self):
        self.complement = ''
        for i in self.code:
            if i == 'A':
                self.complement += 'T'
            elif i == 'T':
                self.complement += 'A'
            elif i == 'C':
                self.complement += 'G'
            elif i == 'G':
                self.complement += 'C'

    def set_lower_bound(self, lower):
        self.lowerBound = lower

    def get_lower_bound(self):
        return self.lowerBound

    def set_upper_bound(self, upper):
        self.upperBound = upper

    def get_upper_bound(self):
        return self.upperBound

    def is_in_atom(self, num):
        return self.lowerBound < num <= self.upperBound
