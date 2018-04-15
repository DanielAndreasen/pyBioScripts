import matplotlib.pyplot as plt


class DNA:
    DNA_DICT = {'N': 0, '-': 0,
                'A': 1, 'C': 2, 'G': 3, 'T': 4, 'U': 4,
                'R': 6, 'Y': 7, 'S': 8, 'W': 9, 'K': 10, 'M': 11}

    def __init__(self, name, sequence):
        assert isinstance(name, str)
        assert isinstance(sequence, str)
        self.name = name
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __sub__(self, other):
        if not isinstance(self, other.__class__):
            raise TypeError("{} cannot be used with DNA type".format(other.__class__))
        if len(self) != len(other):
            raise ValueError("The two sequences must have the same length")
        res = []
        for s1, s2 in zip(self.sequence, other.sequence):
            # Trivial cases
            if s1 == s2:
                res.append(0)
            # TODO: Check with Bea
            elif s1 == 'N' and s2 in '-N':
                res.append(0)
            elif s1 == '-' and s2 in '-N':
                res.append(0)
            # R case
            elif s1 == 'R' and (s2 in 'AG'):
                res.append(0.5)
            elif s2 == 'R' and (s1 in 'AG'):
                res.append(0.5)
            # Y case
            elif s1 == 'Y' and (s2 in 'CT'):
                res.append(0.5)
            elif s2 == 'Y' and (s1 in 'CT'):
                res.append(0.5)
            # S case
            elif s1 == 'S' and (s2 in 'GC'):
                res.append(0.5)
            elif s2 == 'S' and (s1 in 'GC'):
                res.append(0.5)
            # W case
            elif s1 == 'W' and (s2 in 'AT'):
                res.append(0.5)
            elif s2 == 'W' and (s1 in 'AT'):
                res.append(0.5)
            # K case
            elif s1 == 'K' and (s2 in 'GT'):
                res.append(0.5)
            elif s2 == 'K' and (s1 in 'GT'):
                res.append(0.5)
            # M case
            elif s1 == 'M' and (s2 in 'AC'):
                res.append(0.5)
            elif s2 == 'M' and (s1 in 'AC'):
                res.append(0.5)
            # No match at all
            else:
                res.append(1)
        return res

    def __repr__(self):
        return '{}: {}'.format(self.name, self.sequence)

    def equalness(self, other, verbose=False):
        res = self - other
        N = len(res)
        different = len(list(filter(None, res)))  # Remove zeros
        p = round((N - different)/N*100, 2)
        if verbose:
            print('Equalness: {}%'.format(p))
        return p

    def _colors(self):
        res = ['C0' if v in '-N' else 'C1' if v in 'ACGTU' else 'C2' for v in self.sequence]
        res = []
        for v in self.sequence:
            if v in '-N':
                res.append('k')
            elif v == 'A':
                res.append('C0')
            elif v == 'C':
                res.append('C1')
            elif v == 'G':
                res.append('C2')
            elif v in 'TU':
                res.append('C4')
            else:
                res.append('C3')
        return res

    def plot(self, other):
        p = self.equalness(other)
        d = self - other
        x1, y1, c1 = list(range(len(self.sequence))),  [2]*len(self),  self._colors()
        x2, y2, c2 = x1, [1]*len(other), other._colors()
        x3, y3, c3 = x1, [-1]*len(d), ['w' if v==0 else 'C1' if v==0.5 else 'C3' for v in d]
        s = [0 if v==0 else 10 if v==0.5 else 20 for v in d]
        dx = x1[1]-x1[0]
        plt.figure(figsize=(8, 3.5))
        plt.title("Match of {} and {}\nEqualness: {}%".format(self.name, other.name, p))
        plt.yticks([-1, 1, 2], ['Diff', self.name, other.name])
        plt.xticks([], [])
        plt.scatter(x1, y1, c=c1)
        plt.scatter(x2, y2, c=c2)
        plt.hlines(0, x1[0]-dx, x1[-1]+dx, linestyles='dashed')
        plt.scatter(x3, y3, c=c3, s=s)
        plt.ylim(-10, 10)
        plt.show()


if __name__ == '__main__':
    d1 = DNA('DNA1', 'AGTCTGNNN-NTGAAA-A')
    d2 = DNA('DNA2', 'AGTTTGNNN-NTGAAMMA')

    d1.plot(d2)
