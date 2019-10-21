def matrixParser(filename):
    with open(filename, encoding='utf-8') as f:
        matrix = SubstitutionMatrix()
        for line in f:
            if line[0] not in ["#", " ", "\n"]:
                matrix.Append(line[0])
                matrix.addRow(line[1:].split())
    return matrix


def sequencesParser(filename):
    with open(filename, encoding='utf-8') as f:
        mylist = []
        sequence = ''
        for line in f:
            if line[0] not in [">", "\n"]:
                sequence += line.strip("\n")
            else:
                if len(sequence):
                    mylist.append(ADTsequence(sequence))
                    sequence = ''
    return mylist


class ADTsequence(object):

    def __init__(self, sequence):
        self.sequence = sequence

    def __getitem__(self, key):
        return self.sequence[key]

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        output = ""
        for i in range(len(self) // 60 + 1):
            output += self.sequence[i * 60:(i + 1) * 60]
            output += "\n"
        return output[:-1]


class ADTmatrix:

    def __init__(self, n=None, m=None):

        if n and m:
            self.matrix = [[0 for i in range(m)] for j in range(n)]
        else:
            self.matrix = []

    def getItem(self, i, j):
        return self.matrix[i][j]

    def setItem(self, i, j, value):
        self.matrix[i][j] = value

    def addRow(self, row):
        self.matrix.append(row)

    def __str__(self):
        output = ""
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[0])):
                output += str(self.matrix[i][j])
                output += " "
            output += "\n"
        return output[:-1]


class SubstitutionMatrix(ADTmatrix):

    def __init__(self, n=None, m=None):
        super().__init__(n, m)
        self.myList = []

    def __getitem__(self, key):
        return self.myList.index(key)

    def Append(self, value):
        self.myList.append(value)


class ScoringMatrix(ADTmatrix):

    def __init__(self, character, mode, sequence1, sequence2, I=None, E=None):
        super().__init__(len(sequence1) + 1, len(sequence2) + 1)
        self.sequence1 = sequence1
        self.sequence2 = sequence2
        self.n = len(sequence1) + 1
        self.m = len(sequence2) + 1

        if character == "V":
            self.initV()
        elif character == "W":
            self.initW()
        else:
            self.initS(mode, I, E)

    def initS(self, mode, I, E):

        for i in range(1, self.m):
            gap = - I - (i - 1) * E
            if mode == "LOCAL":
                self.setItem(0, i, max(0, gap))
            else:
                self.setItem(0, i, gap)

        for j in range(1, self.n):
            gap = - I - (j - 1) * E
            if mode == "LOCAL":
                self.setItem(j, 0, max(0, gap))
            else:
                self.setItem(j, 0, gap)

    def initV(self):
        for i in range(self.m):
            self.setItem(0, i, - float("inf"))

        for j in range(1, self.n):
            self.setItem(j, 0, 0)

    def initW(self):
        for i in range(1, self.m):
            self.setItem(0, i, 0)

        for j in range(self.n):
            self.setItem(j, 0, - float("inf"))


class NeedlemanWunsch(object):
    def __init__(self, S, V, W, MATSUB, k, I, E):

        self.S = S
        self.V = V
        self.W = W
        self.MATSUB = MATSUB
        self.k = k
        self.I = I
        self.E = E

        for i in range(1, self.S.n):
            for j in range(1, self.S.m):
                self.fill(i, j)
        self.findWay()

    def findWay(self, i=None, j=None, way=[], ALIGN=[], SOLUTIONS=[]):

        if i is None and j is None:
            i = self.S.n - 1
            j = self.S.m - 1

        value = self.S.getItem(i, j)
        way.append(value)

        if i == 0 and j == 0:
            if len(SOLUTIONS) < self.k:
                SOLUTIONS.append((way[::-1], ALIGN[::-1]))

        else:

            if value == self.S.getItem(i - 1, j - 1) + self.t(i, j):
                ALIGN.append((self.S.sequence1[i - 1], self.S.sequence2[j - 1]))
                self.findWay(i - 1, j - 1, way, ALIGN)
                ALIGN.pop()

            if value == self.V.getItem(i, j) or j == 0:
                ALIGN.append(("-", self.S.sequence2[j - 1]))
                self.findWay(i - 1, j, way, ALIGN)
                ALIGN.pop()

            if value == self.W.getItem(i, j) or i == 0:
                ALIGN.append((self.S.sequence1[i - 1], "-"))
                self.findWay(i, j - 1, way, ALIGN)
                ALIGN.pop()

        way.pop()
        return SOLUTIONS

    def t(self, i, j):
        x = self.MATSUB[self.S.sequence2[j - 1]]
        y = self.MATSUB[self.S.sequence1[i - 1]]
        return int(self.MATSUB.getItem(x, y))

    def fill(self, i, j):
        Vij = max(self.S.getItem(i - 1, j) - self.I, self.V.getItem(i - 1, j) - self.E)
        Wij = max(self.S.getItem(i, j - 1) - self.I, self.W.getItem(i, j - 1) - self.E)
        Sij = max(self.S.getItem(i - 1, j - 1) + self.t(i, j), Vij, Wij)

        self.V.setItem(i, j, Vij)
        self.W.setItem(i, j, Wij)
        self.S.setItem(i, j, Sij)


class SmithWaterman(object):
    def __init__(self, S, V, W, MATSUB, l, I, E):

        self.S = S
        self.V = V
        self.W = W
        self.MATSUB = MATSUB
        self.I = I
        self.E = E

        for i in range(1, self.S.n):
            for j in range(1, self.S.m):
                self.fill(i, j)

        i = 0
        while i < l and sum([sum(row) for row in self.S.matrix]):
            indexList = []
            self.topDown(indexList)
            self.fillZero(indexList)
            x, y = indexList[-1]
            self.computeAgain(x, y)
            i += 1

    def t(self, i, j):
        x = self.MATSUB[self.S.sequence2[j - 1]]
        y = self.MATSUB[self.S.sequence1[i - 1]]
        return int(self.MATSUB.getItem(x, y))

    def fill(self, i, j):

        Vij = max(self.S.getItem(i - 1, j) - self.I, self.V.getItem(i - 1, j) - self.E)
        Wij = max(self.S.getItem(i, j - 1) - self.I, self.W.getItem(i, j - 1) - self.E)
        Sij = max(self.S.getItem(i - 1, j - 1) + self.t(i, j), Vij, Wij)

        self.V.setItem(i, j, max(0, Vij))
        self.W.setItem(i, j, max(0, Wij))
        self.S.setItem(i, j, max(0, Sij))

    def maxValue(self):

        x, y = 1, 1
        maximum = self.S.getItem(x, y)
        for i in range(1, self.S.n):
            for j in range(1, self.S.m):
                if maximum < self.S.getItem(i, j):
                    maximum = self.S.getItem(i, j)
                    x, y = i, j
        return x, y

    def topDown(self, indexList):

        i, j = self.maxValue()
        indexList.append((i, j))
        score = self.S.getItem(i, j)
        seq = ""
        seq2 = ""
        seq3 = ""
        pos = (i,j)
        while self.S.getItem(i, j):

            value = self.S.getItem(i, j)
            pos2 = (i,j)

            if value == self.S.getItem(i - 1, j - 1) + self.t(i, j):
                seq += self.S.sequence1[i - 1]
                seq2 += self.S.sequence2[j - 1]
                if self.S.sequence1[i - 1] == self.S.sequence2[j - 1]:
                    seq3 += "|"
                else:
                    seq3 += ":" if self.t(i, j) >= 0 else "."
                i, j = i - 1, j - 1

            elif value == self.V.getItem(i, j):
                seq += self.S.sequence1[i - 1]
                seq2 += "-"
                seq3 += " "
                i = i - 1

            elif value == self.W.getItem(i, j):
                seq += "-"
                seq2 += self.S.sequence2[j - 1]
                seq3 += " "
                j = j - 1

            indexList.append((i, j))

        show(score, seq[::-1], seq2[::-1], seq3[::-1],pos,pos2)
        indexList.pop()

    def fillZero(self, indexList):

        for elem in indexList:
            self.S.setItem(elem[0], elem[1], 0)
            self.V.setItem(elem[0], elem[1], 0)
            self.W.setItem(elem[0], elem[1], 0)

    def computeAgain(self, x, y):
        i, j = x, y
        while i < self.S.n:
            while j < self.S.m:
                if self.S.getItem(i,j) > 0 or self.V.getItem(i,j) > 0 or self.W.getItem(i,j) > 0:
                    self.fill(i, j)
                j += 1
            i += 1
            j = y

def show(score, sequence, sequence2, sequence3,pos,pos2):
    for i in range(len(sequence) // 60 + 1):
        print(sequence[i * 60:(i + 1) * 60] \
            + " - from {0} to {1}".format(pos2[0],pos[0]) \
            + "\n" + sequence3[i * 60:(i + 1) * 60] \
            + "\n" + sequence2[i * 60:(i + 1) * 60] \
            + " - from {0} to {1}".format(pos2[1],pos[1]) +"\n")

    print("Score :", score)
    print("Identity percentage :", round(sequence3.count("|") / len(sequence3) * 100,2), "%")
    print("Similarity rate:", round((sequence3.count("|") + sequence3.count(":")) / len(sequence3) * 100,2), "%")
    print("-" * 60)

MATSUB = matrixParser("pam120.txt")
sequences = sequencesParser("WW-sequence.fasta")
for i in range(len(sequences)):
    for j in range(i+1,len(sequences)) :
            print("Alignment between sequence {0} and {1}".format(i+1,j+1))
            S = ScoringMatrix("S","GLOBAL",sequences[i],sequences[j],6,3)
            V = ScoringMatrix("V","GLOBAL",sequences[i],sequences[j])
            W = ScoringMatrix("W","GLOBAL",sequences[i],sequences[j])
            NeedlemanWunsch(S,V,W,MATSUB,5,6,3)

sequences = sequencesParser("protein-sequences.fasta")
for i in range(len(sequences)):
    for j in range(i+1,len(sequences)) :
            print("Alignment between sequence {0} and {1}".format(i+1,j+1))
            S = ScoringMatrix("S","LOCAL",sequences[i],sequences[j],8,1)
            V = ScoringMatrix("V","LOCAL",sequences[i],sequences[j])
            W = ScoringMatrix("W","LOCAL",sequences[i],sequences[j])
            SmithWaterman(S, V, W, MATSUB, 3, 8, 1)
