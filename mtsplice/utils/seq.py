import numpy as np

def randomseq(length=100):
    alphabet = "ATGC"
    indx = np.random.randint(0,3,length)
    seq = [alphabet[s] for s in indx]
    seq = ''.join(seq)
    return seq

bases = ['A', 'C', 'G', 'T']
def onehot(seq):
    X = np.zeros((len(seq), len(bases)))
    for i, char in enumerate(seq):
        if char == "N":
            pass
        else:
            X[i, bases.index(char.upper())] = 1
    return X
