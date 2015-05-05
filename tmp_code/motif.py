#motif.py
# * 1 indexed
# * end inclusive
#
import numpy as np

NT_ORDER = {'+':{'A':0, 'C':1, 'G':2, 'T':3},'-':{'A':3, 'C':2, 'G':1, 'T':0}}


def read_motif_PWM(file_path):
    '''
    PWM format:
                A  C  G  T
    position1   *  *  *  *
    position2   *  *  *  *
    ...
    @param: file path contains all motif PWM in fasta format 
    '''
    motif_text = ''
    PWM = []
    for line in open(file_path):
        if line[0] == '>':
            if PWM:
                PWM = np.array(PWM)
                yield motif_text,PWM
                PWM = []
                motif_text = line.strip().split()[0][1:]
            else:
                motif_text = line.strip().split()[0][1:]
        else:
            PWM.append([float(it) for it in line.strip().split()[1:]])
    PWM = np.array(PWM)
    yield motif_text,PWM


if __name__ == '__main__':
    f_path = '/Users/luyi/data/motif/motifs.txt'
    for a,b in read_motif_PWM(f_path):
        print a
        print b
        print b[1, 3]
        break
