#motif.py
# * 1 indexed
# * end inclusive
#
import numpy as np

def read_motif_PWM(file_path):
    motif_text = ''
    PWM = []
    for line in open(file_path):
        if line[0] == '>':
            if PWM:
                PWM = np.array(PWM)
                yield motif_text,PWM
                PWM = None
                motif_text = ''
            else:
                motif_text = line.strip()[1:]
        else:
            PWM.append([float(it) for it in line.strip().split()[1:]])
    PWM = np.array(PWM)
    yield motif_text,PWM


if __name__ == '__main__':
    f_path = '/Users/luyi/data/motif/motifs.txt'
    for a,b in read_motif_PWM(f_path):
        print a
        print b
        break
