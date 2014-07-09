##import sys
import os
import re
##sys.path.append("..")
from parameters import TF_peak_dir,data_dir
cell_type = "Gm12878"
file_list = os.listdir(TF_peak_dir)


out_dir = os.path.join(TF_peak_dir, cell_type)
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
for i in range(1,len(file_list)):
    sourceFile = os.path.join(TF_peak_dir, file_list[i])
    targetFile = os.path.join(out_dir,file_list[i].split("_VS_")[0])
    if os.path.isfile(sourceFile):
        if cell_type in file_list[i]:
            if not os.path.exists(targetFile) or\
            (os.path.exists(targetFile) and (os.path.getsize(targetFile) != os.path.getsize(sourceFile))): 
                open(targetFile, "wb").write(open(sourceFile, "rb").read())


def get_antibody_dict(des_file_name):
    res_dict = {}
    ANTIBODY_STATE = 0
    term = None
    for line in open(des_file_name,'r'):
        if line.strip() == "#### Antibodies":
            ANTIBODY_STATE = 1
            continue
        if ANTIBODY_STATE ==1:
            if line[:5] == "term ":
                term = line.strip().split(' ')[1].upper()
            elif line[:7] == "target ":
                if term:
                    res_dict[term] = line.strip().split(' ')[1].upper()
                    term = None
                else:
                    print "no match term for target",line.strip().split(' ')[1].upper()
                    raise Exception
    return res_dict


# vc.ra is from http://genome.ucsc.edu/ENCODE/downloads.html
convert_dict = get_antibody_dict(os.path.join(data_dir, "cv.ra"))

selected_list = os.listdir(out_dir)
pattern = re.compile("Gm12878(\S+?)(Pcr\dx|Std|Iggmus|Iggrab){0,1}AlnRep")
for i in range(1,len(selected_list)):
    new_name = pattern.findall(selected_list[i])[0][0].upper()
    if convert_dict.has_key(new_name):new_name = convert_dict[new_name]
    print new_name
    sourceFile = os.path.join(out_dir, selected_list[i])
    targetFile = os.path.join(out_dir, str(i)+"__"+new_name)
    os.rename(sourceFile,targetFile)
