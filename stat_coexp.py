from parameters import processed_data_dir
def get_in_out(file_name,cutoff=0.29):
    res_dict={}
    for line in open(file_name,'r'):
        n1,n2,val = line.strip().split('\t')
        if abs(float(val))>cutoff:
            res_dict.setdefault(int(n1),[]).append(float(val))
            res_dict.setdefault(int(n2),[]).append(float(val))
    return res_dict
if __name__ =='__main__':
    file_name = processed_data_dir+"out_coexpression"
    res = get_in_out(file_name)
    tmp=[]
    for key,val in res.items():
        tmp.append((key,len(val)))
    tmp = sorted(tmp, key = lambda x:x[1],reverse = True)
    for i in range(20):print tmp[i]
