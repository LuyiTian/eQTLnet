def read_narrowpeak(file_path):
    '''
    return:
        {id:(chromosome, start, end, signalValue)}
    '''
    res_dict = {}
    index = 0
    for line in open(file_path,'r'):
        items = line.strip().split('\t')
        try:
            res_dict[index] = (int(items[0][3:]), int(items[1]), int(items[2]), float(items[6]))
        except ValueError:
            continue
        index+=1
    return res_dict