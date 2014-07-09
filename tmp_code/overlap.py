
def cal_overlap(the_list,ref_list):
    '''
    the_list:
        [(lower_bound,upper_bound,id),..] or
        [(a_float,id),...]
        sorted by lambda x:x[0]
    ref_list [(lower_bound,upper_bound,id),..]
    return :
        {ref_id:[matched_ids],...}
    '''
    if type(the_list) != list:
        return _iter_point_overlap(the_list,ref_list)
    elif len(the_list[0]) == 3:
        return _seg_overlap(the_list,ref_list)
    elif len(the_list[0]) == 2:
        return _point_overlap(the_list,ref_list)
def _seg_overlap(the_list,ref_list):
    STATE = 0
    the_i = ref_i = 0
    res_dict = {}
    while True:
        if ref_i >=len(ref_list):
            print 'all ref used. lefted:',len(the_list)-the_i
            break
        if the_i >=len(the_list):
            #print 'all parsed'
            break
        if STATE == 1:
            if tmp_i >= len(ref_list):
                the_i+=1
                STATE = 0
                continue
            elif the_list[the_i][1] >= ref_list[tmp_i][0] and the_list[the_i][0] <= ref_list[tmp_i][1]:
                res_dict.setdefault(ref_list[tmp_i][2],[]).append(the_list[the_i][2])
                tmp_i+=1
            elif the_list[the_i][1] < ref_list[tmp_i][0]:
                the_i+=1
                STATE = 0
                continue
            elif the_list[the_i][0] > ref_list[tmp_i][1]:
                tmp_i+=1
        elif STATE == 0:
            if the_list[the_i][1] >= ref_list[ref_i][0] and the_list[the_i][0] <= ref_list[ref_i][1]:
                STATE = 1
                tmp_i = ref_i
                res_dict.setdefault(ref_list[ref_i][2],[]).append(the_list[the_i][2])
                tmp_i +=1
            elif the_list[the_i][1] < ref_list[ref_i][0]: the_i+=1
            elif the_list[the_i][0] > ref_list[ref_i][1]: ref_i+=1
    return res_dict
def _point_overlap(the_list,ref_list):
    STATE = 0
    the_i, ref_i = 0
    res_dict = {}
    while True:
        if ref_i >=len(ref_list):
            #print 'all ref used. lefted:',len(the_list)-the_i
            break
        if the_i >=len(the_list):
            #print 'all parsed'
            break
        if STATE == 1:
            if tmp_i >= len(ref_list):
                the_i+=1
                STATE = 0
                continue
            elif the_list[the_i][0] >= ref_list[tmp_i][0] and the_list[the_i][0] <= ref_list[tmp_i][1]:
                res_dict.setdefault(ref_list[tmp_i][2],[]).append(the_list[the_i][1])
                tmp_i+=1
            elif the_list[the_i][0] < ref_list[tmp_i][0]:
                the_i+=1
                STATE = 0
                continue
            elif the_list[the_i][0] > ref_list[tmp_i][1]:
                tmp_i+=1
        elif STATE == 0:
            if the_list[the_i][0] >= ref_list[ref_i][0] and the_list[the_i][0] <= ref_list[ref_i][1]:
                STATE = 1
                tmp_i = ref_i
                res_dict.setdefault(ref_list[ref_i][2],[]).append(the_list[the_i][1])
                tmp_i +=1
            elif the_list[the_i][0] < ref_list[ref_i][0]: the_i+=1
            elif the_list[the_i][0] > ref_list[ref_i][1]: ref_i+=1
    return res_dict
def _iter_point_overlap(iter_list,ref_list):
    res_dict = {}
    mapping_dict = {}
    ref_i = 0
    ahaha = 0
    for items in iter_list:
        if ref_i >=len(ref_list):
            print 'all ref used. '
            break
        while True:
            pos = float(items[0])
            if pos >= ref_list[ref_i][0] and pos <= ref_list[ref_i][1]:
                res_dict[pos] = [float(ff) for ff in items[1:]]
                mapping_dict.setdefault(ref_list[ref_i][2],[]).append(pos)
                tmp_i = ref_i
                tmp_i+=1
                if tmp_i >= len(ref_list):break
                while pos >= ref_list[tmp_i][0]:
                    if pos <= ref_list[tmp_i][1]:
                        mapping_dict.setdefault(ref_list[tmp_i][2],[]).append(pos)
                    tmp_i+=1
                    if tmp_i >= len(ref_list):break
                break
            elif pos < ref_list[ref_i][0]: 
                break
            elif pos > ref_list[ref_i][1]:
                ref_i+=1      
                if ref_i >=len(ref_list):
                    break 
        ahaha += 1   
    return res_dict,mapping_dict
if __name__ == "__main__":
    pass