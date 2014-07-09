def read_gff_by_chr(chr):
    f = open("AnnotatedFeatures.gff",'r')
    res_dict = {}
    for line in f:
        items = line.strip().split('\t')
        
        if items[0] ==chr:
            res_dict[(int(items[3]),int(items[4]))] = items[8].split(';')[0].split('=')[1]
    return res_dict
def read_eQTL(f_name):
    res = []
    f = open(f_name,'r')
    f.readline()
    for line in f:
        snp = line.strip().split('\t')[0].split('_')[1]
        res.append(float(snp))
    res = list(set(res))
    print len(res),"eQTLs detacted in file",f_name
    res=sorted(res)
    return res
def cal_overlap(snp_list, anno_list, dis_to_TSS):
    STATE = 0
    anno_i = 0
    snp_i = 0
    res_dict={}.fromkeys(snp_list,0)
    
    while True:
        if anno_i >= len(anno_list):
            print 'all annotation used'
            print 'lefted snps:',len(snp_list)-snp_i
            break
        if snp_i >= len(snp_list):
            print 'all snp parsed'
            break
        if snp_i%1000 == 0:print '---',snp_i,snp_list[snp_i]
        #print '~~~~~~',snp_list[snp_i],anno_list[anno_i]
        if STATE == 1:
            if tmp_anno_i >= len(anno_list):
                snp_i+=1
                STATE = 0
                #print 'change state',tmp_anno_i
                continue
            elif snp_list[snp_i]>=anno_list[tmp_anno_i][0] and snp_list[snp_i]<=anno_list[tmp_anno_i][1]:
                res_dict[snp_list[snp_i]]+=1
                dis_to_TSS.append(snp_list[snp_i]-anno_list[tmp_anno_i][2])
                tmp_anno_i+=1
            elif snp_list[snp_i]<anno_list[tmp_anno_i][0]:
                snp_i+=1
                STATE = 0
                continue
            elif snp_list[snp_i]>anno_list[tmp_anno_i][1]:
                tmp_anno_i+=1
            #print 'tmp_state',snp_list[snp_i],anno_list[tmp_anno_i-1],tmp_anno_i
        elif STATE == 0 :
            #print 'main_state',snp_list[snp_i],anno_list[anno_i]
            if snp_list[snp_i]>=anno_list[anno_i][0] and snp_list[snp_i]<=anno_list[anno_i][1]:
                STATE = 1
                tmp_anno_i = anno_i
                res_dict[snp_list[snp_i]]+=1
                dis_to_TSS.append(snp_list[snp_i] - anno_list[anno_i][2])
                tmp_anno_i+=1
                continue
            elif snp_list[snp_i]<anno_list[anno_i][0]:
                snp_i+=1
            elif snp_list[snp_i]>anno_list[anno_i][1]:
                anno_i+=1
    return res_dict
def cal_dis_to_TSS(snp_list,TSS_list,distance = 10000):
    pass
def get_gene_TSS_in_chr(f_name,chr=None,distance = 3000):
    res = []
    for line in open(f_name,'r'):
        items = line.strip().split('\t')
        if items[3] == str(chr):
            res.append((int(items[4])-distance,int(items[4])+distance,int(items[0])))
    return res
if __name__ == "__main__":
    a = read_eQTL('out_3')
    #b = read_gff_by_chr('chr2')
    b = get_gene_TSS_in_chr("Norm_raw_gene.txt",3)
    print len(a),len(b)
    for i in range(10):print b[i]
    dis_to_TSS = []
    c = cal_overlap(a,sorted(b,key = lambda x:x[0]),dis_to_TSS)
    cot = 0
    #for key,val in c.items():
    #    if val>0:cot +=1
    #print '-'*5,len(a),cot
    import pylab as pl
    pl.hist(dis_to_TSS,bins = 60)
    pl.show()