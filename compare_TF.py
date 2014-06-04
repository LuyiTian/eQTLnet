from parameters import processed_data_dir

TF_network_file = open(processed_data_dir+"enets1.Proximal_raw.txt",'r')
id_convert_file = open(processed_data_dir+"mart_export.txt",'r')
eqtl_network_file = open(processed_data_dir+"direct_network",'r')

def read_id_convert(id_convert_file):
    res_dict = {}
    id_convert_file.readline()
    for line in id_convert_file:
        id1,id2 = line.strip().split('\t')
        #print id1,id2
        res_dict[id2] = int(id1[4:])
    print res_dict["CTCF"]
    return res_dict

def read_TF_network(TF_network_file, id_dict):
    TF_network = {}
    for line in TF_network_file:
        fm,_,to = line.strip().split(' ')
        if id_dict.has_key(fm) and id_dict.has_key(to):
            TF_network.setdefault(id_dict[fm],[]).append(id_dict[to])
        else:
            pass
            #print "cannot convert id",fm,to
    return TF_network

def read_eQTL_network(eqtl_network_file, TF_network):
    '''
    if score <0, M->A<-B conditioning on A activities M-B
    if score >0, M->A->B M B is independent conditioning on A 
    '''
    eqtl_network = {}
    for line in eqtl_network_file:
        item = line.strip().split('\t')
        score = float(item[8])-float(item[6])#p_MB - p_partial_coeff
        n1 = int(item[2].split('.')[0][5:])
        n2 = int(item[3].split('.')[0][5:])
        if TF_network.has_key(n1):
            eqtl_network.setdefault(n1, {})
            if eqtl_network[n1].has_key(n2):
                if score*eqtl_network[n1][n2]<0:
                    print "edge direction conflict",n1,n2,score,eqtl_network[n1][n2]
                    if (score>0 and score>-eqtl_network[n1][n2]) or (score<0 and score<-eqtl_network[n1][n2]):
                        eqtl_network[n1][n2] = score#update score
            else:
                eqtl_network[n1][n2]=score
        elif TF_network.has_key(n2):
            score = -score
            eqtl_network.setdefault(n2, {})
            if eqtl_network[n2].has_key(n1):
                if score*eqtl_network[n2][n1]<0:
                    print "edge direction conflict",n2,n1,score,eqtl_network[n2][n1]
                    if (score>0 and score>-eqtl_network[n2][n1]) or (score<0 and score<-eqtl_network[n2][n1]):
                        eqtl_network[n2][n1] = score#update score
            else:
                eqtl_network[n2][n1]=score
    return eqtl_network
def compare_net(eqtl_network, TF_network, threshold = 0.1):
    True_dict = {}
    False_dict = {}
    Unknown_dict= {}
    for TF in TF_network.keys():
        if not eqtl_network.has_key(TF):
            print TF,"not found in eQTL network"
        else:
            for value in TF_network[TF]:
                if eqtl_network[TF].has_key(value):
                    if eqtl_network[TF][value]>threshold: True_dict.setdefault(TF,[]).append(eqtl_network[TF][value])
                    elif eqtl_network[TF][value]<-threshold: False_dict.setdefault(TF,[]).append(eqtl_network[TF][value])
                    else:Unknown_dict.setdefault(TF,[]).append(eqtl_network[TF][value])
    for TF in TF_network.keys():
        if True_dict.has_key(TF):T = len(True_dict[TF])
        else:T = 0
        if False_dict.has_key(TF):F = len(False_dict[TF])
        else:F = 0
        if Unknown_dict.has_key(TF):U = len(Unknown_dict[TF])
        else:U = 0
        print TF,T,F,U
        
if __name__ == "__main__":
    id_dict = read_id_convert(id_convert_file)
    TF_network = read_TF_network(TF_network_file, id_dict)
    eqtl_network = read_eQTL_network(eqtl_network_file, TF_network)
    compare_net(eqtl_network, TF_network)
            