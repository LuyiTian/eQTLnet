
def read_previous_eqtl(file_name):
    '''
    return dictionary:
    key: gene name i.e: ENSG00000137411.10
    value: position of the variant
    '''
    result_dict={}
    for line in open(file_name, 'r'):
        items = line.strip().split('\t')
        result_dict[items[2]] = int(items[6])#assume cis-qtl
    return result_dict
def read_gene_mapping(file_name):
    res_dict={}
    for line in open(file_name,'r'):
        items = line.strip().split('\t')
        res_dict[int(items[0])] = (items[1],items[3])
    return res_dict    
def read_undirected_n(file_name,gene_mapping):
    res_dict={}
    for line in open(file_name,'r'):
        n1,n2 = [int(i) for i in line.strip().split('\t')[:2]]
        res_dict.setdefault(gene_mapping[n1][0],[]).append(gene_mapping[n2][0])
        res_dict.setdefault(gene_mapping[n2][0],[]).append(gene_mapping[n1][0])
    return res_dict
if __name__ == '__main__':
    file_name = "G:\\res_files\\EUR373.gene.cis.FDR5.best.rs137.txt"
    res = read_previous_eqtl(file_name)