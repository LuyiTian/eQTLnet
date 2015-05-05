import sqlite3
import cPickle
#from parameters import TF_peak_dir,data_dir,chrom_list
from read_processed_data import read_gene_mapping,get_gene_TSS_seg, read_expression_value
from overlap import cal_overlap
import os
from motif import NT_ORDER, read_motif_PWM
class Db_SNP():
    def __init__(self, database_path):
        self.con = sqlite3.connect(database_path)
    def get_all_pos_by_chrm(self, chrm):
        cur=self.con.cursor()
        cur.execute("SELECT position FROM {chromosome}".format(chromosome = chrm))
        res = [i[0] for i in cur.fetchall()]
        return res
    def get_pos_by_seg(self,seg,chrm):
        '''
        @param:TODO
        '''
        cur = self.con.cursor()
        cur.execute("SELECT * FROM {}_base WHERE position >= ? and position <= ?".format(chrm),seg)
        return cur.fetchall()
    def find_motif_snp(self,motif_dict, match_file, out_file,thr = 0.8):
        '''
        @param:TODO
        '''
        all_res = []
        #out_f = open(out_file,'w')
        site_dist_dict = {}
        for key in motif_dict:
            site_dist_dict[key] = [0]*motif_dict[key].shape[0]
        for line in open(match_file):
            motif_id, chrm, sta,end, strand = line.strip().split()
            if chrm != 'chr1':
                break
            snp_in_seg = self.get_pos_by_seg((sta,end),chrm)
            sta,end = int(sta),int(end)
            if len(snp_in_seg) == 0:
                continue
            else:
                for snp in snp_in_seg:
                    if len(snp[1])>1 or len(snp[2])>1:
                        continue
                    score = motif_dict[motif_id][int(snp[0])-sta,\
                    NT_ORDER[strand][snp[1]]]-motif_dict[motif_id][int(snp[0])-sta,NT_ORDER[strand][snp[2]]]
                    if score >= -1.9:
                        site_dist_dict[motif_id][int(snp[0])-sta]+=1
                        #print score,snp,motif_id,sta,end,strand
                        #all_res.append(motif_id.split('_')[0])
        return site_dist_dict


if __name__ == "__main__":
    #gene_file_path = os.path.join('/Users/luyi/data',"Norm_raw_gene.txt")
    #id_to_gene, gene_to_id, _ = read_gene_mapping(gene_file_path)
    #for chrm in range(1,23):
    #    TSS_in_chr = get_gene_TSS_seg(gene_file_path,chrm)
    ##############
    import numpy as np
    motif_dict = {}
    f_path = '/Users/luyi/data/motif/motifs.txt'
    dist = []
    for name, PWM in read_motif_PWM(f_path):
        motif_dict[name] = PWM
    the_SNP = Db_SNP('/Users/luyi/data/Database/SNPdatabase.db')
    ssss = the_SNP.find_motif_snp(motif_dict,'/Users/luyi/data/motif/matches.txt', 'aaa')
    site_dist_tuple = [(key,max(ssss[key])) for key in ssss]
    site_dist_tuple.sort(key=lambda x: x[1], reverse=True)
    for it in site_dist_tuple[:20]:
        print it[0]
        tmp = np.array(ssss[it[0]])
        tmp = np.atleast_2d(tmp).T
        print np.hstack((tmp,motif_dict[it[0]]))
        print '###############'
    print '~'*20
    for it in site_dist_tuple[-10+len(site_dist_tuple)/2:10+len(site_dist_tuple)/2]:
        print it[0]
        tmp = np.array(ssss[it[0]])
        tmp = np.atleast_2d(tmp).T
        print np.hstack((tmp,motif_dict[it[0]]))
        print '###############'
    #aaaa = [(the_id,ssss.count(the_id)) for the_id in list(set(ssss))]
    #aaaa.sort(key=lambda x: x[1])
    #print aaaa[-20:]
