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
                    if score == 1.:
                        print score,snp,motif_id,sta,end,strand
                        print motif_dict[motif_id]
                        raise EOFError
                    all_res.append(score)
        return all_res



if __name__ == "__main__":
    #gene_file_path = os.path.join('/Users/luyi/data',"Norm_raw_gene.txt")
    #id_to_gene, gene_to_id, _ = read_gene_mapping(gene_file_path)
    #for chrm in range(1,23):
    #    TSS_in_chr = get_gene_TSS_seg(gene_file_path,chrm)
    test_data = [(100000,500000),(1000000,1100000),(5000000,6000000),(6000000,7000000)]
    database_path = '/Users/luyi/data/Database/SNPdatabase.db'
    the_SNP = Db_SNP(database_path)
    print the_SNP.get_pos_by_seg(test_data[0],'chr1')
    ##############
    motif_dict = {}
    f_path = '/Users/luyi/data/motif/motifs.txt'
    dist = []
    for name,PWM in read_motif_PWM(f_path):
        motif_dict[name] = PWM
    ssss = the_SNP.find_motif_snp(motif_dict,'/Users/luyi/data/motif/matches.txt','aaa')