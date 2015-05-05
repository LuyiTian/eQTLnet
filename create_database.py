'''
create sqlite database from tab delimited files.
'''
import sqlite3
import vcf

def add_chr_snp(database_path, file_path, sample_list):
    '''
    add sno info by chr in databases
    input:
        database_path:
            can be opened by sqlite3.connect(database_path)
        file_path:
            file name include file path that can be opened by 'open(file_path)'
            file format:
                position genotypes1 genotypes2 ...
                ...
        sample_list:
            a list of sample id, the order should be the same in file_path
    '''
    con = sqlite3.connect(database_path)
    cur = con.cursor()
    chrm = file_path.split('_')[-1]
    cur.execute("CREATE TABLE IF NOT EXISTS {chromosome} (position FLOAT PRIMARY KEY,{samples})".format(
        chromosome = 'chr'+chrm, samples = ', '.join([i+' INTEGER' for i in sample_list])))
    con.commit()
    for num,line in enumerate(open(file_path,'r')):
        items = line.strip().split(' ')
        cur.execute('INSERT INTO {chromosome} VALUES ({position}, {genotypes})'.format(chromosome = 'chr'+chrm, position = items[0],genotypes = ', '.join(items[1:])))
        if num>0 and num%10000==0:
            con.commit()
            print 'INSERT {} lines'.format(num)
    cur.close()
    con.close()

def add_chr_snp(database_path, file_path, chrm):
    '''
    add snp info by chr in databases
    input:
        database_path:
            can be opened by sqlite3.connect(database_path)
        file_path:
            file name include file path that can be opened by 'open(file_path)'
            file format:
                position genotypes1 genotypes2 ...
                ...
        sample_list:
            a list of sample id, the order should be the same in file_path
    '''
    con = sqlite3.connect(database_path)
    cur = con.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS {chromosome}_base (position FLOAT PRIMARY KEY, REF CHAR(20), ALT CHAR(20))".format(
        chromosome = chrm))
    con.commit()
    vcf_reader = vcf.Reader(open(file_path, "rb"))
    for num,record in enumerate(vcf_reader):
        if len(str(record.REF))<20 and len(str(record.ALT[0]))<20:
            cur.execute("INSERT INTO {chromosome}_base VALUES ({position}, '{ref}', '{alt}')"\
            .format(chromosome = chrm, position = record.POS, ref = str(record.REF), alt = str(record.ALT[0])))
        if num>0 and num%20000==0:
            con.commit()
            print 'INSERT {} lines'.format(num)
    cur.close()
    con.close()

if __name__== "__main__":
    snp_path = '/Users/luyi/data/snp_processed/'
    sample_file = '/Users/luyi/data/Norm_raw_sample_id.txt'
    database_path = '/Users/luyi/data/Database/SNPdatabase.db'
    add_chr_snp(database_path, '/Users/luyi/data/snp/GEUVADIS.chr1.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz', 'chr1')
    '''
    sample_list = []
    for line in open(sample_file,'r'):
        sample_list.append(line.strip())
    for ch in range(1,23):
        print 'processing chromosome: {}'.format(ch)
        snp_file = snp_path+'tmp_snp_'+str(ch)
        add_chr_snp(database_path, snp_file, sample_list)
    '''


