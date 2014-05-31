import multiprocessing
from read_snp import read_snp_by_chr
from parameters import chrom_list
if __name__ == "__main__":
    pool = multiprocessing.Pool(processes=4)
    result = []
    for i in chrom_list:
        msg = 'chr'+str(i)
        print msg
        result.append(pool.apply_async(read_snp_by_chr, (msg, )))
    pool.close()
    pool.join()
    for res in result:
        print res.get()
    print "Sub-process(es) done."