from random import gauss, random, sample
from scipy.stats import linregress
import multiprocessing
import cPickle

def sampling_r(N,lo_N,hi_N, sample_size,res_name):
    '''
    return in file {subset_size:{r_val:(low_boundry,upper_boundry)}, ...}
    '''
    ind_list = range(N)
    res = dict([(i,{}) for i in range(lo_N,hi_N)])
    for coeff in [0.01+0.01*a_i for a_i in range(150)]:
        for n in xrange(lo_N,hi_N):
            X = [gauss(0,1) for i in range(N)]
            Y = [coeff*X[i]+gauss(0,1) for i in range(N)]
            _, _, all_r, _, _ = linregress(X,Y)
            tmp_l = []
            for ith_sam in xrange(sample_size):
                ind_l = sample(ind_list,n)
                _, _, sub_r, _, _ = linregress([x for i,x in enumerate(X) if i in ind_l],\
                    [y for i,y in enumerate(Y) if i in ind_l])
                tmp_l.append(sub_r)
            tmp_l.sort()
            res[n][all_r]=(tmp_l[0],tmp_l[-1])
            print n,all_r,tmp_l[0],tmp_l[-1]
    f = open(res_name,"w")
    cPickle.dump(res,f)
    f.close()
    return "OK"

if __name__ == "__main__":
    N = 420
    min_N = int(N*0.1)
    max_N = int(N*0.9)
    step = int((max_N-min_N)/4)
    sample_size = 1000
    result = []
    pool = multiprocessing.Pool(processes=4)
    for i in range(1,5):
        if i < 4:msg = (N,min_N+(i-1)*step,min_N+i*step, sample_size,"out_"+str(i))
        if i == 4:msg = (N,min_N+(i-1)*step,max_N, sample_size,"out_"+str(i))
        result.append(pool.apply_async(sampling_r, msg))
    pool.close()
    pool.join()
    print "Sub-process(es) done."