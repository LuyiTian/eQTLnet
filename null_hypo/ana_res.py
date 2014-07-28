import cPickle
from scipy.stats import norm, linregress
from numpy import linspace
tmp_data_dir = "/Users/luyi/tmp_data/"
result = {}
for i in range(1,5):
    tmp_dict = cPickle.load(open(tmp_data_dir+"out_"+str(i)))
    result=dict(result, **tmp_dict)
reg_result = {}
for N in result.keys():
    X = []
    Y_min = []
    Y_max = []
    thr = 1e-10
    for r_val,min_max in result[N].items():
        X.append(r_val)
        Y_min.append(min_max[0])
        Y_max.append(min_max[1])
    slope_lo, intercept_lo, _, p_val, _ = linregress(X,Y_min)
    assert p_val<thr
    slope_hi, intercept_hi, _, p_val, _ = linregress(X,Y_max)
    assert p_val<thr
    reg_result[N] = (slope_lo, intercept_lo, slope_hi, intercept_hi)
cPickle.dump(reg_result,open(tmp_data_dir+"reg_result",'w'))