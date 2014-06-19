import cpp_eqtl
a = cpp_eqtl.PyeQTL_Matrix(400,20000,2000)
for i in xrange(20000):
    for j in xrange(400):
        a.add_exp(float(j),i*400+j)
for i in xrange(2000):
    for j in xrange(400):
        a.add_SNP(float(j),i+j*2000)
a.exp_norm()
a.SNP_norm()
a.matrix_multi()
res = [0.]*(20000*2000)
a.get_r(res)
for i in range(20):
    print res[i]