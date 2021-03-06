import rpy2.robjects as robjects
import numpy as np
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
import e_divisive as pye
import e_agglomerative

utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)
packnames = ['ecp', 'MASS']
names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
print(names_to_install)
if len(names_to_install) > 0:
    utils.install_packages(StrVector(names_to_install))

data = robjects.r("""
    library(ecp)
    set.seed(100)
    x1 = matrix(c(rnorm(100),rnorm(100,3),rnorm(100,0,2)))
    x2 = rbind(MASS::mvrnorm(100,c(0,0),diag(2)), MASS::mvrnorm(100,c(2,2),diag(2)))
    y1 = e.divisive(X=x1,sig.lvl=0.05,R=199,k=NULL,min.size=30,alpha=1)
    y2 = e.divisive(X=x2,sig.lvl=0.05,R=499,k=NULL,min.size=30,alpha=1)
    y1_no_perm = e.divisive(X = x1, k = 2, min.size = 30, alpha = 1)
    y2_no_perm = e.divisive(X = x2, k = 1, min.size = 30, alpha = 1)
    
    mem = rep(c(1,2,3,4),times=c(10,10,10,10))
    x_agg = as.matrix(c(rnorm(10,0,1),rnorm(20,2,1),rnorm(10,-1,1)))
    y_agg = e.agglo(X=x_agg,member=mem,alpha=1,penalty=function(cp,Xts) 0)
    """)

py_agglo_X1 = np.array(robjects.r["x_agg"])
py_agglo_mem1 = np.array(robjects.r["mem"])
py_agglo_Y1 = \
    e_agglomerative.e_agglo(X = py_agglo_X1, \
                            member=py_agglo_mem1, \
                            alpha = 1, \
                            penalty = (lambda cp: 0))
print(py_agglo_Y1)
rY1 = robjects.r["y_agg"]
print(rY1)

pyX1 = np.array(robjects.r["x1"])
pyY1_no_perm = pye.e_divisive(X = pyX1, k = 2, min_size = 30, alpha = 1)
print(pyY1_no_perm)
rY1_no_perm = robjects.r["y1_no_perm"]
print(rY1_no_perm)

pyX2 = np.array(robjects.r["x2"])
pyY2_no_perm = pye.e_divisive(X = pyX2, k = 1, min_size = 30, alpha = 1)
print(pyY2_no_perm)
rY2_no_perm = robjects.r["y2_no_perm"]
print(rY2_no_perm)

pyY1 = pye.e_divisive(X = pyX1, sig_lvl = 0.05, R = 199, k = None, min_size = 30, alpha = 1)
print(pyY1)
rY1 = robjects.r["y1"]
print(rY1)

pyY2 = pye.e_divisive(X = pyX2, sig_lvl = 0.05, R = 499, k = None, min_size = 30, alpha = 1)
print(pyY2)
rY2 = robjects.r["y2"]
print(rY2)
