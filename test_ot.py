import numpy as np
import ot_py
import pickle

file = open('measures_small', 'rb')
images=pickle.load(file)
file.close()
mu=(images['cat'][0], images['cat'][1])
nu=(images['trefle'][0], images['trefle'][1])

# f, g, P, niter, errP, obj =ot_py.ot(mu,nu,eps=1e-3,itmax=1000,tol=5e-4,verb=1)
Sdiv =ot_py.Sdiv(mu,nu,eps=1e-3,itmax=1000,tol=5e-4,verb=1)