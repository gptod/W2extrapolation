import numpy as np
import extra_py as extra
import pickle
import timeit

T=4 # final time
Nt=9 # number of steps
dt=(T-1)/Nt
d=2 # dimension

eps=1e-3 # temperature value
tau=5 # stepsize
itmax=50000 # maximum number of iterations
tol=5e-4 # stop tolerance

filein='measures_small'
fileout='trefle2thinspiral_small'
# filein='measures'
# fileout='trefle2thinspiral'

file = open(filein, 'rb')
images=pickle.load(file)
file.close()

extrapolations={}
start = timeit.default_timer()

X=images['trefle'][0]
a=images['trefle'][1]
Y=images['thinspiral'][0]
b=images['thinspiral'][1]

extrapolations={'0':((X,a),), '1':((Y,b),)}
for i in range(1,Nt+1):
    t=1+i*dt
    if i==1:
        f0=np.ones(np.size(X,0))
        Z0=Y
    f, g, mu, P, nu0bar, niter, errP, errG, obj =extra.extra((X,a),(Y,b),f0,Z0,t,eps,False,itmax,tol,tau,0)

    print('iterations:',niter,'err.marg.:',errP[niter-1],'err.grad.:',errG[niter-1])
    if max(errP[niter-1],errG[niter-1])>tol:
        print('Warning: ','t=',f"{t:1.4}",'error bigger than tolerance')
    extrapolations[f"{t:1.4}"]=(mu,errP[niter-1],errG[niter-1])
    f0=f
    Z0=mu[0]


stop = timeit.default_timer()
print('Time: ', stop - start)

file = open(fileout, 'wb')
pickle.dump(extrapolations, file)
file.close()

        




