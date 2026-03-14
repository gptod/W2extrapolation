import numpy as np
import extra_py as extra
import pickle
import timeit

# T=4 # final time
# Nt=9 # number of steps
T=2 # final time
Nt=3 # number of steps
dt=(T-1)/Nt
d=2 # dimension

eps=1e-3 # temperature value
tau=5 # stepsize
itmax=50000 # maximum number of iterations
tol=5e-4 # stop tolerance

filein='measures_small'
fileout='Nto1_small'
# filein='measures'
# fileout='Nto1'

file = open(filein, 'rb')
images=pickle.load(file)
file.close()

Nm=3
num=( (images['cat'][0],images['cat'][1]), (images['star8'][0],images['star8'][1]), (images['trefle'][0],images['trefle'][1]) )
nup=(images['thinspiral'][0],images['thinspiral'][1])
lm=np.array([1.0/3, 1.0/3, 1.0/3])

extrapolations={}
start = timeit.default_timer()

extrapolations["0-1"]=(num[0],)
extrapolations["0-2"]=(num[1],)
extrapolations["0-3"]=(num[2],)
extrapolations["1"]=(nup,)


for i in range(1,Nt+1):
    t=1+i*dt
    if i==1:
        f0=tuple([np.ones(num[k][0].shape[0]) for k in range(Nm)])
        Z0=nup[0]
    f, g, mu, P, numbar, niter, errP, errG, obj =extra.extra_Nto1(num,nup,f0,Z0,t,lm,eps=eps,deb=True,itmax=itmax,tol=tol,verb=0,bb=True)

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

        




