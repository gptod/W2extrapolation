import numpy as np
import extra_py as extra
import pickle
import timeit

T=3 # final time
Nt=3 # number of steps
dt=(T-1)/Nt
d=2 # dimension

eps=1e-3 # temperature value
tau=5 # stepsize
itmax=50000 # maximum number of iterations
tol=5e-4 # stop tolerance

filein='measures_small'
fileout='extrapolations_small'
#filein='measures'
#fileout='extrapolations'

file = open(filein, 'rb')
images=pickle.load(file)
file.close()

couples=[['cat','star8'],
         ['cat','trefle'],
         ['cat','thinspiral'],
         ['star8','trefle'],
         ['star8','thinspiral'],
         ['star8','cat'],
         ['thinspiral','cat'],
         ['thinspiral','star8'],
         ['thinspiral','trefle'],
         ['trefle','thinspiral'],
         ['trefle','cat'],
         ['trefle','star8']]

extrapolations={}
start = timeit.default_timer()
for ind in couples:
     
    X=images[ind[0]][0]
    mu=images[ind[0]][1]
    Y=images[ind[1]][0]
    nu=images[ind[1]][1]

    extra_i={'0':(X,mu), '1':(Y,nu)}
    for i in range(1,Nt+1):
        t=1+i*dt
        if i==1:
            f0=np.ones(np.size(X,0))
            Z0=Y
        f, g, Z, P, niter, errP, errG, obj =extra.extra((X,mu),(Y,nu),f0,Z0,t,eps,itmax,tol,tau,0)

        print(ind[0]+'2'+ind[1],'iterations:',niter,'err.marg.:',errP[niter-1],'err.grad.:',errG[niter-1])
        if max(errP[niter-1],errG[niter-1])>tol:
            print('Warning: ',ind[0]+'2'+ind[1],'t=',f"{t:1.4}",'error bigger than tolerance')
        extra_i[f"{t:1.4}"]=(Z,nu,errP[niter-1],errG[niter-1])
        f0=f
        Z0=Z

    extrapolations[ind[0]+'2'+ind[1]]=extra_i

stop = timeit.default_timer()
print('Time: ', stop - start)

file = open(fileout, 'wb')
pickle.dump(extrapolations, file)
file.close()

        




