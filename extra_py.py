import numpy as np
import w2extrapolation

def extra(nu0,nu1,f0,Z0,t,eps=1e-3,deb=False,itmax=10000,tol=5e-4,tau=1,verb=0):

    X,a=nu0[0],nu0[1]
    Y,b=nu1[0],nu1[1]
    n=a.shape[0]
    m=b.shape[0]
    d=X.shape[1]

    f, g, Z, P, niter, errP, errG, obj = w2extrapolation.extra_fortran(n,m,d,a,b,X,Y,f0,Z0,t,eps,deb,itmax,tol,tau,verb) 
    mu=(Z,b)
    dt=1/(1-t)
    Xb=dt*Z+(1-dt)*Y
    nu0bar=(Xb,b)
   
    return f, g, mu, P, nu0bar, niter, errP, errG, obj

def extra_Nto1(num,nup,f0,Z0,t,lm,eps=1e-3,deb=False,itmax=10000,tol=5e-4,tau=1,verb=0):

    Nm=len(num)
    X=[num[k][0] for k in range(Nm)]
    n=[X[k].shape[0] for k in range(Nm)]
    nb=max(n)
    a=tuple([num[k][1] for k in range(Nm)]) #useless step 
    Y, b=nup[0], nup[1]
    m=b.shape[0]
    d=Y.shape[1]
    
    a_all=np.zeros((nb,Nm))
    f0all=np.zeros((nb,Nm))
    Xall=np.zeros((nb,d,Nm))
    for k in range(Nm):
        a_all[0:n[k],k]=a[k]
        Xall[0:n[k],:,k]=X[k]
        f0all[0:n[k],k]=f0[k]
    n=np.array(n)
    lm=np.array(lm)

    fall, g, Z, Pall, niter, errP, errG, obj = w2extrapolation.extra_nto1_fortran(Nm,nb,n,m,d,a_all,b,Xall,Y,f0all,Z0,t,lm,eps,deb,itmax,tol,tau,verb)
    f=tuple( [fall[0:n[k],k] for k in range(Nm)] )
    P=tuple( [Pall[0:n[k],:,k] for k in range(Nm)] )
    mu=(Z,b)
    dt=1/(1-t)
    Xb=dt*Z+(1-dt)*Y
    numbar=(Xb,b)

    return f, g, mu, P, numbar, niter, errP, errG, obj