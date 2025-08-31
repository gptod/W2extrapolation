import numpy as np
import w2extrapolation

def extra(nu0,nu1,f0,Z0,t,eps=1e-3,itmax=10000,tol=5e-4,tau=1,verb=0):

    X,mu=nu0[0],nu0[1]
    Y,nu=nu1[0],nu1[1]
    n=mu.shape[0]
    m=nu.shape[0]
    d=X.shape[1]
    return w2extrapolation.extra_fortran(n,m,d,mu,nu,X,Y,f0,Z0,t,eps,itmax,tol,tau,verb) 
