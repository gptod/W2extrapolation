import numpy as np
import ot_sinkhorn

def ot(mu,nu,cost=ot_sinkhorn.ot_sinkhorn.cost_l2,eps=1e-3,itmax=10000,tol=5e-4,verb=0):

    X,a=mu[0],mu[1]
    Y,b=nu[0],nu[1]
    n=a.shape[0]
    m=b.shape[0]
    d=X.shape[1]
    C=cost(n,m,d,X,Y)
    
    f0=np.ones(mu[0].shape[0])
    f, g, P, niter, errP, obj =ot_sinkhorn.ot_sinkhorn.ote(n,m,a,b,C,f0,eps,itmax,tol,verb)
    errP=errP[:niter]
    obj=obj[:niter]
    if errP[-1]> tol: print('Warning: error above tolerance. Error:',errP[-1])
    return f, g, P, niter, errP, obj

def ot_self(mu,cost=ot_sinkhorn.ot_sinkhorn.cost_l2,eps=1e-3,itmax=10000,tol=5e-4,verb=0):

    X,a=mu[0],mu[1]
    n=a.shape[0]
    d=X.shape[1]
    C=cost(n,n,d,X,X)
    
    f0=np.ones(n)
    f, P, niter, errP, obj =ot_sinkhorn.ot_sinkhorn.ote_self(n,a,C,f0,eps,itmax,tol,verb)
    if errP> tol: print('Warning: error above tolerance. Error:',errP)

    return f, P, niter, errP, obj

def Sdiv(mu,nu,cost=ot_sinkhorn.ot_sinkhorn.cost_l2,eps=1e-3,itmax=10000,tol=5e-4,verb=0):

    _,_,_,_,_,res12 =ot(mu,nu,cost=cost,eps=eps,itmax=itmax,tol=tol,verb=verb)
    _,_,_,_,res11 =ot_self(mu,cost=cost,eps=eps,itmax=itmax,tol=tol,verb=verb)
    _,_,_,_,res22 =ot_self(nu,cost=cost,eps=eps,itmax=itmax,tol=tol,verb=verb)

    return res12[-1]-0.5*res11-0.5*res22
