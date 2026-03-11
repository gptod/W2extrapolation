import numpy as np
import w2extrapolation

def extra(nu0,nu1,f0,Z0,t,eps=1e-3,deb=False,itmax=10000,tol=5e-4,tau=1,verb=0):

    # input:
    # nu0 -> initial measure; tuple containing two numpy arrays, X the particles' positions of size (n,d) and mu the weights of size n
    # nu1 -> final measure; same structure of nu0
    # f0 -> initial value for the potential f; numpy array of size n
    # Z0 -> initial value for the extrapolation's particle positions; numpy array of size (m,d), m corresponding to the number of 
    #       particles of nu1 
    # t -> time of the extrapolation
    # eps -> temperature parameter in the Entropic regularization approach
    # deb -> logical; True (add debiasing term) or False (no debiasing)
    # itmax -> maximum number of iterations for the SISTA algorithm
    # tol -> tolerance for the error
    # tau -> stepsize for the gradient descent step in the SISTA algorithm
    # verb -> verbosity level: 0 no output, 1 print iterations' details

    # output: 
    # f, g -> optimal potentials; numpy arrays of size respextively n and m
    # mu -> extrapolated discrete measure in the form (Z,b), with Z numpy array of size (m,d) corresponding to particles' position 
    #       and b numpy array of size m corresponding to the weights (which coincides with those of the measure nu1)
    # P -> optimal plan in the weak optimal transport barycentric formulation; numpy array of dimension (n,m)
    # nu0bar -> discrete measure in convex order with nu0, with the same form of mu
    # niter -> number of iterations
    # errP, errG -> arrays of iterations' errors on the marginal condition on the plan P and on the norm of the gradient
    # obj -> array of iterations' objective values
    
    X,a=nu0[0],nu0[1]
    Y,b=nu1[0],nu1[1]
    n=a.shape[0]
    m=b.shape[0]
    d=X.shape[1]

    f, g, Z, P, niter, errP, errG, obj = w2extrapolation.extrapolation.extra_fortran(n,m,d,a,b,X,Y,f0,Z0,t,eps,deb,itmax,tol,tau,verb) 
    mu=(Z,b)
    dt=1/(1-t)
    Xb=dt*Z+(1-dt)*Y
    nu0bar=(Xb,b)
   
    if min(errP[-1],errG[-1])> tol: print('Warning: error above tolerance. Error:',min(errP[-1],errG[-1]))

    return f, g, mu, P, nu0bar, niter, errP, errG, obj

def extra_Nto1(num,nup,f0,Z0,t,lm,eps=1e-3,deb=False,itmax=10000,tol=5e-4,tau=1,verb=0):

    # input:
    # num -> initial measures; tuple containing the N measures in the form of tuples (X,a), 
    #        X numpy array of size (nk,d) corresponding to particles' positions and a numpy array of weights of size nk,
    #        nk being the number of particles of the k-th measure in num, d the space dimension
    # nup -> final measure; tuple of the same structure as those in num
    # f0 -> initial value for the potential f; tuple of length len(num) of numpy arrays of size nk
    # Z0 -> initial value for the extrapolation's particle positions; numpy array of size (m,d), m corresponding to the number of 
    #       particles of nup
    # t -> time of the extrapolation
    # eps -> temperature parameter in the Entropic regularization approach
    # deb -> logical; True (add debiasing term) or False (no debiasing)
    # itmax -> maximum number of iterations for the SISTA algorithm
    # tol -> tolerance for the error
    # tau -> stepsize for the gradient descent step in the SISTA algorithm
    # verb -> verbosity level: 0 no output, 1 print iterations' details

    # output: 
    # f, g -> optimal potentials; tuple of length len(num) of numpy arrays of size respextively nk and m
    # mu -> extrapolated discrete measure in the form (Z,b), with Z numpy array of size (m,d) corresponding to particles' position 
    #       and b numpy array of size m corresponding to the weights (which coincides with those of the measure nu1)
    # P -> ptimal plan in the weak optimal transport barycentric formulation; tuple of length len(num) of numpy arrays of dimension (nk,m)
    # nu0bar -> discrete measure in convex order with nu0, with the same form of mu
    # niter -> number of iterations
    # errP, errG -> arrays of iterations' errors on the marginal condition on the plan P and on the norm of the gradient
    # obj -> array of iterations' objective values

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

    fall, g, Z, Pall, niter, errP, errG, obj = w2extrapolation.extrapolation.extra_nto1_fortran(Nm,nb,n,m,d,a_all,b,Xall,Y,f0all,Z0,t,lm,eps,deb,itmax,tol,tau,verb)
    f=tuple( [fall[0:n[k],k] for k in range(Nm)] )
    P=tuple( [Pall[0:n[k],:,k] for k in range(Nm)] )
    mu=(Z,b)
    dt=1/(1-t)
    Xb=dt*Z+(1-dt)*Y
    numbar=(Xb,b)

    if min(errP[-1],errG[-1])> tol: print('Warning: error above tolerance. Error:',min(errP[-1],errG[-1]))
    
    return f, g, mu, P, numbar, niter, errP, errG, obj
