import numpy as np
import w2extrapolation

def extra(nu0,nu1,f0,Z0,t,eps=1e-3,itmax=10000,tol=5e-4,tau=1,verb=0):

    # input:
    # nu0 -> initial measure; tuple containing two numpy arrays, X the particles' positions of size (n,d) and mu the weights of size d
    # nu1 -> final measure; same structure of nu0
    # f0 -> initial value for the potential f; numpy array of size n
    # Z0 -> initial value for the extrapolation's particle positions; numpy array of size (m,d)
    # t -> time of the extrapolation
    # eps -> temperature parameter in the Entropic regularization approach
    # itmax -> maximum number of iterations for the SISTA algorithm
    # tol -> tolerance for the error
    # tau -> stepsize for the gradient descent step in the SISTA algorithm
    # verb -> verbosity level: 0 no output, 1 print iterations' details

    # output: 
    # f, g -> optimal potentials; numpy arrays of size respextively n and m
    # Z -> extrapolation's particle positions
    # P -> optimal plan in the weak optimal transport barycentric formulation
    # niter -> number of iterations
    # errP, errG -> arrays of iterations' errors on the marginal condition on the plan P and on the norm of the gradient
    # obj -> array of iterations' objective values

    # The extrapolation is provided as the couple of Z and nu, where Z is the computed array of particles' position and
    # nu is the array of weights of the final measure nu1
    
    X,mu=nu0[0],nu0[1]
    Y,nu=nu1[0],nu1[1]
    n=mu.shape[0]
    m=nu.shape[0]
    d=X.shape[1]
    return w2extrapolation.extra_fortran(n,m,d,mu,nu,X,Y,f0,Z0,t,eps,itmax,tol,tau,verb) 
