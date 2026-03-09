import numpy as np
import matplotlib.pyplot as plt
import pickle


filein='trefle2thinspiral_small'
fileout='images_trefle2thinspiral_small/'
# filein='trefle2thinspiral'
# fileout='images_trefle2thinspiral/'

file = open(filein, 'rb')
extrapolations=pickle.load(file)
file.close()

# colorscale
c0=np.array([59/255,83/255,182/255])
c1=np.array([1,0,0])

# find maximum size of the domain 
a,b,c,d=(0,0,0,0)
for time in extrapolations:
    X=extrapolations[time][0][0]
    a=min(a,min(X[:,0]))
    b=max(b,max(X[:,0]))
    c=min(c,min(X[:,1]))
    d=max(d,max(X[:,1]))
    
Nc=len(extrapolations)
col=np.zeros((Nc,3))
col[0,:]=c0
col[1,:]=c0
for i in range(1,Nc-1):
    col[i+1,:]=(1-i/(Nc-2))*c0+(i/(Nc-2))*c1

ind=-1
for time in extrapolations:
    ind+=1
    X=extrapolations[time][0][0]
    fig, ax = plt.subplots()
    ax.plot(X[:, 0], X[:, 1], 'or', markerfacecolor=col[ind,:], markersize=5, markeredgecolor=col[ind,:])
    ax.set_xlim(a, b)
    ax.set_ylim(c, d)
    ax.set_aspect('equal')
    ax.axis('off')
    # plt.savefig(fileout+couple+'time'+time+'.jpg')
    plt.savefig(fileout+'trefle2thinspiral'+'time'+time+'.png', transparent=True, bbox_inches='tight', dpi=600)
    plt.close()