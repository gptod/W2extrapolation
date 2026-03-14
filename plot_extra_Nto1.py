import numpy as np
import matplotlib.pyplot as plt
import pickle


filein='Nto1_small'
fileout='images_extra_Nto1_small/'
# filein='Nto1'
# fileout='images_extra_Nto1/'

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
a=a-0.05*(b-a)
b=b+0.05*(b-a)
c=c-0.05*(d-c)
d=d+0.05*(d-c)
    
Nc=len(extrapolations)
col=np.zeros((Nc,3))
col[0,:]=c0
col[1,:]=c0
col[2,:]=c0
col[3,:]=c0
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
    # if ind<=2:
    #    plt.savefig(fileout+'Nto1_'+f"{ind:1}"+'.png', transparent=True, bbox_inches='tight', dpi=600)
    # else:
    plt.savefig(fileout+'Nto1'+'time'+time+'.png', transparent=True, bbox_inches='tight', dpi=600)
    plt.close()