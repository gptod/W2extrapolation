import numpy as np
import matplotlib.pyplot as plt
import pickle


filein='extrapolations_small'
fileout='images_extra_n_small/'
#filein='extrapolations'
#fileout='images_extra_n/'

file = open(filein, 'rb')
extrapolations=pickle.load(file)
file.close()

# colorscale
c0=np.array([59/255,83/255,182/255])
c1=np.array([1,0,0])

# find maximum size of the domain 
a,b,c,d=(0,0,0,0)
for couple in extrapolations:
    for time in extrapolations[couple]:
        X=extrapolations[couple][time][0]
        a=min(a,min(X[:,0]))
        b=max(b,max(X[:,0]))
        c=min(c,min(X[:,1]))
        d=max(d,max(X[:,1]))

for couple in extrapolations:
    
    Nc=len(extrapolations[couple])
    col=np.zeros((Nc,3))
    col[0,:]=c0
    col[1,:]=c0
    for i in range(1,Nc-1):
        col[i+1,:]=(1-i/(Nc-2))*c0+(i/(Nc-2))*c1

    ind=-1
    for time in extrapolations[couple]:
        ind+=1
        X=extrapolations[couple][time][0]
        fig, ax = plt.subplots()
        ax.plot(X[:, 0], X[:, 1], 'or', markerfacecolor=col[ind,:], markersize=5, markeredgecolor=col[ind,:])
        ax.set_xlim(a, b)
        ax.set_ylim(c, d)
        ax.set_aspect('equal')
        ax.axis('off')
        # plt.savefig(fileout+couple+'time'+time+'.jpg')
        plt.savefig(fileout+couple+'time'+time+'.png', transparent=True, bbox_inches='tight', dpi=600)
        plt.close()